package net.sf.picard.sam;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.util.Histogram;
import net.sf.picard.util.Log;
import net.sf.picard.util.PeekableIterator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.SequenceUtil;
import net.sf.samtools.util.SortingCollection;

import java.io.*;
import java.util.*;

/**
 * <p>Attempts to estimate library complexity from sequence alone. Does so by sorting all reads
 * by the first N bases (5 by default) of each read and then comparing reads with the first
 * N bases identical to each other for duplicates.  Reads are considered to be duplicates if
 * they match each other with no gaps and an overall mismatch rate less than or equal to
 * MAX_DIFF_RATE (0.03 by default).</p>
 *
 * <p>Reads of poor quality are filtered out so as to provide a more accurate estimate. The filtering
 * removes reads with any no-calls in the first N bases or with a mean base quality lower than
 * MIN_MEAN_QUALITY across either the first or second read.</p>
 *
 * <p>The algorithm attempts to detect optical duplicates separately from PCR duplicates and excludes
 * these in the calculation of library size. Also, since there is no alignment to screen out technical
 * reads one further filter is applied on the data.  After examining all reads a histogram is built of
 * [#reads in duplicate set -> #of duplicate sets]; all bins that contain exactly one duplicate set are
 * then removed from the histogram as outliers before library size is estimated.</p>
 *
 * @author Tim Fennell
 */
public class EstimateLibraryComplexity extends AbstractDuplicateFindingAlgorithm {
    @Usage public final String USAGE =
            "Attempts to estimate library complexity from sequence of read pairs alone. Does so by sorting all reads " +
            "by the first N bases (5 by default) of each read and then comparing reads with the first " +
            "N bases identical to each other for duplicates.  Reads are considered to be duplicates if " +
            "they match each other with no gaps and an overall mismatch rate less than or equal to " +
            "MAX_DIFF_RATE (0.03 by default).\n\n" +
            "Reads of poor quality are filtered out so as to provide a more accurate estimate. The filtering " +
            "removes reads with any no-calls in the first N bases or with a mean base quality lower than " +
            "MIN_MEAN_QUALITY across either the first or second read.\n\n" +
            "Unpaired reads are ignored in this computation.\n\n" +
            "The algorithm attempts to detect optical duplicates separately from PCR duplicates and excludes " +
            "these in the calculation of library size. Also, since there is no alignment to screen out technical " +
            "reads one further filter is applied on the data.  After examining all reads a histogram is built of " +
            "[#reads in duplicate set -> #of duplicate sets]; all bins that contain exactly one duplicate set are " +
            "then removed from the histogram as outliers before library size is estimated.";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="One or more files to combine and " +
            "estimate library complexity from. Reads can be mapped or unmapped.")
    public List<File> INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc="Output file to writes per-library metrics to.")
    public File OUTPUT;

    @Option(doc="The minimum number of bases at the starts of reads that must be identical for reads to " +
            "be grouped together for duplicate detection.  In effect total_reads / 4^max_id_bases reads will " +
            "be compared at a time, so lower numbers will produce more accurate results but consume " +
            "exponentially more memory and CPU.")
    public int MIN_IDENTICAL_BASES = 5;

    @Option(doc="The maximum rate of differences between two reads to call them identical.")
    public double MAX_DIFF_RATE = 0.03;

    @Option(doc="The minimum mean quality of the bases in a read pair for the read to be analyzed. Reads with " +
            "lower average quality are filtered out and not considered in any calculations.")
    public int MIN_MEAN_QUALITY = 20;

    private final Log log = Log.getInstance(EstimateLibraryComplexity.class);

    /**
     * Little class to hold the sequence of a pair of reads and tile location information.
     */
    static class PairedReadSequence implements PhysicalLocation {
        static int size_in_bytes = 2 + 1 + 4 + 1 + 300; // rough guess at memory footprint
        short readGroup = -1;
        byte tile = -1;
        short x = -1, y = -1;
        boolean qualityOk = true;
        byte[] read1;
        byte[] read2;

        public short getReadGroup() { return this.readGroup; }
        public void setReadGroup(final short readGroup) { this.readGroup = readGroup; }

        public byte getTile() { return this.tile; }
        public void setTile(final byte tile) { this.tile = tile; }

        public short getX() { return this.x; }
        public void setX(final short x) { this.x = x; }

        public short getY() { return this.y; }
        public void setY(final short y) { this.y = y; }
    }

    /**
     * Codec class for writing and read PairedReadSequence objects.
     */
    static class PairedReadCodec implements SortingCollection.Codec<PairedReadSequence> {
        private DataOutputStream out;
        private DataInputStream in;

        public void setOutputStream(final OutputStream out) {
            this.out = new DataOutputStream(out);
        }

        public void setInputStream(final InputStream in) {
            this.in= new DataInputStream(in);
        }

        public void encode(final PairedReadSequence val) {
            try {
                this.out.writeShort(val.readGroup);
                this.out.writeByte(val.tile);
                this.out.writeShort(val.x);
                this.out.writeShort(val.y);
                this.out.writeInt(val.read1.length);
                this.out.write(val.read1);
                this.out.writeInt(val.read2.length);
                this.out.write(val.read2);
            }
            catch (IOException ioe) {
                throw new PicardException("Error write out read pair.", ioe);
            }
        }

        public PairedReadSequence decode() {
            try {
                final PairedReadSequence val = new PairedReadSequence();
                try {
                    val.readGroup = this.in.readShort();
                }
                catch (EOFException eof) {
                    return null;
                }

                val.tile = this.in.readByte();
                val.x = this.in.readShort();
                val.y = this.in.readShort();

                int length = this.in.readInt();
                val.read1  = new byte[length];
                if (this.in.read(val.read1) != length) {
                    throw new PicardException("Could not read " + length + " bytes from temporary file.");
                }

                length = this.in.readInt();
                val.read2  = new byte[length];
                if (this.in.read(val.read2) != length) {
                    throw new PicardException("Could not read " + length + " bytes from temporary file.");
                }

                return val;
            }
            catch (IOException ioe) {
                throw new PicardException("Exception reading read pair.", ioe);
            }
        }

        @Override
        public SortingCollection.Codec<PairedReadSequence> clone() { return new PairedReadCodec(); }
    }

    /**
     * Comparator that orders read pairs on the first N bases of both reads.
     */
    class PairedReadComparator implements Comparator<PairedReadSequence> {
        final int BASES = EstimateLibraryComplexity.this.MIN_IDENTICAL_BASES;

        public int compare(final PairedReadSequence lhs, final PairedReadSequence rhs) {
            // First compare the first N bases of the first read
            for (int i=0; i<BASES; ++i) {
                final int retval = lhs.read1[i] - rhs.read1[i];
                if (retval != 0) return retval;
            }

            // Then compare the first N bases of the second read
            for (int i=0; i<BASES; ++i) {
                final int retval = lhs.read2[i] - rhs.read2[i];
                if (retval != 0) return retval;
            }

            return System.identityHashCode(lhs) - System.identityHashCode(rhs);
        }
    }

    /** Stock main method. */
    public static void main(final String[] args) {
        new EstimateLibraryComplexity().instanceMainWithExit(args);
    }

    /**
     * Method that does most of the work.  Reads through the input BAM file and extracts the
     * read sequences of each read pair and sorts them via a SortingCollection.  Then traverses
     * the sorted reads and looks at small groups at a time to find duplicates.
     */
    @Override protected int doWork() {
        for (final File f : INPUT) IoUtil.assertFileIsReadable(f);

        final int maxInMemory = (int) (Runtime.getRuntime().maxMemory() / PairedReadSequence.size_in_bytes) / 2;
        log.info("Will store " + maxInMemory + " read pairs in memory before sorting.");

        final List<SAMReadGroupRecord> readGroups = new ArrayList<SAMReadGroupRecord>();
        int recordsRead = 0;
        final SortingCollection<PairedReadSequence> sorter = SortingCollection.newInstance(PairedReadSequence.class,
                                                                                           new PairedReadCodec(),
                                                                                           new PairedReadComparator(),
                                                                                           maxInMemory,
                                                                                           TMP_DIR);

        // Loop through the input files and pick out the read sequences etc.
        for (final File f : INPUT) {
            final Map<String,PairedReadSequence> pendingByName = new HashMap<String, PairedReadSequence>();
            final SAMFileReader in = new SAMFileReader(f);
            readGroups.addAll(in.getFileHeader().getReadGroups());

            for (final SAMRecord rec : in) {
                if (!rec.getReadPairedFlag()) continue;
                if (!rec.getFirstOfPairFlag() && !rec.getSecondOfPairFlag()) {
                    continue;
                }

                PairedReadSequence prs = pendingByName.remove(rec.getReadName());
                if (prs == null) {
                    // Make a new paired read object and add RG and physical location information to it
                    prs = new PairedReadSequence();
                    if (addLocationInformation(rec.getReadName(), prs)) {
                        final SAMReadGroupRecord rg = rec.getReadGroup();
                        if (rg != null) prs.setReadGroup((short) readGroups.indexOf(rg));
                    }

                    pendingByName.put(rec.getReadName(), prs);
                }

                // Read passes quality check if both ends meet the mean quality criteria
                final boolean passesQualityCheck = passesQualityCheck(rec.getReadBases(),
                                                                      rec.getBaseQualities(),
                                                                      MIN_IDENTICAL_BASES,
                                                                      MIN_MEAN_QUALITY);
                prs.qualityOk = prs.qualityOk && passesQualityCheck;

                // Get the bases and restore them to their original orientation if necessary
                final byte[] bases = rec.getReadBases();
                if (rec.getReadNegativeStrandFlag()) SequenceUtil.reverseComplement(bases);

                if (rec.getFirstOfPairFlag()) {
                    prs.read1 = bases;
                }
                else {
                    prs.read2 = bases;
                }

                if (prs.read1 != null && prs.read2 != null && prs.qualityOk) {
                    sorter.add(prs);
                }

                if (++recordsRead % 1000000 == 0) log.info("Read " + recordsRead + " records.");
            }
        }

        log.info("Finished reading - moving on to scanning for duplicates.");

        // Now go through the sorted reads and attempt to find duplicates
        final PeekableIterator<PairedReadSequence> iterator = new PeekableIterator<PairedReadSequence>(sorter.iterator());

        final Map<String, Histogram<Integer>> duplicationHistosByLibrary = new HashMap<String,Histogram<Integer>>();
        final Map<String,Histogram<Integer>> opticalHistosByLibrary = new HashMap<String,Histogram<Integer>>();

        int groupsProcessed = 0;
        final int groupChunkSize = Math.min((int) (Math.pow(4, MIN_IDENTICAL_BASES*2) / 25), 1000000);

        while (iterator.hasNext()) {
            // Get the next group and split it apart by library
            final List<PairedReadSequence> group = getNextGroup(iterator);
            final Map<String,List<PairedReadSequence>> sequencesByLibrary = splitByLibrary(group, readGroups);

            // Now process the reads by library
            for (final Map.Entry<String,List<PairedReadSequence>> entry : sequencesByLibrary.entrySet()) {
                final String library = entry.getKey();
                final List<PairedReadSequence> seqs = entry.getValue();

                Histogram<Integer> duplicationHisto = duplicationHistosByLibrary.get(library);
                Histogram<Integer> opticalHisto     = opticalHistosByLibrary.get(library);
                if (duplicationHisto == null) {
                    duplicationHisto = new Histogram<Integer>("duplication_group_count", library);
                    opticalHisto     = new Histogram<Integer>("duplication_group_count", "optical_duplicates");
                    duplicationHistosByLibrary.put(library, duplicationHisto);
                    opticalHistosByLibrary.put(library, opticalHisto);
                }

                // Figure out if any reads within this group are duplicates of one another
                for (int i=0; i<seqs.size(); ++i) {
                    final PairedReadSequence lhs = seqs.get(i);
                    if (lhs == null) continue;
                    final List<PairedReadSequence> dupes = new ArrayList<PairedReadSequence>();

                    for (int j=i+1; j<seqs.size(); ++j) {
                        final PairedReadSequence rhs = seqs.get(j);
                        if (rhs == null) continue;

                        if (matches(lhs, rhs, MAX_DIFF_RATE)) {
                            dupes.add(rhs);
                            seqs.set(j, null);
                        }
                    }

                    if (dupes.size() > 0) {
                        dupes.add(lhs);
                        final int duplicateCount = dupes.size();
                        duplicationHisto.increment(duplicateCount);

                        final boolean[] flags = findOpticalDuplicates(dupes, OPTICAL_DUPLICATE_PIXEL_DISTANCE);
                        for (final boolean b : flags) {
                            if (b) opticalHisto.increment(duplicateCount);
                        }
                    }
                    else {
                        duplicationHisto.increment(1);
                    }
                }
            }

            if (++groupsProcessed % groupChunkSize == 0) {
                log.info("Processed " + groupsProcessed + " groups.");
            }
        }

        iterator.close();
        sorter.cleanup();

        final MetricsFile<DuplicationMetrics,Integer> file = getMetricsFile();
        for (final String library : duplicationHistosByLibrary.keySet()) {
            final Histogram<Integer> duplicationHisto = duplicationHistosByLibrary.get(library);
            final Histogram<Integer> opticalHisto     = opticalHistosByLibrary.get(library);
            final DuplicationMetrics metrics = new DuplicationMetrics();
            metrics.LIBRARY = library;

            // Filter out any bins that have only a single entry in them and calcu
            for (final Integer bin : duplicationHisto.keySet()) {
                final double duplicateGroups   = duplicationHisto.get(bin).getValue();
                final double opticalDuplicates = opticalHisto.get(bin) == null ? 0 : opticalHisto.get(bin).getValue();

                if (duplicateGroups > 1) {
                    metrics.READ_PAIRS_EXAMINED          += (bin * duplicateGroups);
                    metrics.READ_PAIR_DUPLICATES         += ((bin-1) * duplicateGroups);
                    metrics.READ_PAIR_OPTICAL_DUPLICATES += opticalDuplicates;
                }
            }

            metrics.calculateDerivedMetrics();
            file.addMetric(metrics);
            file.addHistogram(duplicationHisto);

        }

        file.write(OUTPUT);

        return 0;
    }

    /**
     * Checks to see if two reads pairs have sequence that are the same, give or take a few
     * errors/diffs as dictated by the maxDiffRate.
     */
    private boolean matches(final PairedReadSequence lhs, final PairedReadSequence rhs, final double maxDiffRate) {
        final int read1Length = Math.min(lhs.read1.length, rhs.read1.length);
        final int read2Length = Math.min(lhs.read2.length, rhs.read2.length);
        final int maxErrors = (int) Math.floor((read1Length + read2Length) * maxDiffRate);
        int errors = 0;

        // The loop can start from MIN_IDENTICAL_BASES because we've already confirmed that
        // at least those first few bases are identical when sorting.
        for (int i=MIN_IDENTICAL_BASES; i<read1Length; ++i) {
            if (lhs.read1[i] != rhs.read1[i]) {
                if (++errors > maxErrors) return false;
            }
        }

        for (int i=MIN_IDENTICAL_BASES; i<read2Length; ++i) {
            if (lhs.read2[i] != rhs.read2[i]) {
                if (++errors > maxErrors) return false;
            }
        }

        return true;
    }

    /**
     * Pulls out of the iterator the next group of reads that can be compared to each other to
     * identify duplicates.
     */
    List<PairedReadSequence> getNextGroup(final PeekableIterator<PairedReadSequence> iterator) {
        final List<PairedReadSequence> group = new ArrayList<PairedReadSequence>();
        final PairedReadSequence first = iterator.next();
        group.add(first);

        outer: while(iterator.hasNext()) {
            final PairedReadSequence next = iterator.peek();
            for (int i=0; i<MIN_IDENTICAL_BASES; ++i) {
                if (first.read1[i] != next.read1[i] || first.read2[i] != next.read2[i]) break outer;
            }

            group.add(iterator.next());

        }

        return group;
    }

    /**
     * Takes a list of PairedReadSequence objects and splits them into lists by library.
     */
    Map<String,List<PairedReadSequence>> splitByLibrary(final List<PairedReadSequence> input,
                                                        final List<SAMReadGroupRecord> rgs) {

        final Map<String,List<PairedReadSequence>> out = new HashMap<String,List<PairedReadSequence>>();
        for (final PairedReadSequence seq : input) {
            String library = null;
            if (seq.getReadGroup() != -1) {
                library = rgs.get(seq.getReadGroup()).getLibrary();
                if (library == null) library = "Unknown";
            }
            else {
                library = "Unknown";
            }

            List<PairedReadSequence> librarySeqs = out.get(library);
            if (librarySeqs == null) {
                librarySeqs = new ArrayList<PairedReadSequence>();
                out.put(library, librarySeqs);
            }
            librarySeqs.add(seq);
        }

        return out;
    }

    /**
     * Checks that the average quality over the entire read is >= min, and that the first N bases do
     * not contain any no-calls.
     */
    boolean passesQualityCheck(final byte[] bases, final byte[] quals, final int seedLength, final int minQuality) {
        if (bases.length < seedLength) return false;

        for (int i=0; i<seedLength; ++i) {
            if (SequenceUtil.isNoCall(bases[i])) return false;
        }

        int total = 0;
        for (final byte b : quals) total += b;
        return total / quals.length >= minQuality;
    }
}
