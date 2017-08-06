/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.sam.markduplicates;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;
import picard.sam.DuplicationMetrics;
import picard.sam.markduplicates.util.AbstractOpticalDuplicateFinderCommandLineProgram;
import picard.sam.util.PhysicalLocationShort;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.*;

import static java.lang.Math.pow;

/**
 * <p>Attempts to estimate library complexity from sequence alone. Does so by sorting all reads
 * by the first N bases (5 by default) of each read and then comparing reads with the first
 * N bases identical to each other for duplicates.  Reads are considered to be duplicates if
 * they match each other with no gaps and an overall mismatch rate less than or equal to
 * MAX_DIFF_RATE (0.03 by default).</p>
 * <p/>
 * <p>Reads of poor quality are filtered out so as to provide a more accurate estimate. The filtering
 * removes reads with any no-calls in the first N bases or with a mean base quality lower than
 * MIN_MEAN_QUALITY across either the first or second read.</p>
 * <p/>
 * <p>The algorithm attempts to detect optical duplicates separately from PCR duplicates and excludes
 * these in the calculation of library size. Also, since there is no alignment to screen out technical
 * reads one further filter is applied on the data.  After examining all reads a Histogram is built of
 * [#reads in duplicate set -> #of duplicate sets]; all bins that contain exactly one duplicate set are
 * then removed from the Histogram as outliers before library size is estimated.</p>
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = EstimateLibraryComplexity.USAGE_SUMMARY + EstimateLibraryComplexity.USAGE_DETAILS,
        oneLineSummary = EstimateLibraryComplexity.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class EstimateLibraryComplexity extends AbstractOpticalDuplicateFinderCommandLineProgram {
    static final String USAGE_SUMMARY = "Estimates the numbers of unique molecules in a sequencing library.  ";
    static final String USAGE_DETAILS = "<p>This tool outputs quality metrics for a sequencing library preparation." +
            "Library complexity refers to the number of unique DNA fragments present in a given library.  Reductions in complexity " +
            "resulting from PCR amplification during library preparation will ultimately compromise downstream analyses " +
            "via an elevation in the number of duplicate reads.  PCR-associated duplication artifacts can result from: inadequate amounts " +
            "of starting material (genomic DNA, cDNA, etc.), losses during cleanups, and size selection issues.  " +
            "Duplicate reads can also arise from optical duplicates resulting from sequencing-machine optical sensor artifacts.</p>  " +

            "<p>This tool attempts to estimate library complexity from sequence of read pairs alone.  Reads are sorted by the first N bases " +
            "(5 by default) of the first read and then the first N bases of the second read of a pair.   Read pairs are considered to " +
            "be duplicates if they match each other with no gaps and an overall mismatch rate less than or equal to MAX_DIFF_RATE " +
            "(0.03 by default).  Reads of poor quality are filtered out to provide a more accurate estimate.  The filtering removes reads" +
            " with any poor quality bases as defined by a read's MIN_MEAN_QUALITY (20 is the default value) across either the first or " +
            "second read.  Unpaired reads are ignored in this computation.</p> " +
            "" +
            "<p>The algorithm attempts to detect optical duplicates separately from PCR duplicates and excludes these in the calculation " +
            "of library size.  Also, since there is no alignment information used in this algorithm, an additional filter is applied to " +
            "the data as follows.  After examining all reads, a histogram is built in which the number of reads in a duplicate set is " +
            "compared with the number of of duplicate sets.   All bins that contain exactly one duplicate set are then removed from the " +
            "histogram as outliers prior to the library size estimation.  </p>" +

            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar EstimateLibraryComplexity \\<br />" +
            "     I=input.bam \\<br />" +
            "     O=est_lib_complex_metrics.txt" +
            "</pre>" +
            "Please see the documentation for the companion " +
            "<a href='https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates'>MarkDuplicates</a> tool." +
            "<hr />";
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "One or more files to combine and " +
            "estimate library complexity from. Reads can be mapped or unmapped.")
    public List<File> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file to writes per-library metrics to.")
    public File OUTPUT;

    @Argument(doc = "The minimum number of bases at the starts of reads that must be identical for reads to " +
            "be grouped together for duplicate detection.  In effect total_reads / 4^max_id_bases reads will " +
            "be compared at a time, so lower numbers will produce more accurate results but consume " +
            "exponentially more memory and CPU.")
    public int MIN_IDENTICAL_BASES = 5;

    @Argument(doc = "The maximum rate of differences between two reads to call them identical.")
    public double MAX_DIFF_RATE = 0.03;

    @Argument(doc = "The minimum mean quality of the bases in a read pair for the read to be analyzed. Reads with " +
            "lower average quality are filtered out and not considered in any calculations.")
    public int MIN_MEAN_QUALITY = 20;

    @Argument(doc = "Do not process self-similar groups that are this many times over the mean expected group size. " +
            "I.e. if the input contains 10m read pairs and MIN_IDENTICAL_BASES is set to 5, then the mean expected " +
            "group size would be approximately 10 reads.")
    public int MAX_GROUP_RATIO = 500;

    @Argument(doc = "Barcode SAM tag (ex. BC for 10X Genomics)", optional = true)
    public String BARCODE_TAG = null;

    @Argument(doc = "Read one barcode SAM tag (ex. BX for 10X Genomics)", optional = true)
    public String READ_ONE_BARCODE_TAG = null;

    @Argument(doc = "Read two barcode SAM tag (ex. BX for 10X Genomics)", optional = true)
    public String READ_TWO_BARCODE_TAG = null;

    @Argument(doc = "The maximum number of bases to consider when comparing reads (0 means no maximum).", optional = true)
    public int MAX_READ_LENGTH = 0;

    @Argument(doc = "Minimum number group count.  On a per-library basis, we count the number of groups of duplicates " +
            "that have a particular size.  Omit from consideration any count that is less than this value.  For " +
            "example, if we see only one group of duplicates with size 500, we omit it from the metric calculations if " +
            "MIN_GROUP_COUNT is set to two.  Setting this to two may help remove technical artifacts from the library " +
            "size calculation, for example, adapter dimers.", optional = true)
    public int MIN_GROUP_COUNT = 2;

    private final Log log = Log.getInstance(EstimateLibraryComplexity.class);

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errorMsgs = new ArrayList<String>();
        if (0 < MAX_READ_LENGTH && MAX_READ_LENGTH < MIN_IDENTICAL_BASES) {
            errorMsgs.add("MAX_READ_LENGTH must be greater than MIN_IDENTICAL_BASES");
        }
        if (MIN_IDENTICAL_BASES <= 0) {
            errorMsgs.add("MIN_IDENTICAL_BASES must be greater than 0");
        }
        return errorMsgs.isEmpty() ? super.customCommandLineValidation() : errorMsgs.toArray(new String[errorMsgs.size()]);
    }

    /**
     * Little class to hold the sequence of a pair of reads and tile location information.
     */
    static class PairedReadSequence extends PhysicalLocationShort {

        //just for rough estimate size of reads, does not affect the fundamental operation of the algorithm
        static final int NUMBER_BASES_IN_READ = 150;

        short readGroup = -1;
        boolean qualityOk = true;
        byte[] read1;
        byte[] read2;
        short libraryId;

        // Hashes corresponding to read1 and read2
        int[] hashes1;
        int[] hashes2;

        public static int getSizeInBytes() {
            // rough guess at memory footprint, summary size of all fields
            return 16 + 4 + (2 * 4) + 1 + 2 * (24 + 8 + NUMBER_BASES_IN_READ) + 2 + (2 * (24 + 8)) + 8 + 4;
        }

        public short getReadGroup() { return this.readGroup; }

        public void setReadGroup(final short readGroup) { this.readGroup = readGroup; }

        public short getLibraryId() { return this.libraryId; }

        public void setLibraryId(final short libraryId) { this.libraryId = libraryId; }

        public static SortingCollection.Codec<PairedReadSequence> getCodec() {
            return new PairedReadCodec();
        }

        void initHashes(int numberOfHashes, int skippedBases, int minReadLength) {
            hashes1 = getHashes(read1, numberOfHashes, skippedBases, minReadLength);
            hashes2 = getHashes(read2, numberOfHashes, skippedBases, minReadLength);
        }

        // Split read by numberOfHashes parts and hash each part
        // For instance:
        //        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
        // read = A C G T A C G T A C G  T  A  C  G  T  A  C  G  T
        // numberOfHashes = 5
        // skippedBases = 1
        // minReadLength = 15
        // So, method returns hashValues with 5 hash value
        // first value calculated from read[1], read[6], read[11]
        // second value calculated from read[2], read[7], read[12]
        // etc.
        // chars from 16 to 19 position is a tail, see compareTails() in ElcHashBasedDuplicatesFinder
        private int[] getHashes(byte[] read, int numberOfHashes, int skippedBases, int minReadLength) {
            final int[] hashValues = new int[numberOfHashes];
            for (int i = 0; i < numberOfHashes; ++i) {
                hashValues[i] = 1;
                int position = skippedBases + i;
                while (position < minReadLength) {
                    hashValues[i] = 31 * hashValues[i] + read[position];
                    position += numberOfHashes;
                }
            }
            return hashValues;
        }
    }

    static class PairedReadSequenceWithBarcodes extends PairedReadSequence {
        int barcode; // primary barcode for this read (and pair)
        int readOneBarcode; // read one barcode, 0 if not present
        int readTwoBarcode; // read two barcode, 0 if not present or not paired

        public PairedReadSequenceWithBarcodes() {
            super();
        }

        public PairedReadSequenceWithBarcodes(final PairedReadSequence val) {
            if (null == val) throw new PicardException("val was null");
            this.readGroup = val.getReadGroup();
            this.tile = val.getTile();
            this.x = val.getX();
            this.y = val.getY();
            this.qualityOk = val.qualityOk;
            this.read1 = val.read1.clone();
            this.read2 = val.read2.clone();
            this.libraryId = val.getLibraryId();
        }

        public static int getSizeInBytes() {
            return PairedReadSequence.getSizeInBytes() + (3 * 4); // rough guess at memory footprint
        }
    }

    /**
     * Codec class for writing and read PairedReadSequence objects.
     */
    static class PairedReadCodec implements SortingCollection.Codec<PairedReadSequence> {
        protected DataOutputStream out;
        protected DataInputStream in;

        public void setOutputStream(final OutputStream out) {
            this.out = new DataOutputStream(out);
        }

        public void setInputStream(final InputStream in) {
            this.in = new DataInputStream(in);
        }

        public void encode(final PairedReadSequence val) {
            try {
                this.out.writeShort(val.readGroup);
                this.out.writeShort(val.tile);
                this.out.writeShort(val.x);
                this.out.writeShort(val.y);
                this.out.writeInt(val.read1.length);
                this.out.write(val.read1);
                this.out.writeInt(val.read2.length);
                this.out.write(val.read2);
            } catch (final IOException ioe) {
                throw new PicardException("Error write out read pair.", ioe);
            }
        }

        public PairedReadSequence decode() {
            try {
                final PairedReadSequence val = new PairedReadSequence();
                try {
                    val.readGroup = this.in.readShort();
                } catch (final EOFException eof) {
                    return null;
                }

                val.tile = this.in.readShort();
                val.x = this.in.readShort();
                val.y = this.in.readShort();

                int length = this.in.readInt();
                val.read1 = new byte[length];
                if (this.in.read(val.read1) != length) {
                    throw new PicardException("Could not read " + length + " bytes from temporary file.");
                }

                length = this.in.readInt();
                val.read2 = new byte[length];
                if (this.in.read(val.read2) != length) {
                    throw new PicardException("Could not read " + length + " bytes from temporary file.");
                }

                return val;
            } catch (final IOException ioe) {
                throw new PicardException("Exception reading read pair.", ioe);
            }
        }

        @Override
        public SortingCollection.Codec<PairedReadSequence> clone() { return new PairedReadCodec(); }
    }


    /**
     * Codec class for writing and read PairedReadSequence objects.
     */
    static class PairedReadWithBarcodesCodec extends PairedReadCodec {
        @Override
        public void encode(final PairedReadSequence val) {
            if (!(val instanceof PairedReadSequenceWithBarcodes)) {
                throw new PicardException("Val was not a PairedReadSequenceWithBarcodes");
            }
            final PairedReadSequenceWithBarcodes data = (PairedReadSequenceWithBarcodes) val;

            super.encode(val);

            try {
                this.out.writeInt(data.barcode);
                this.out.writeInt(data.readOneBarcode);
                this.out.writeInt(data.readTwoBarcode);
            } catch (final IOException ioe) {
                throw new PicardException("Error write out read pair.", ioe);
            }
        }

        @Override
        public PairedReadSequence decode() {
            try {
                final PairedReadSequence parentVal = super.decode();
                if (null == parentVal) return null; // EOF
                final PairedReadSequenceWithBarcodes val = new PairedReadSequenceWithBarcodes(parentVal);
                val.barcode = this.in.readInt();
                val.readOneBarcode = this.in.readInt();
                val.readTwoBarcode = this.in.readInt();

                return val;
            } catch (final IOException ioe) {
                throw new PicardException("Exception reading read pair.", ioe);
            }
        }

        @Override
        public SortingCollection.Codec<PairedReadSequence> clone() { return new PairedReadWithBarcodesCodec(); }
    }

    /**
     * Comparator that orders read pairs on the first N bases of both reads.
     * There is no tie-breaking, so any sort is stable, not total.
     */
    private class PairedReadComparator implements Comparator<PairedReadSequence> {
        final int BASES = EstimateLibraryComplexity.this.MIN_IDENTICAL_BASES;

        public int compare(final PairedReadSequence lhs, final PairedReadSequence rhs) {
            // First compare the first N bases of the first read
            for (int i = 0; i < BASES; ++i) {
                final int retval = lhs.read1[i] - rhs.read1[i];
                if (retval != 0) return retval;
            }

            // Then compare the first N bases of the second read
            for (int i = 0; i < BASES; ++i) {
                final int retval = lhs.read2[i] - rhs.read2[i];
                if (retval != 0) return retval;
            }

            return 0;
        }
    }

    public int getBarcodeValue(final SAMRecord record) {
        return getReadBarcodeValue(record, BARCODE_TAG);
    }

    public static int getReadBarcodeValue(final SAMRecord record, final String tag) {
        if (null == tag) return 0;
        final String attr = record.getStringAttribute(tag);
        if (null == attr) return 0;
        else return attr.hashCode();
    }

    private int getReadOneBarcodeValue(final SAMRecord record) {
        return getReadBarcodeValue(record, READ_ONE_BARCODE_TAG);
    }

    private int getReadTwoBarcodeValue(final SAMRecord record) {
        return getReadBarcodeValue(record, READ_TWO_BARCODE_TAG);
    }

    /** Stock main method. */
    public static void main(final String[] args) {
        new EstimateLibraryComplexity().instanceMainWithExit(args);
    }

    public EstimateLibraryComplexity() {
        final int sizeInBytes;
        if (null != BARCODE_TAG || null != READ_ONE_BARCODE_TAG || null != READ_TWO_BARCODE_TAG) {
            sizeInBytes = PairedReadSequenceWithBarcodes.getSizeInBytes();
        } else {
            sizeInBytes = PairedReadSequence.getSizeInBytes();
        }
        MAX_RECORDS_IN_RAM = (int) (Runtime.getRuntime().maxMemory() / sizeInBytes) / 2;
    }

    /**
     * Method that does most of the work.  Reads through the input BAM file and extracts the
     * read sequences of each read pair and sorts them via a SortingCollection.  Then traverses
     * the sorted reads and looks at small groups at a time to find duplicates.
     */
    @Override
    protected int doWork() {
        for (final File f : INPUT) IOUtil.assertFileIsReadable(f);

        log.info("Will store " + MAX_RECORDS_IN_RAM + " read pairs in memory before sorting.");

        final List<SAMReadGroupRecord> readGroups = new ArrayList<SAMReadGroupRecord>();
        final SortingCollection<PairedReadSequence> sorter;
        final boolean useBarcodes = (null != BARCODE_TAG || null != READ_ONE_BARCODE_TAG || null != READ_TWO_BARCODE_TAG);

        if (!useBarcodes) {
            sorter = SortingCollection.newInstance(PairedReadSequence.class,
                    new PairedReadCodec(),
                    new PairedReadComparator(),
                    MAX_RECORDS_IN_RAM,
                    TMP_DIR);
        } else {
            sorter = SortingCollection.newInstance(PairedReadSequence.class,
                    new PairedReadWithBarcodesCodec(),
                    new PairedReadComparator(),
                    MAX_RECORDS_IN_RAM,
                    TMP_DIR);
        }

        // Loop through the input files and pick out the read sequences etc.
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Read");
        for (final File f : INPUT) {
            final Map<String, PairedReadSequence> pendingByName = new HashMap<String, PairedReadSequence>();
            final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(f);
            readGroups.addAll(in.getFileHeader().getReadGroups());

            for (final SAMRecord rec : in) {
                if (!rec.getReadPairedFlag()) continue;
                if (!rec.getFirstOfPairFlag() && !rec.getSecondOfPairFlag()) {
                    continue;
                }
                if (rec.isSecondaryOrSupplementary()) continue;

                PairedReadSequence prs = pendingByName.remove(rec.getReadName());
                if (prs == null) {
                    // Make a new paired read object and add RG and physical location information to it
                    prs = useBarcodes ? new PairedReadSequenceWithBarcodes() : new PairedReadSequence();
                    if (opticalDuplicateFinder.addLocationInformation(rec.getReadName(), prs)) {
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

                final PairedReadSequenceWithBarcodes prsWithBarcodes = (useBarcodes) ? (PairedReadSequenceWithBarcodes) prs : null;

                if (rec.getFirstOfPairFlag()) {
                    prs.read1 = bases;
                    if (useBarcodes) {
                        prsWithBarcodes.barcode = getBarcodeValue(rec);
                        prsWithBarcodes.readOneBarcode = getReadOneBarcodeValue(rec);
                    }
                } else {
                    prs.read2 = bases;
                    if (useBarcodes) {
                        prsWithBarcodes.readTwoBarcode = getReadTwoBarcodeValue(rec);
                    }
                }

                if (prs.read1 != null && prs.read2 != null && prs.qualityOk) {
                    sorter.add(prs);
                }

                progress.record(rec);
            }
            CloserUtil.close(in);
        }

        log.info(String.format("Finished reading - read %d records - moving on to scanning for duplicates.", progress.getCount()));

        // Now go through the sorted reads and attempt to find duplicates
        final PeekableIterator<PairedReadSequence> iterator = new PeekableIterator<PairedReadSequence>(sorter.iterator());

        final Map<String, Histogram<Integer>> duplicationHistosByLibrary = new HashMap<String, Histogram<Integer>>();
        final Map<String, Histogram<Integer>> opticalHistosByLibrary = new HashMap<String, Histogram<Integer>>();

        int groupsProcessed = 0;
        long lastLogTime = System.currentTimeMillis();
        final int meanGroupSize = (int) (Math.max(1, (progress.getCount() / 2) / (int) pow(4, MIN_IDENTICAL_BASES * 2)));

        ElcDuplicatesFinderResolver algorithmResolver = new ElcDuplicatesFinderResolver(
                MAX_DIFF_RATE,
                MAX_READ_LENGTH,
                MIN_IDENTICAL_BASES,
                useBarcodes,
                opticalDuplicateFinder
        );

        while (iterator.hasNext()) {
            // Get the next group and split it apart by library
            final List<PairedReadSequence> group = getNextGroup(iterator);

            if (group.size() > meanGroupSize * MAX_GROUP_RATIO) {
                final PairedReadSequence prs = group.get(0);
                log.warn("Omitting group with over " + MAX_GROUP_RATIO + " times the expected mean number of read pairs. " +
                        "Mean=" + meanGroupSize + ", Actual=" + group.size() + ". Prefixes: " +
                        StringUtil.bytesToString(prs.read1, 0, MIN_IDENTICAL_BASES) +
                        " / " +
                        StringUtil.bytesToString(prs.read2, 0, MIN_IDENTICAL_BASES));
            } else {
                final Map<String, List<PairedReadSequence>> sequencesByLibrary = splitByLibrary(group, readGroups);

                // Now process the reads by library
                for (final Map.Entry<String, List<PairedReadSequence>> entry : sequencesByLibrary.entrySet()) {
                    final String library = entry.getKey();
                    final List<PairedReadSequence> seqs = entry.getValue();

                    Histogram<Integer> duplicationHisto = duplicationHistosByLibrary.get(library);
                    Histogram<Integer> opticalHisto = opticalHistosByLibrary.get(library);
                    if (duplicationHisto == null) {
                        duplicationHisto = new Histogram<>("duplication_group_count", library);
                        opticalHisto = new Histogram<>("duplication_group_count", "optical_duplicates");
                        duplicationHistosByLibrary.put(library, duplicationHisto);
                        opticalHistosByLibrary.put(library, opticalHisto);
                    }

                    algorithmResolver.resolveAndSearch(seqs, duplicationHisto, opticalHisto);
                }

                ++groupsProcessed;
                if (lastLogTime < System.currentTimeMillis() - 60000) {
                    log.info("Processed " + groupsProcessed + " groups.");
                    lastLogTime = System.currentTimeMillis();
                }
            }
        }

        iterator.close();
        sorter.cleanup();

        final MetricsFile<DuplicationMetrics, Integer> file = getMetricsFile();
        for (final String library : duplicationHistosByLibrary.keySet()) {
            final Histogram<Integer> duplicationHisto = duplicationHistosByLibrary.get(library);
            final Histogram<Integer> opticalHisto = opticalHistosByLibrary.get(library);
            final DuplicationMetrics metrics = new DuplicationMetrics();
            metrics.LIBRARY = library;

            // Filter out any bins that have fewer than MIN_GROUP_COUNT entries in them and calculate derived metrics
            for (final Integer bin : duplicationHisto.keySet()) {
                final double duplicateGroups = duplicationHisto.get(bin).getValue();
                final double opticalDuplicates = opticalHisto.get(bin) == null ? 0 : opticalHisto.get(bin).getValue();

                if (duplicateGroups >= MIN_GROUP_COUNT) {
                    metrics.READ_PAIRS_EXAMINED += (bin * duplicateGroups);
                    metrics.READ_PAIR_DUPLICATES += ((bin - 1) * duplicateGroups);
                    metrics.READ_PAIR_OPTICAL_DUPLICATES += opticalDuplicates;
                }
            }

            metrics.calculateDerivedFields();
            file.addMetric(metrics);
            file.addHistogram(duplicationHisto);

        }

        file.write(OUTPUT);

        return 0;
    }

    /**
     * Pulls out of the iterator the next group of reads that can be compared to each other to
     * identify duplicates.
     */
    List<PairedReadSequence> getNextGroup(final PeekableIterator<PairedReadSequence> iterator) {
        final List<PairedReadSequence> group = new ArrayList<PairedReadSequence>();
        final PairedReadSequence first = iterator.next();
        group.add(first);

        outer:
        while (iterator.hasNext()) {
            final PairedReadSequence next = iterator.peek();
            for (int i = 0; i < MIN_IDENTICAL_BASES; ++i) {
                if (first.read1[i] != next.read1[i] || first.read2[i] != next.read2[i]) break outer;
            }
            group.add(iterator.next());
        }
        return group;
    }

    /**
     * Takes a list of PairedReadSequence objects and splits them into lists by library.
     */
    Map<String, List<PairedReadSequence>> splitByLibrary(final List<PairedReadSequence> input,
                                                         final List<SAMReadGroupRecord> rgs) {

        final Map<String, List<PairedReadSequence>> out = new HashMap<>();
        for (final PairedReadSequence seq : input) {
            String library;
            if (seq.getReadGroup() != -1) {
                library = rgs.get(seq.getReadGroup()).getLibrary();
                if (library == null) library = "Unknown";
            } else {
                library = "Unknown";
            }

            List<PairedReadSequence> librarySeqs = out.get(library);
            if (librarySeqs == null) {
                librarySeqs = new ArrayList<>();
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

        for (int i = 0; i < seedLength; ++i) {
            if (SequenceUtil.isNoCall(bases[i])) return false;
        }

        final int maxReadLength = (MAX_READ_LENGTH <= 0) ? Integer.MAX_VALUE : MAX_READ_LENGTH;
        final int readLength = Math.min(bases.length, maxReadLength);
        int total = 0;
        for (int i = 0; i < readLength; i++) total += quals[i];
        return total / readLength >= minQuality;
    }
}