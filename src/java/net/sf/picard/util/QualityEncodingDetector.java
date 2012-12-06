package net.sf.picard.util;

import net.sf.picard.PicardException;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import java.util.*;

import static java.util.Arrays.asList;

/**
 * Utility for determining the type of quality encoding/format (see FastqQualityFormat) used in a SAM/BAM or Fastq.
 * <p/>
 * To use this class, invoke the detect() method with a SAMFileReader or FastqReader, as appropriate.  The consumer is
 * responsible for closing readers.
 *
 * @author mccowan@broadinstitute.org
 */
public class QualityEncodingDetector {

    private QualityRecordAggregator qualityAggregator = new QualityRecordAggregator();

    /**
     * The maximum number of records over which the detector will iterate before making a determination, by default.
     */
    public static final long DEFAULT_MAX_RECORDS_TO_ITERATE = 10000;
    private static final Log log = Log.getInstance(QualityEncodingDetector.class);

    public enum FileContext {FASTQ, SAM}

    static class Range {
        final int low, high;

        Range(final int low, final int high) {
            this.low = low;
            this.high = high;
        }

        boolean contains(final int value) {
            return value <= high && value >= low;
        }
    }

    /**
     * Collection of data about the different quality formats and how they are interpreted.
     */
    enum QualityScheme {
        Phred(
                new Range(0, 93),           // Raw value range
                new Range(33, 126),         // ASCII value range
                asList(new Range(33, 58)),  // Ranges into which we expect at least one ASCII value to fall
                FastqQualityFormat.Standard // Associated quality format
        ),
        Solexa(
                new Range(-5, 62),
                new Range(59, 126),
                new ArrayList<Range>(),
                FastqQualityFormat.Solexa
        ),
        Illumina(
                new Range(0, 62),
                new Range(64, 126),
                new ArrayList<Range>(),
                FastqQualityFormat.Illumina
        );
        final Range rawRange, asciiRange;
        /**
         * Ranges into which we expect at least one value to fall if this formatting is being used.  For example, for
         * Standard encoding, we expect to at least one ASCII value between 33 and 58 (0 and 25); otherwise, it's
         * probably not Standard-encoded.
         */
        final List<Range> expectedAsciiRanges;
        final FastqQualityFormat qualityFormat;

        QualityScheme(final Range rawRange, final Range asciiRange, final List<Range> expectedAsciiRanges, final FastqQualityFormat qualityFormat) {
            this.rawRange = rawRange;
            this.asciiRange = asciiRange;
            this.expectedAsciiRanges = expectedAsciiRanges;
            this.qualityFormat = qualityFormat;
        }
    }

    /**
     * Collecting reads and their quality scores for later analysis.  Uses ASCII values since those are the inherently
     * "raw-est", read unmodified from the file.
     */
    private static class QualityRecordAggregator {
        private Set<Integer> observedAsciiQualities = new HashSet<Integer>();

        public Set<Integer> getObservedAsciiQualities() {
            return Collections.unmodifiableSet(observedAsciiQualities);
        }

        /**
         * Adds the FastqRecord's quality scores.
         */
        public void add(final FastqRecord fastqRecord) {
            addAsciiQuality(fastqRecord.getBaseQualityString().getBytes());
        }

        /**
         * Adds the SAMRecord's quality scores.
         * <p/>
         * Does not assume Phred quality encoding (for obvious reasons); getBaseQualityString() is used to read the
         * unmodified ASCII score.  To elaborate, SAMFileReader, which is generating these SAMRecords, builds the
         * SAMRecord by subtracting a value from each quality score and storing that transformed value internally.
         * Since we desire original scores here (whatever was in the file to begin with), we effectively undo this
         * transformation by asking SAMRecord to convert the quality back into the ASCII that was read in the file.
         */
        public void add(final SAMRecord samRecord) {
            addAsciiQuality(samRecord.getBaseQualityString().getBytes());
        }

        private void addAsciiQuality(final byte... asciiQualities) {
            for (final byte asciiQuality : asciiQualities) {
                observedAsciiQualities.add((int) asciiQuality);
            }
        }
    }

    /**
     * Adds the provided reader's records to the detector.
     * @return The number of records read
     */
    public long add(final long maxRecords, final FastqReader... readers) {
        final Iterator<FastqRecord> iterator = generateInterleavedFastqIterator(readers);
        long recordCount = 0;
        while (iterator.hasNext() && recordCount++ != maxRecords) {
            this.add(iterator.next());
        }
        log.debug(String.format("Read %s records from %s.", recordCount, Arrays.toString(readers)));
        return recordCount;
    }

    /**
     * Adds the provided reader's records to the detector.
     * @return The number of records read
     */
    public long add(final long maxRecords, final SAMFileReader reader) {
        final SAMRecordIterator iterator = reader.iterator();
        long recordCount = 0;
        try {
            while (iterator.hasNext() && recordCount++ != maxRecords) {
                this.add(iterator.next());
            }
            
            return recordCount;
        } finally {
            iterator.close();
        }
    }

    /**
     * Adds the provided record's qualities to the detector.
     */
    public void add(final FastqRecord fastqRecord) {
        this.qualityAggregator.add(fastqRecord);
    }

    /**
     * Adds the provided record's qualities to the detector.
     */
    public void add(final SAMRecord samRecord) {
        this.qualityAggregator.add(samRecord);
    }

    /**
     * Tests whether or not the detector can make a determination without guessing (i.e., if all but one quality format
     * can be excluded using established exclusion conventions).
     *
     * @return True if more than one format is possible after exclusions; false otherwise
     */
    public boolean isDeterminationAmbiguous() {
        return this.generateCandidateQualities().size() > 1;
    }

    /**
     * Processes collected quality data and applies rules to determine which quality formats are possible.
     * <p/>
     * Specifically, for each format's known range of possible values (its "quality scheme"), exclude formats if any
     * observed values fall outside of that range.  Additionally, exclude formats for which we expect to see at
     * least one quality in a range of values, but do not.  (For example, for Phred, we expect to eventually see
     * a value below 58.  If we never see such a value, we exclude Phred as a possible format.)
     */
    public EnumSet<FastqQualityFormat> generateCandidateQualities() {
        final EnumSet<FastqQualityFormat> candidateFormats = EnumSet.allOf(FastqQualityFormat.class);
        final Set<Integer> observedAsciiQualities = this.qualityAggregator.getObservedAsciiQualities();
        if (observedAsciiQualities.isEmpty())
            throw new PicardException("Cannot determine candidate qualities: no qualities found.");

        for (final QualityScheme scheme : QualityScheme.values()) {
            final Iterator<Integer> qualityBinIterator = observedAsciiQualities.iterator();
            final Collection<Range> remainingExpectedValueRanges = new ArrayList<Range>(scheme.expectedAsciiRanges);
            while (qualityBinIterator.hasNext()) {
                final int quality = qualityBinIterator.next();
                if (!scheme.asciiRange.contains(quality)) {
                    candidateFormats.remove(scheme.qualityFormat);
                }

                final Iterator<Range> expectedValueRangeIterator = remainingExpectedValueRanges.iterator();
                while (expectedValueRangeIterator.hasNext()) {
                    if (expectedValueRangeIterator.next().contains(quality)) {
                        expectedValueRangeIterator.remove();
                    }
                }
            }

            /**
             * We remove elements from this list as we observe values in the corresponding range; if the list isn't 
             * empty, we haven't seen a value in that range.  In other words, we haven't seen a value we expected.
             * Consequently, we remove the corresponding format from the running possibilities.
             */
            if (!remainingExpectedValueRanges.isEmpty()) {
                candidateFormats.remove(scheme.qualityFormat);
            }
        }

        return candidateFormats;
    }

    /**
     * Based on the quality scores accumulated in the detector as well as the context in which this guess applies (a BAM
     * or fastq), determines the best guess as to the quality format.
     * <p/>
     * This method does not exclude any quality formats based on observed quality values; even if there remain multiple
     * candidate qualities (see generateCandidateQualities()), picks a single one based on a high-level logic.
     *
     * @param context The type of file for which the guess is being made, Fastq or Bam
     */
    public FastqQualityFormat generateBestGuess(final FileContext context) {
        final EnumSet<FastqQualityFormat> possibleFormats = this.generateCandidateQualities();
        switch (possibleFormats.size()) {
            case 1:
                return possibleFormats.iterator().next();
            case 2:
                if (possibleFormats.equals(EnumSet.of(FastqQualityFormat.Illumina, FastqQualityFormat.Solexa))) {
                    return FastqQualityFormat.Illumina;
                } else if (possibleFormats.equals(EnumSet.of(FastqQualityFormat.Illumina, FastqQualityFormat.Standard))) {
                    switch (context) {
                        case FASTQ:
                            return FastqQualityFormat.Illumina;
                        case SAM:
                            return FastqQualityFormat.Standard;
                    }
                } else if (possibleFormats.equals(EnumSet.of(FastqQualityFormat.Standard, FastqQualityFormat.Solexa))) {
                    throw new PicardException("The quality format cannot be determined: both Phred and Solexa formats are possible; this application's logic does not handle this scenario.");
                } else throw new PicardException("Unreachable code.");
            case 3:
                throw new PicardException("The quality format cannot be determined: no formats were excluded.");
            case 0:
                throw new PicardException("The quality format cannot be determined: all formats were excluded.");
            default:
                throw new PicardException("Unreachable code.");
        }
    }

    /**
     * Interleaves FastqReader iterators so that serial-iteration of the result cycles between the constituent iterators.
     */
    private static Iterator<FastqRecord> generateInterleavedFastqIterator(final FastqReader... readers) {
        return new Iterator<FastqRecord>() {
            private Queue<Iterator<FastqRecord>> queue = new LinkedList<Iterator<FastqRecord>>();

            {
                for (final FastqReader reader : readers) {
                    queue.add(reader.iterator());
                }
            }

            public boolean hasNext() {
                // If this returns true, the head of the queue will have a next element
                while (!queue.isEmpty()) {
                    if (queue.peek().hasNext()) {
                        return true;
                    }
                    queue.poll();
                }
                return false;
            }

            public FastqRecord next() {
                if (!hasNext()) throw new NoSuchElementException();
                final Iterator<FastqRecord> i = queue.poll();
                final FastqRecord result = i.next();
                queue.offer(i);
                return result;
            }

            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    /**
     * Reads through the records in the provided fastq reader and uses their quality scores to determine the quality
     * format used in the fastq.
     *
     * @param readers    The fastq readers from which qualities are to be read; at least one must be provided
     * @param maxRecords The maximum number of records to read from the reader before making a determination (a guess,
     *                   so more records is better)
     * @return The determined quality format
     */
    public static FastqQualityFormat detect(final long maxRecords, final FastqReader... readers) {
        final QualityEncodingDetector detector = new QualityEncodingDetector();
        final long recordCount = detector.add(maxRecords, readers);
        log.debug(String.format("Read %s records from %s.", recordCount, Arrays.toString(readers)));
        return detector.generateBestGuess(FileContext.FASTQ);
    }

    public static FastqQualityFormat detect(final FastqReader... readers) {
        return detect(DEFAULT_MAX_RECORDS_TO_ITERATE, readers);
    }

    /**
     * Reads through the records in the provided SAM reader and uses their quality scores to determine the quality
     * format used in the SAM.
     *
     * @param reader     The SAM reader from which records are to be read
     * @param maxRecords The maximum number of records to read from the reader before making a determination (a guess,
     *                   so more records is better)
     * @return The determined quality format
     */
    public static FastqQualityFormat detect(final long maxRecords, final SAMFileReader reader) {
        final QualityEncodingDetector detector = new QualityEncodingDetector();
        final long recordCount = detector.add(maxRecords, reader);
        log.debug(String.format("Read %s records from %s.", recordCount, reader));
        return detector.generateBestGuess(FileContext.SAM);
    }

    public static FastqQualityFormat detect(final SAMFileReader reader) {
        return detect(DEFAULT_MAX_RECORDS_TO_ITERATE, reader);
    }
}