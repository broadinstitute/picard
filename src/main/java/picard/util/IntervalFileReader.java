package picard.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import picard.PicardException;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * Utility class for reading interval files in BED or Picard interval_list format.
 *
 * <p>BED parsing delegates to htsjdk's {@link BEDCodec}, decoding each line individually
 * over a {@link BufferedReader}. Driving the codec ourselves (instead of going through
 * {@link htsjdk.tribble.AbstractFeatureReader}, which is path-based) lets us read directly
 * from {@code /dev/stdin} and other streams without staging through a temp file. {@link BEDCodec}
 * still does the actual parsing, so {@code track}/{@code browser}/{@code #} header lines and
 * other quirks are handled exactly as the rest of htsjdk handles them.
 *
 * <p>BED coordinate conventions: BED uses 0-based half-open coordinates. {@link BEDCodec}
 * converts them to the 1-based closed intervals used throughout Picard / htsjdk.
 */
public final class IntervalFileReader {

    /** Recognized interval file formats. */
    public enum IntervalFileFormat { INTERVAL_LIST, BED, UNKNOWN }

    /** Result of format detection: the detected format and, when format is UNKNOWN, the first data line. */
    public record FormatDetectionResult(IntervalFileFormat format, String firstLine) {}

    private static final Log log = Log.getInstance(IntervalFileReader.class);

    private IntervalFileReader() {} // utility class — no instances

    /**
     * Opens a {@link BufferedReader} for {@code file}. Recognizes {@code /dev/stdin} as a
     * request to read from {@link System#in} (this respects {@link System#setIn} for tests
     * and avoids re-opening stdin as a file). Other paths go through
     * {@link IOUtil#openFileForBufferedReading(File)} so gzipped inputs are handled transparently.
     */
    private static BufferedReader openBuffered(final File file) {
        if ("/dev/stdin".equals(file.getPath())) {
            return new BufferedReader(new InputStreamReader(System.in));
        }
        return IOUtil.openFileForBufferedReading(file);
    }

    /**
     * Detects whether a reader contains interval_list or BED content by inspecting the first
     * significant (non-empty, non-comment) line.  The reader's position is NOT reset after this
     * call — callers must {@link BufferedReader#mark} and {@link BufferedReader#reset} as needed.
     *
     * <ul>
     *   <li>Lines starting with {@code @} indicate interval_list (SAM-style header).</li>
     *   <li>Lines with ≥3 tab-separated fields indicate BED.</li>
     *   <li>Comment lines starting with {@code #} are skipped.</li>
     * </ul>
     */
    public static FormatDetectionResult detectIntervalFormat(final BufferedReader reader) throws IOException {
        String line;
        while ((line = reader.readLine()) != null) {
            final String trimmed = line.trim();
            if (trimmed.isEmpty() || trimmed.startsWith("#")) {
                continue;
            }
            if (trimmed.startsWith("@")) {
                return new FormatDetectionResult(IntervalFileFormat.INTERVAL_LIST, null);
            } else if (trimmed.split("\t").length >= 3) {
                return new FormatDetectionResult(IntervalFileFormat.BED, null);
            } else {
                return new FormatDetectionResult(IntervalFileFormat.UNKNOWN, trimmed);
            }
        }
        return new FormatDetectionResult(IntervalFileFormat.UNKNOWN, null);
    }

    /**
     * Loads an {@link IntervalList} from {@code file}, auto-detecting whether the content is
     * BED or interval_list format.  Returns a uniqued {@link IntervalList}.
     *
     * <p>The file is opened once and sniffed via {@link BufferedReader#mark}/{@link BufferedReader#reset},
     * so {@code /dev/stdin} works without temp-file staging.
     *
     * @param file       BED or interval_list file to read
     * @param dictionary sequence dictionary used when parsing BED-format content
     * @throws PicardException on read error or unrecognized format
     */
    public static IntervalList loadIntervals(final File file, final SAMSequenceDictionary dictionary) {
        try (final BufferedReader reader = openBuffered(file)) {
            // Sniffing only needs to peek as far as the first significant line; 8 KB is plenty.
            reader.mark(8 * 1024);
            final FormatDetectionResult detected = detectIntervalFormat(reader);
            reader.reset();
            return switch (detected.format()) {
                case INTERVAL_LIST -> IntervalList.fromReader(reader).uniqued();
                case BED -> {
                    final SAMFileHeader header = new SAMFileHeader();
                    header.setSequenceDictionary(dictionary);
                    yield fromBed(reader, header, false, false).uniqued();
                }
                case UNKNOWN -> throw new PicardException(
                        "Unrecognized interval file format. Expected interval_list (lines starting with @) " +
                        "or BED (≥3 tab-separated fields). First data line: " + detected.firstLine());
            };
        } catch (final IOException e) {
            throw new PicardException("Error reading intervals from " + file, e);
        }
    }

    /** Calls {@link #fromBed(File, SAMFileHeader, boolean, boolean)} with both flags false. */
    public static IntervalList fromBed(final File file, final SAMFileHeader header) {
        return fromBed(file, header, false, false);
    }

    /**
     * Parses BED records from {@code file} into an {@link IntervalList}.  Recognizes
     * {@code /dev/stdin} as a request to read from {@link System#in}; other paths are
     * opened via {@link IOUtil#openFileForBufferedReading(File)} (which handles gzip).
     *
     * @see #fromBed(BufferedReader, SAMFileHeader, boolean, boolean)
     */
    public static IntervalList fromBed(final File file, final SAMFileHeader header,
                                       final boolean dropMissingContigs,
                                       final boolean keepLengthZeroIntervals) {
        try (final BufferedReader reader = openBuffered(file)) {
            return fromBed(reader, header, dropMissingContigs, keepLengthZeroIntervals);
        } catch (final IOException e) {
            throw new PicardException("Error reading BED file: " + file, e);
        }
    }

    /**
     * Parses BED records from {@code reader} into an {@link IntervalList} by feeding each line
     * to {@link BEDCodec#decode(String)}.  Empty lines and BED header lines (lines starting with
     * {@code #}, {@code track}, or {@code browser}) are silently skipped, matching the codec's
     * built-in behavior.  Coordinates are returned in 1-based closed format as required by
     * Picard / htsjdk.
     *
     * @param reader                  reader positioned at the start of the BED data
     * @param header                  SAMFileHeader whose sequence dictionary is used for validation
     * @param dropMissingContigs      if true, records whose contig is absent from the dictionary are silently dropped;
     *                                if false, such records throw a {@link PicardException}
     * @param keepLengthZeroIntervals if true, BED records where chromStart == chromEnd are included;
     *                                if false, they are silently dropped
     * @return an unsorted, non-uniqued {@link IntervalList}
     * @throws PicardException if a record references a contig absent from the dictionary and dropMissingContigs is false
     */
    public static IntervalList fromBed(final BufferedReader reader, final SAMFileHeader header,
                                       final boolean dropMissingContigs,
                                       final boolean keepLengthZeroIntervals) {
        final BEDCodec codec = new BEDCodec();
        final IntervalList intervalList = new IntervalList(header);
        int lengthZeroIntervals = 0;

        try {
            String line;
            while ((line = reader.readLine()) != null) {
                final BEDFeature bedFeature = codec.decode(line);
                if (bedFeature == null) continue; // empty / # / track / browser

                final String sequenceName = bedFeature.getContig();
                // BEDCodec converts to 1-based start and 1-based inclusive end
                final int start = bedFeature.getStart();
                final int end   = bedFeature.getEnd();
                // if the interval has no name, use null to avoid confusion with empty-string names
                final String name = bedFeature.getName().isEmpty() ? null : bedFeature.getName();
                final boolean isNegativeStrand = bedFeature.getStrand() == Strand.NEGATIVE;

                final SAMSequenceRecord sequenceRecord = header.getSequenceDictionary().getSequence(sequenceName);

                if (null == sequenceRecord) {
                    if (dropMissingContigs) {
                        log.info(String.format("Dropping interval with missing contig: %s:%d-%d", sequenceName, start, end));
                        continue;
                    }
                    throw new PicardException(String.format("Sequence '%s' was not found in the sequence dictionary", sequenceName));
                } else if (start < 1) {
                    throw new PicardException(String.format("Start on sequence '%s' was less than one: %d", sequenceName, start));
                } else if (sequenceRecord.getSequenceLength() < start) {
                    throw new PicardException(String.format("Start on sequence '%s' was past the end: %d < %d", sequenceName, sequenceRecord.getSequenceLength(), start));
                } else if ((end == 0 && start != 1) || end < 0) {
                    throw new PicardException(String.format("End on sequence '%s' was less than one: %d", sequenceName, end));
                } else if (sequenceRecord.getSequenceLength() < end) {
                    throw new PicardException(String.format("End on sequence '%s' was past the end: %d < %d", sequenceName, sequenceRecord.getSequenceLength(), end));
                } else if (end < start - 1) {
                    throw new PicardException(String.format("On sequence '%s', end < start-1: %d <= %d", sequenceName, end, start));
                }

                if (start == end + 1) {
                    lengthZeroIntervals++;
                    if (!keepLengthZeroIntervals) {
                        log.info(String.format("Skipping writing length zero interval at %s:%d-%d.", sequenceName, start, end));
                        continue;
                    }
                }

                intervalList.add(new Interval(sequenceName, start, end, isNegativeStrand, name));
            }
        } catch (final IOException e) {
            throw new PicardException("Error reading BED data", e);
        }

        if (lengthZeroIntervals == 0) {
            log.info("No input regions had length zero, so none were skipped.");
        } else if (!keepLengthZeroIntervals) {
            log.info(String.format("Skipped writing a total of %d entries with length zero in the input file.", lengthZeroIntervals));
        } else {
            log.warn(String.format("Input file had %d entries with length zero.", lengthZeroIntervals));
        }

        return intervalList;
    }
}
