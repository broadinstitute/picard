package picard.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import picard.PicardException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * Utility class for reading interval files in BED or Picard interval_list format.
 *
 * <p>This class provides format-detection and generic loading (auto-detecting BED vs interval_list).
 * BED parsing delegates to htsjdk's {@link BEDCodec} / {@link AbstractFeatureReader} (Tribble),
 * preserving the same robust behavior used by {@link BedToIntervalList}.
 *
 * <p>BED coordinate conventions: BED uses 0-based half-open coordinates.  {@link BEDCodec}
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
     * BED or interval_list format.
     *
     * <p>Returns a sorted, uniqued {@link IntervalList}.
     *
     * @param file       BED or interval_list file to read
     * @param dictionary sequence dictionary used when parsing BED-format content
     * @throws PicardException on read error or unrecognized format
     */
    public static IntervalList loadIntervals(final File file, final SAMSequenceDictionary dictionary) {
        final FormatDetectionResult detected;
        try (final BufferedReader reader = new BufferedReader(new FileReader(file))) {
            reader.mark(8 * 1024);
            detected = detectIntervalFormat(reader);
        } catch (final IOException e) {
            throw new PicardException("Error detecting format of " + file, e);
        }

        return switch (detected.format()) {
            case INTERVAL_LIST -> {
                try (final BufferedReader reader = new BufferedReader(new FileReader(file))) {
                    yield IntervalList.fromReader(reader).uniqued();
                } catch (final IOException e) {
                    throw new PicardException("Error reading interval_list from " + file, e);
                }
            }
            case BED -> {
                final SAMFileHeader header = new SAMFileHeader();
                header.setSequenceDictionary(dictionary);
                yield fromBed(file, header).uniqued();
            }
            case UNKNOWN -> throw new PicardException(
                    "Unrecognized interval file format. Expected interval_list (lines starting with @) " +
                    "or BED (≥3 tab-separated fields). First data line: " + detected.firstLine());
        };
    }

    /** Calls {@link #fromBed(File, SAMFileHeader, boolean, boolean)} with both flags false. */
    public static IntervalList fromBed(final File file, final SAMFileHeader header) {
        return fromBed(file, header, false, false);
    }

    /**
     * Parses BED records from {@code file} into an {@link IntervalList}.
     *
     * <p>Uses htsjdk's {@link BEDCodec} / {@link AbstractFeatureReader} for parsing, which
     * silently skips {@code track} and {@code browser} lines and handles real-world BED quirks.
     * Coordinates are returned in 1-based closed format as required by Picard / htsjdk.
     *
     * @param file                    BED file to parse
     * @param header                  SAMFileHeader whose sequence dictionary is used for validation
     * @param dropMissingContigs      if true, records whose contig is absent from the dictionary are silently dropped;
     *                                if false, such records throw a {@link PicardException}
     * @param keepLengthZeroIntervals if true, BED records where chromStart == chromEnd are included;
     *                                if false, they are silently dropped
     * @return an unsorted, non-uniqued {@link IntervalList}
     * @throws PicardException if a record references a contig absent from the dictionary and dropMissingContigs is false
     */
    public static IntervalList fromBed(final File file, final SAMFileHeader header,
                                       final boolean dropMissingContigs,
                                       final boolean keepLengthZeroIntervals) {
        final IntervalList intervalList = new IntervalList(header);
        int lengthZeroIntervals = 0;

        final FeatureReader<BEDFeature> bedReader = AbstractFeatureReader.getFeatureReader(
                file.getAbsolutePath(), new BEDCodec(), false);
        try {
            final CloseableTribbleIterator<BEDFeature> iterator = bedReader.iterator();
            while (iterator.hasNext()) {
                final BEDFeature bedFeature = iterator.next();
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
            throw new PicardException("Error reading BED file: " + file, e);
        } finally {
            CloserUtil.close(bedReader);
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
