/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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

package picard.sam.markduplicates.util;

import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;
import htsjdk.samtools.MergingSamRecordIterator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamFileHeaderMerger;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Histogram;
import picard.PicardException;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.sam.DuplicationMetrics;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Abstract class that holds parameters and methods common to classes that perform duplicate
 * detection and/or marking within SAM/BAM files.
 *
 * @author Nils Homer
 */
public abstract class AbstractMarkDuplicatesCommandLineProgram extends AbstractOpticalDuplicateFinderCommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "One or more input SAM or BAM files to analyze. Must be coordinate sorted.")
    public List<String> INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "The output file to write marked records to")
    public File OUTPUT;

    @Option(shortName = "M",
            doc = "File to write duplication metrics to")
    public File METRICS_FILE;

    @Option(doc = "If true do not write duplicates to the output file instead of writing them with appropriate flags set.")
    public boolean REMOVE_DUPLICATES = false;

    @Option(shortName = StandardOptionDefinitions.ASSUME_SORTED_SHORT_NAME,
            doc = "If true, assume that the input file is coordinate sorted even if the header says otherwise.")
    public boolean ASSUME_SORTED = false;

    @Option(shortName = "DS", doc = "The scoring strategy for choosing the non-duplicate among candidates.")
    public ScoringStrategy DUPLICATE_SCORING_STRATEGY = ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH;

    
    @Option(shortName = StandardOptionDefinitions.PROGRAM_RECORD_ID_SHORT_NAME,
            doc = "The program record ID for the @PG record(s) created by this program. Set to null to disable " +
                    "PG record creation.  This string may have a suffix appended to avoid collision with other " +
                    "program record IDs.",
            optional = true)
    public String PROGRAM_RECORD_ID = "MarkDuplicates";

    @Option(shortName = "PG_VERSION",
            doc = "Value of VN tag of PG record to be created. If not specified, the version will be detected automatically.",
            optional = true)
    public String PROGRAM_GROUP_VERSION;

    @Option(shortName = "PG_COMMAND",
            doc = "Value of CL tag of PG record to be created. If not supplied the command line will be detected automatically.",
            optional = true)
    public String PROGRAM_GROUP_COMMAND_LINE;

    @Option(shortName = "PG_NAME",
            doc = "Value of PN tag of PG record to be created.")
    public String PROGRAM_GROUP_NAME = getClass().getSimpleName();

    @Option(shortName = "CO",
            doc = "Comment(s) to include in the output file's header.",
            optional = true)
    public List<String> COMMENT = new ArrayList<String>();

    /** The program groups that have been seen during the course of examining the input records. */
    protected final Set<String> pgIdsSeen = new HashSet<String>();

    /**
     * We have to re-chain the program groups based on this algorithm.  This returns the map from existing program group ID
     * to new program group ID.
     */
    protected Map<String, String> getChainedPgIds(final SAMFileHeader outputHeader) {
        final Map<String, String> chainedPgIds;
        // Generate new PG record(s)
        if (PROGRAM_RECORD_ID != null) {
            final SAMFileHeader.PgIdGenerator pgIdGenerator = new SAMFileHeader.PgIdGenerator(outputHeader);
            if (PROGRAM_GROUP_VERSION == null) {
                PROGRAM_GROUP_VERSION = this.getVersion();
            }
            if (PROGRAM_GROUP_COMMAND_LINE == null) {
                PROGRAM_GROUP_COMMAND_LINE = this.getCommandLine();
            }
            chainedPgIds = new HashMap<String, String>();
            for (final String existingId : this.pgIdsSeen) {
                final String newPgId = pgIdGenerator.getNonCollidingId(PROGRAM_RECORD_ID);
                chainedPgIds.put(existingId, newPgId);
                final SAMProgramRecord programRecord = new SAMProgramRecord(newPgId);
                programRecord.setProgramVersion(PROGRAM_GROUP_VERSION);
                programRecord.setCommandLine(PROGRAM_GROUP_COMMAND_LINE);
                programRecord.setProgramName(PROGRAM_GROUP_NAME);
                programRecord.setPreviousProgramGroupId(existingId);
                outputHeader.addProgramRecord(programRecord);
            }
        } else {
            chainedPgIds = null;
        }
        return chainedPgIds;
    }

    /**
     * Writes the metrics given by the libraryIdGenerator to the METRICS_FILE.
     *
     * @param libraryIdGenerator
     */
    protected void finalizeAndWriteMetrics(final LibraryIdGenerator libraryIdGenerator) {
        final Map<String, DuplicationMetrics> metricsByLibrary = libraryIdGenerator.getMetricsByLibraryMap();
        final Histogram<Short> opticalDuplicatesByLibraryId = libraryIdGenerator.getOpticalDuplicatesByLibraryIdMap();
        final Map<String, Short> libraryIds = libraryIdGenerator.getLibraryIdsMap();

        // Write out the metrics
        final MetricsFile<DuplicationMetrics, Double> file = getMetricsFile();
        for (final Map.Entry<String, DuplicationMetrics> entry : metricsByLibrary.entrySet()) {
            final String libraryName = entry.getKey();
            final DuplicationMetrics metrics = entry.getValue();

            metrics.READ_PAIRS_EXAMINED = metrics.READ_PAIRS_EXAMINED / 2;
            metrics.READ_PAIR_DUPLICATES = metrics.READ_PAIR_DUPLICATES / 2;

            // Add the optical dupes to the metrics
            final Short libraryId = libraryIds.get(libraryName);
            if (libraryId != null) {
                final Histogram<Short>.Bin bin = opticalDuplicatesByLibraryId.get(libraryId);
                if (bin != null) {
                    metrics.READ_PAIR_OPTICAL_DUPLICATES = (long) bin.getValue();
                }
            }
            metrics.calculateDerivedMetrics();
            file.addMetric(metrics);
        }

        if (metricsByLibrary.size() == 1) {
            file.setHistogram(metricsByLibrary.values().iterator().next().calculateRoiHistogram());
        }

        file.write(METRICS_FILE);
    }

    /** Little class used to package up a header and an iterable/iterator. */
    public static final class SamHeaderAndIterator {
        public final SAMFileHeader header;
        public final CloseableIterator<SAMRecord> iterator;

        public SamHeaderAndIterator(final SAMFileHeader header, final CloseableIterator<SAMRecord> iterator) {
            this.header = header;
            this.iterator = iterator;
        }
    }

    /**
     * Since this may read it's inputs more than once this method does all the opening
     * and checking of the inputs.
     */
    protected SamHeaderAndIterator openInputs() {
        final List<SAMFileHeader> headers = new ArrayList<SAMFileHeader>(INPUT.size());
        final List<SamReader> readers = new ArrayList<SamReader>(INPUT.size());

        for (final String input : INPUT) {
            SamReader reader = SamReaderFactory.makeDefault()
                .enable(SamReaderFactory.Option.EAGERLY_DECODE)
                .open(SamInputResource.of(input));
            final SAMFileHeader header = reader.getFileHeader();

            if (!ASSUME_SORTED && header.getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
                throw new PicardException("Input file " + input + " is not coordinate sorted.");
            }

            headers.add(header);
            readers.add(reader);
        }

        if (headers.size() == 1) {
            return new SamHeaderAndIterator(headers.get(0), readers.get(0).iterator());
        } else {
            final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate, headers, false);
            final MergingSamRecordIterator iterator = new MergingSamRecordIterator(headerMerger, readers, ASSUME_SORTED);
            return new SamHeaderAndIterator(headerMerger.getMergedHeader(), iterator);
        }
    }

    /**
     * Looks through the set of reads and identifies how many of the duplicates are
     * in fact optical duplicates, and stores the data in the instance level histogram.
     */
    public static void trackOpticalDuplicates(List<? extends ReadEnds> ends,
                                              final OpticalDuplicateFinder opticalDuplicateFinder,
                                              final LibraryIdGenerator libraryIdGenerator) {
        boolean hasFR = false, hasRF = false;

        // Check to see if we have a mixture of FR/RF
        for (final ReadEnds end : ends) {
            if (ReadEnds.FR == end.orientationForOpticalDuplicates) {
                hasFR = true;
            } else if (ReadEnds.RF == end.orientationForOpticalDuplicates) {
                hasRF = true;
            }
        }

        // Check if we need to partition since the orientations could have changed
        if (hasFR && hasRF) { // need to track them independently
            // Variables used for optical duplicate detection and tracking
            final List<ReadEnds> trackOpticalDuplicatesF = new ArrayList<ReadEnds>();
            final List<ReadEnds> trackOpticalDuplicatesR = new ArrayList<ReadEnds>();

            // Split into two lists: first of pairs and second of pairs, since they must have orientation and same starting end
            for (final ReadEnds end : ends) {
                if (ReadEnds.FR == end.orientationForOpticalDuplicates) {
                    trackOpticalDuplicatesF.add(end);
                } else if (ReadEnds.RF == end.orientationForOpticalDuplicates) {
                    trackOpticalDuplicatesR.add(end);
                } else {
                    throw new PicardException("Found an unexpected orientation: " + end.orientation);
                }
            }

            // track the duplicates
            trackOpticalDuplicates(trackOpticalDuplicatesF, opticalDuplicateFinder, libraryIdGenerator.getOpticalDuplicatesByLibraryIdMap());
            trackOpticalDuplicates(trackOpticalDuplicatesR, opticalDuplicateFinder, libraryIdGenerator.getOpticalDuplicatesByLibraryIdMap());
        } else { // No need to partition
            AbstractMarkDuplicatesCommandLineProgram.trackOpticalDuplicates(ends, opticalDuplicateFinder, libraryIdGenerator.getOpticalDuplicatesByLibraryIdMap());
        }
    }

    /**
     * Looks through the set of reads and identifies how many of the duplicates are
     * in fact optical duplicates, and stores the data in the instance level histogram.
     * 
     * We expect only reads with FR or RF orientations, not a mixture of both. 
     * 
     * In PCR duplicate detection, a duplicates can be a have FR and RF when fixing the orientation order to the first end of the mate.  In
     * optical duplicate detection, we do not consider them duplicates if one read as FR ann the other RF when we order orientation by the
     * first mate sequenced (read #1 of the pair).
     */
    private static void trackOpticalDuplicates(final List<? extends OpticalDuplicateFinder.PhysicalLocation> list,
                                               final OpticalDuplicateFinder opticalDuplicateFinder,
                                               final Histogram<Short> opticalDuplicatesByLibraryId) {
        final boolean[] opticalDuplicateFlags = opticalDuplicateFinder.findOpticalDuplicates(list);

        int opticalDuplicates = 0;
        for (final boolean b : opticalDuplicateFlags) if (b) ++opticalDuplicates;
        if (opticalDuplicates > 0) {
            opticalDuplicatesByLibraryId.increment(list.get(0).getLibraryId(), opticalDuplicates);
        }
    }
}
