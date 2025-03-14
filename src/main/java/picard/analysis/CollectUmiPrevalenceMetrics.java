/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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

package picard.analysis;

import htsjdk.samtools.*;
import htsjdk.samtools.filter.*;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import picard.filter.CountingFilterWrapper;
import picard.filter.CountingMapQFilter;
import picard.filter.CountingPairedFilter;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

import static picard.cmdline.StandardOptionDefinitions.MINIMUM_MAPPING_QUALITY_SHORT_NAME;

/**
 * A CLP that, tallies the number of UMIs in a duplicate set across the entire bam and produces a histogram.
 *
 * @author Yossi Farjoun
 */
@DocumentedFeature
@ExperimentalFeature
@CommandLineProgramProperties(
        summary = CollectUmiPrevalenceMetrics.USAGE_SUMMARY + CollectUmiPrevalenceMetrics.USAGE_DETAILS,
        oneLineSummary = CollectUmiPrevalenceMetrics.USAGE_SUMMARY,
        programGroup = DiagnosticsAndQCProgramGroup.class
)
public class CollectUmiPrevalenceMetrics extends CommandLineProgram {

    static final String USAGE_SUMMARY = "Tally the counts of UMIs in duplicate sets within a bam. \n";
    static final String USAGE_DETAILS = "<p>" +
            "This tool collects the Histogram of the number of duplicate sets that contain a given number of UMIs. " +
            "Understanding this distribution can help understand the role that the UMIs have in the determination of " +
            "consensus sets, the risk of UMI collisions, and of spurious reads that result from uncorrected UMIs.";


    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input (indexed) BAM/CRAM file.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Write metrics to this file")
    public File OUTPUT;

    @Argument(shortName = MINIMUM_MAPPING_QUALITY_SHORT_NAME, doc = "minimal value for the mapping quality of the reads to be used in the estimation.", optional = true)
    public Integer MINIMUM_MQ = 30;

    @Argument(doc = "Barcode SAM tag.", optional = true)
    public String BARCODE_TAG = "RX";

    @Argument(doc = "Barcode Quality SAM tag.", optional = true)
    public String BARCODE_BQ = "QX";

    @Argument(shortName = "BQ", doc = "minimal value for the base quality of all the bases in a molecular barcode, for it to be used.", optional = true)
    public Integer MINIMUM_BARCODE_BQ = 30;

    @Argument(shortName = "FUR", doc = "Whether to filter unpaired reads from the input.", optional = true)
    public boolean FILTER_UNPAIRED_READS = true;

    @Argument(
            fullName = "PROGRESS_STEP_INTERVAL",
            doc = "The interval between which progress will be displayed.",
            optional = true
    )
    public int PROGRESS_STEP_INTERVAL = 1_000_000;


    private static final Log log = Log.getInstance(CollectUmiPrevalenceMetrics.class);

    @Override
    protected int doWork() {

        IOUtil.assertFileIsReadable(INPUT);

        // get an iterator to reads that overlap the heterozygous sites
        final Histogram<Integer> umiCount = new Histogram<>("numUmis", "duplicateSets");
        final CountingPairedFilter countingPairedFilter = new CountingPairedFilter();
        final CountingFilterWrapper countingAlignedFilter = new CountingFilterWrapper(new AlignedFilter(true));
        final CountingMapQFilter countingMapQFilter = new CountingMapQFilter(MINIMUM_MQ);
        final CountingFilterWrapper countingSecondaryOrSupplementaryFilter =
                new CountingFilterWrapper(new SecondaryOrSupplementaryFilter());
        final ProgressLogger progress = new ProgressLogger(log, PROGRESS_STEP_INTERVAL, "examined", "duplicate sets");

        try (SamReader in = SamReaderFactory.makeDefault()
                .referenceSequence(REFERENCE_SEQUENCE)
                .open(INPUT)) {

            IOUtil.assertFileIsWritable(OUTPUT);

            final SAMRecordIterator samRecordIterator = in.iterator();
            final List<SamRecordFilter> samFilters = CollectionUtil.makeList(
                    countingAlignedFilter,
                    countingMapQFilter,
                    countingSecondaryOrSupplementaryFilter
            );
            if (FILTER_UNPAIRED_READS) {
                samFilters.add(countingPairedFilter);
            }

            final FilteringSamIterator filteredSamRecordIterator = new FilteringSamIterator(samRecordIterator, new AggregateFilter(samFilters));
            log.info("Queried BAM, getting duplicate sets.");

            // get duplicate iterator from iterator above
            final DuplicateSetIterator duplicateSets = new DuplicateSetIterator(filteredSamRecordIterator, in.getFileHeader(), false, null, log);

            log.info("Starting iteration on duplicate sets");

            set:
            while (duplicateSets.hasNext()) {

                final DuplicateSet set = duplicateSets.next();
                final SAMRecord setRep = set.getRepresentative();

                progress.record(setRep);
                Set<String> barcodes = new HashSet<>();
                for (final SAMRecord read : set.getRecords()) {
                    if (!read.hasAttribute(BARCODE_TAG)) {
                        log.warn("no barcode tag!");
                        continue;
                    }
                    if (read.hasAttribute(BARCODE_BQ)) {
                        final String barcodeBQ = read.getStringAttribute(BARCODE_BQ).replace(" ", "");
                        ;
                        final byte[] bytes = SAMUtils.fastqToPhred(barcodeBQ);
                        final boolean badQuality = IntStream.range(0, bytes.length).map(i -> bytes[i]).anyMatch(q -> q < MINIMUM_BARCODE_BQ);
                        if (badQuality) {
                            log.warn("bad quality barcode");
                            continue;
                        }
                    }
                    barcodes.add(read.getStringAttribute(BARCODE_TAG));
                }
                umiCount.increment(barcodes.size(), 1);
            }
        } catch (IOException e) {
            throw new RuntimeException("Problem while reading file: " + INPUT, e);
        }
        log.info("Iteration done. Emitting metrics.");
        log.info(String.format("Processed %d sets", progress.getCount()));
        log.info(String.format("Filtered %d unpaired reads", countingPairedFilter.getFilteredRecords()));
        log.info(String.format("Filtered %d unaligned reads", countingAlignedFilter.getFilteredRecords()));
        log.info(String.format("Filtered %d low mapQ reads", countingMapQFilter.getFilteredRecords()));
        log.info(String.format("Filtered %d Secondary or Supplementary reads", countingSecondaryOrSupplementaryFilter.getFilteredRecords()));
        // Emit metrics
        final MetricsFile<?, Integer> metricsFile = getMetricsFile();
        metricsFile.addHistogram(umiCount);
        metricsFile.write(OUTPUT);
        return 0;
    }
}
