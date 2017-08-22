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
package picard.analysis.directed;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.analysis.MetricAccumulationLevel;
import picard.analysis.SinglePassSamProgram;
import picard.cmdline.programgroups.Metrics;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;
import picard.sam.markduplicates.util.ReadEndsForMateCigar;
import picard.sam.markduplicates.util.SamRecordWithOrdinalAndSetDuplicateReadFlag;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Collects a set of insert coverage metrics from a sam or bam file.
 *
 * @author Yossi Farjoun
 */
@CommandLineProgramProperties(
        summary = CollectInsertCoverageMetrics.USAGE_SUMMARY + CollectInsertCoverageMetrics.USAGE_DETAILS,
        oneLineSummary = CollectInsertCoverageMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
@DocumentedFeature
public class CollectInsertCoverageMetrics extends SinglePassSamProgram {
    static final String USAGE_SUMMARY = "blah";
    static final String USAGE_DETAILS = "blah";
    private SAMFileHeader header = null;
    private long recordCount = 0;
    private Map<Interval, Integer> coverageByInterval = new HashMap<>();
    final private OpticalDuplicateFinder opticalDuplicateFinder = new OpticalDuplicateFinder();

    IntervalCoverageMetrics metrics = new IntervalCoverageMetrics();
    // Overlap detector for finding overlaps between the reads and the baits (and the near bait space)
    private final OverlapDetector<Interval> intervalDetector = new OverlapDetector<>(0, 0);

    private Map<Interval, Double> intervalToGc = null;

    @Argument(shortName = "INTS", doc = "An interval list file that contains the locations of the baits used.", minElements = 1)
    public List<File> INTERVALS;

    @Argument(shortName = "MI", doc = "Maximal InsertSize of inserts to be used for the analysis")
    private int MAX_INSERT_SIZE = 750;

    @Override
    protected boolean requiresReference() {
        return true;
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(OUTPUT);
        for (final File interval : INTERVALS) IOUtil.assertFileIsReadable(interval);

        final IntervalList intervals = IntervalList.fromFiles(INTERVALS);
        intervalDetector.addAll(intervals.getIntervals(), intervals.getIntervals());

        intervalToGc = new HashMap<>();
        try(final ReferenceSequenceFile refFile = new IndexedFastaSequenceFile(REFERENCE_SEQUENCE)) {
            for (final Interval interval : intervals) {
                coverageByInterval.put(interval, 0);
                final ReferenceSequence rs = refFile.getSubsequenceAt(interval.getContig(), interval.getStart(), interval.getEnd());
                intervalToGc.put(interval, SequenceUtil.calculateGc(rs.getBases()));
            }
        }catch (final Exception e) {
            throw new RuntimeException(e);
        }
        this.header = header;
    }

    @Override
    protected void acceptRead(final SAMRecord record, final ReferenceSequence ref) {
        /** Adds information about an individual SAMRecord to the statistics. */
        // Just ignore secondary alignments altogether
        if (record.getNotPrimaryAlignmentFlag()) return;

        // Cache some things, and compute the total number of bases aligned in the record.
        final boolean mappedInPair = record.getReadPairedFlag() &&
                !record.getReadUnmappedFlag() &&
                !record.getMateUnmappedFlag() &&
                !record.getSupplementaryAlignmentFlag() &&
                record.getProperPairFlag();

        if (!mappedInPair) return;

        if (record.getInferredInsertSize() == 0 || Math.abs(record.getInferredInsertSize()) > MAX_INSERT_SIZE) return;

        // READ Based Metrics
        this.metrics.TOTAL_READS++;
        SamPairUtil.PairOrientation po = SamPairUtil.getPairOrientation(record);
        if (!record.getReadFailsVendorQualityCheckFlag()) { // only reads that pass vendor's filters
            this.metrics.PF_READS++;
            if (!record.getDuplicateReadFlag()) { // ignore duplicates for unique reads/bases
                this.metrics.PF_UNIQUE_READS++;
                switch (po) {
                    case FR:
                        this.metrics.PF_UNIQUE_INNI_READS++;
                        break;
                    case RF:
                        this.metrics.PF_UNIQUE_OUTIE_READS++;
                        break;
                    case TANDEM:
                        this.metrics.PF_UNIQUE_TANDEM_READS++;
                        break;
                    default:
                        throw new RuntimeException("Unpossible");
                }
            }
        }

        ///////////////////////////////////////////////////////////////////
        // Non-PF reads can be totally ignored beyond this point
        ///////////////////////////////////////////////////////////////////
        if (record.getReadFailsVendorQualityCheckFlag()) return;
        if (po != SamPairUtil.PairOrientation.FR) return;
        if (!record.getFirstOfPairFlag()) return;

        // create an interval consisting of insert between the two reads
        final ReadEndsForMateCigar readEndsforRecord = new ReadEndsForMateCigar(this.header, new SamRecordWithOrdinalAndSetDuplicateReadFlag(record, recordCount++),
                opticalDuplicateFinder, (short) 0);

        final Interval insert = new Interval(record.getReferenceName(),
                Math.min(readEndsforRecord.read1Coordinate, readEndsforRecord.read2Coordinate),
                Math.max(readEndsforRecord.read1Coordinate, readEndsforRecord.read2Coordinate));

        // Prefetch the list of target and bait overlaps here as they're needed multiple times.
        final Collection<Interval> overlappingIntervals = intervalDetector.getOverlaps(insert);

        for (final Interval interval : overlappingIntervals) {
            final Integer currentCount = coverageByInterval.get(interval);
            coverageByInterval.put(interval, currentCount + 1);
        }

    }

    @Override
    protected void finish() {
        final FormatUtil fmt = new FormatUtil();
        try (final PrintWriter out = new PrintWriter(OUTPUT)) {
            out.println("chrom\tstart\tend\tlength\tname\tpct_gc\tinserts");

            for (final Map.Entry<Interval, Integer> entry : this.coverageByInterval.entrySet()) {
                final Interval interval = entry.getKey();
                final Integer coverage = entry.getValue();
                final double intervalGc = intervalToGc.get(interval);

                out.println(interval.getContig() + "\t" +
                        interval.getStart() + "\t" +
                        interval.getEnd() + "\t" +
                        interval.length() + "\t" +
                        interval.getName() + "\t" +
                        fmt.format(intervalGc) + "\t" +
                        fmt.format(coverage)
                );
            }

        } catch (final IOException ioe) {
            throw new RuntimeIOException(ioe);
        }
    }
}
