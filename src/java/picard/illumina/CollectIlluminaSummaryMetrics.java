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

package picard.illumina;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import picard.PicardException;
import picard.analysis.SinglePassSamProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Illumina;

import java.io.File;

/**
 * Program to collect Coverage Summary Metrics of SAM files for data sequenced by Illumina products.
 * It uses our best approximation of the filters that Illumina uses which means:
 * <p/>
 * 1. Examine PF reads only
 * 2. Examine reads that are not marked as duplicate only
 * 3. Examine mapped reads only (it is unclear if Illumina does this or not, but without this we cannot do the next point)
 * 4. Only count bases that are present in two of mated reads, once. (For this we need TLEN from the sam record, which required that reads are mapped)
 *
 * Program assumes that all reads in the SAM file are of equal length.
 *
 * @author farjoun@broadinstitute.org
 */
@CommandLineProgramProperties(
        usage = "Program to collect Coverage Summary Metrics of SAM files for data sequenced by Illumina products. " +
                "It uses our best approximation of the filters that Illumina uses which means:\n\n" +
                "1. Examine PF reads only" +
                "2. Examine reads that are not marked as duplicate only" +
                "3. Examine mapped reads only (it is unclear if Illumina does this or not, but without this we cannot do the next point)" +
                "4. Only count bases that are present in two of mated reads, once. (For this we need TLEN from the sam record, which required that reads are mapped)" +
                "\nProgram assumes that all reads in the SAM file are of equal length, and that if a read is marked as \"mated\", its mate will really exist.\n",
        usageShort = "Collects summary metrics according to Illumina specifications.",
        programGroup = Illumina.class
)
public class CollectIlluminaSummaryMetrics extends SinglePassSamProgram {

    @Option(doc = "IntervalList describing the target interval of the sequencing experiment (for calculating average coverage). " +
            "Use null to use reference from BAM as interval (whole genome)", optional = true)
    public File TARGET_REGION = null;

    private final IlluminaSummaryMetrics metrics = new IlluminaSummaryMetrics();

    private static final Log log = Log.getInstance(CollectIlluminaSummaryMetrics.class);

    private Integer readLength = 0;

    /** Stock main method for a command line program. */
    public static void main(final String[] argv) {
        new CollectIlluminaSummaryMetrics().instanceMainWithExit(argv);
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {

        if (TARGET_REGION == null) {
            metrics.TARGET_TERRITORY = header.getSequenceDictionary().getReferenceLength();
        } else {
            IOUtil.assertFileIsReadable(TARGET_REGION);
            IOUtil.assertFileIsWritable(OUTPUT);
            final IntervalList targetInterval = IntervalList.fromFile(TARGET_REGION);
            metrics.TARGET_TERRITORY = targetInterval.getUniqueBaseCount();
        }
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {

        if (readLength == 0) {
            readLength = rec.getReadLength();
        } else {
            if (rec.getReadLength() != readLength) {
                final String message=String.format("This program only works with uniform read lengths. First record had length %d. Current record %s, has length %d", readLength, rec.getReadName(), rec.getReadLength());
                log.error(message);
                throw new PicardException(message);
            }
        }

        // Not interested in counting unmapped, secondary or supplemental reads.
        if (rec.getReadUnmappedFlag() || rec.isSecondaryOrSupplementary()) {
            return;
        }

        metrics.TOTAL_READS++;

        final int length = rec.getReadLength();
        metrics.TOTAL_BASES += length;

        final boolean isPfRead = !rec.getReadFailsVendorQualityCheckFlag();

        if (isPfRead) {
            metrics.PF_READS++;
            metrics.PF_BASES += length;
        } else {
            return;
        }

        // From here on read must have passed PF
        if (rec.getDuplicateReadFlag()) {
            metrics.DUPLICATE_READS++;
            return;
        }

        // Reads are not duplicated from here on.
        final int tLen = Math.abs(rec.getInferredInsertSize());
        if (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag() && !rec.getReadUnmappedFlag() && tLen != 0) {
            // If mated,both reads are mapped, and insert-size can be calculated, only examine first in pair, assume equal read lengths, and account for overlaps
            if (rec.getFirstOfPairFlag()) {
                // Since read is mated with both mate mapped, and insert-size available, we will only consider pair on first mate.
                if (rec.getReadNegativeStrandFlag() ^ rec.getMateNegativeStrandFlag()) {
                    // If mated reads are in opposite directions:
                    metrics.TOTAL_ILLUMINA_BASES += Math.min(2 * length, tLen);
                } else {
                    // If mated reads are in same direction:
                    metrics.TOTAL_ILLUMINA_BASES += length + Math.min(length, tLen);
                }
            }
        } else {
            // Otherwise add each read on its own.
            metrics.TOTAL_ILLUMINA_BASES += rec.getReadLength();
        }
    }

    // Do not need unmapped reads at end of file.
    @Override
    protected boolean usesNoRefReads() { return false; }

    @Override
    protected void finish() {
        // Calculate some derived metrics
        metrics.READ_LENGTH = readLength;
        metrics.AVERAGE_ILMN_DEPTH = metrics.TOTAL_ILLUMINA_BASES / (double) metrics.TARGET_TERRITORY;

        // Output the file
        final MetricsFile<IlluminaSummaryMetrics, Integer> metricsFile = getMetricsFile();
        metricsFile.addMetric(metrics);
        metricsFile.write(OUTPUT);
    }

    /** A set of metrics used to describe the general quality of a BAM file */
    public static class IlluminaSummaryMetrics extends MetricBase {

        /** The total number of reads in the input file */
        public int TOTAL_READS = 0;

        /** The number of reads that are PF - pass filter */
        public int PF_READS = 0;

        /** The number of duplicate, PF reads */
        public int DUPLICATE_READS = 0;

        /** The total number of bases in all reads */
        public long TOTAL_BASES;

        /** The total number of bases in all PF reads */
        public long PF_BASES = 0;

        /** The number of bases from PF, unique (non-duplicate) reads, only counting for overlapping bases once */
        public long TOTAL_ILLUMINA_BASES = 0;

        /** The average depth considering only bases from PF, unique (non-duplicate) reads, only counting for overlapping bases once */
        public double AVERAGE_ILMN_DEPTH = 0;

        /** The read length */
        public long READ_LENGTH = 0;

        /** The size of the target region */
        public long TARGET_TERRITORY = 0;

        public IlluminaSummaryMetrics() {}

        public IlluminaSummaryMetrics(final int TOTAL_READS, final int PF_READS, final int DUPLICATE_READS, final long TOTAL_BASES, final long PF_BASES, final long TOTAL_ILLUMINA_BASES, final double AVERAGE_ILMN_DEPTH, final long READ_LENGTH, final long TARGET_TERRITORY) {
            this.TOTAL_READS = TOTAL_READS;
            this.PF_READS = PF_READS;
            this.DUPLICATE_READS = DUPLICATE_READS;
            this.TOTAL_BASES = TOTAL_BASES;
            this.PF_BASES = PF_BASES;
            this.TOTAL_ILLUMINA_BASES = TOTAL_ILLUMINA_BASES;
            this.AVERAGE_ILMN_DEPTH = AVERAGE_ILMN_DEPTH;
            this.READ_LENGTH = READ_LENGTH;
            this.TARGET_TERRITORY = TARGET_TERRITORY;
        }
    }
}
