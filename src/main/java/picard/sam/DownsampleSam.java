/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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
package picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.DownsamplingIteratorFactory.Strategy;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.analysis.CollectQualityYieldMetrics.QualityYieldMetrics;
import picard.analysis.CollectQualityYieldMetrics.QualityYieldMetricsCollector;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.argumentcollections.ReferenceArgumentCollection;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Random;

/**
 * Class to randomly downsample a BAM file while respecting that we should either retain or discard
 * all of the reads for a template - i.e. all reads with the same name, whether first or second of
 * pair, secondary or supplementary, all travel together.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = DownsampleSam.USAGE_SUMMARY + DownsampleSam.USAGE_DETAILS,
        oneLineSummary = DownsampleSam.USAGE_SUMMARY,
        programGroup = SamOrBam.class)
@DocumentedFeature
public class DownsampleSam extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Downsample a SAM or BAM file.  ";
    static final String USAGE_DETAILS = "This tool applies a random downsampling algorithm to a SAM or BAM file to retain " +
            "only a random subset of the reads. Reads in a mate-pair are either both kept or both discarded. Reads marked as not primary " +
            "alignments are all discarded. Each read is given a probability P of being retained so that runs performed with the exact " +
            "same input in the same order and with the same value for RANDOM_SEED will produce the same results." +
            "All reads for a template are kept or discarded as a unit, with the goal of retaining reads" +
            "from PROBABILITY * input templates. While this will usually result in approximately " +
            "PROBABILITY * input reads being retained also, for very small PROBABILITIES this may not " +
            "be the case.\n" +
            "A number of different downsampling strategies are supported using the STRATEGY option:\n\n" +
            "ConstantMemory: " + DownsamplingIteratorFactory.CONSTANT_MEMORY_DESCRPTION + "\n\n" +
            "HighAccuracy: " + DownsamplingIteratorFactory.HIGH_ACCURACY_DESCRIPTION + "\n\n" +
            "Chained: " + DownsamplingIteratorFactory.CHAINED_DESCRIPTION + "\n\n" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar DownsampleSam \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=downsampled.bam" +
            "</pre>" +
            "<hr />";
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to downsample.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output, downsampled, SAM or BAM file to write.")
    public File OUTPUT;

    @Argument(shortName="S", doc="The downsampling strategy to use. See usage for discussion.")
    public Strategy STRATEGY = Strategy.ConstantMemory;

    @Argument(shortName = "R", doc = "Random seed to use if deterministic behavior is desired.  " +
            "Setting to null will cause multiple invocations to produce different results.")
    public Integer RANDOM_SEED = 1;

    @Argument(shortName = "P", doc = "The probability of keeping any individual read, between 0 and 1.")
    public double PROBABILITY = 1;

    @Argument(shortName = "A", doc = "The accuracy that the downsampler should try to achieve if the selected strategy supports it. " +
            "Note that accuracy is never guaranteed, but some strategies will attempt to provide accuracy within the requested bounds." +
            "Higher accuracy will generally require more memory.")
    public double ACCURACY = 0.0001;

    @Argument(shortName = "M", doc = "The file to write metrics to (QualityYieldMetrics)", optional=true)
    public File METRICS_FILE;

    private final Log log = Log.getInstance(DownsampleSam.class);

    public static void main(final String[] args) {
        new DownsampleSam().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        // Warn the user if they are running with P=1; 0 <= P <= 1 is checked by the DownsamplingIteratorFactory
        if (PROBABILITY == 1) {
            log.warn("Running DownsampleSam with PROBABILITY=1! This will likely just recreate the input file.");
        }

        final Integer seed = RANDOM_SEED == null ? new Random().nextInt() : RANDOM_SEED;
        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(SamInputResource.of(INPUT));
        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(in.getFileHeader(), true, OUTPUT);
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e7, "Wrote");
        final DownsamplingIterator iterator = DownsamplingIteratorFactory.make(in, STRATEGY, PROBABILITY, ACCURACY, seed);
        final QualityYieldMetricsCollector metricsCollector = new QualityYieldMetricsCollector(true, false, false);

        while (iterator.hasNext()) {
            final SAMRecord rec = iterator.next();
            out.addAlignment(rec);
            if (METRICS_FILE != null) metricsCollector.acceptRecord(rec, null);
            progress.record(rec);
        }

        out.close();
        CloserUtil.close(in);
        final NumberFormat fmt = new DecimalFormat("0.00%");
        log.info("Finished downsampling.");
        log.info("Kept ", iterator.getAcceptedCount(), " out of ", iterator.getSeenCount(), " reads (", fmt.format(iterator.getAcceptedFraction()), ").");

        if (METRICS_FILE != null) {
            final MetricsFile<QualityYieldMetrics, Integer> metricsFile = getMetricsFile();
            metricsCollector.finish();
            metricsCollector.addMetricsToFile(metricsFile);
            metricsFile.write(METRICS_FILE);
        }

        return 0;
    }

    @Override
    protected ReferenceArgumentCollection makeReferenceArgumentCollection() {
        // Override to allow "R" to be hijacked for "RANDOM_SEED"
        return new ReferenceArgumentCollection() {
            @Argument(doc = "The reference sequence file.", optional=true, common=false)
            public File REFERENCE_SEQUENCE;

            @Override
            public File getReferenceFile() {
                return REFERENCE_SEQUENCE;
            }
        };
    }

}
