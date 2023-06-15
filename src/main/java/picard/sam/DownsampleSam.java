/*
 * The MIT License
 *
 * Copyright (c) 2015-2020 The Broad Institute
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

import htsjdk.samtools.DownsamplingIterator;
import htsjdk.samtools.DownsamplingIteratorFactory;
import htsjdk.samtools.DownsamplingIteratorFactory.Strategy;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.analysis.CollectQualityYieldMetrics.QualityYieldMetrics;
import picard.analysis.CollectQualityYieldMetrics.QualityYieldMetricsCollector;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.argumentcollections.ReferenceArgumentCollection;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

/**
 * <h3>Summary</h3>
 * This tool applies a downsampling algorithm to a SAM or BAM file to retain only a (deterministically random) subset of
 * the reads. Reads from the same template (e.g. read-pairs, secondary and supplementary reads) are all either kept or
 * discarded as a unit, with the goal of retaining reads from <code>PROBABILITY * (input templates)</code>. The results
 * will contain approximately <code>PROBABILITY * (input reads)</code>, however for very small
 * probabilities this may not be the case.
 *
 * A number of different downsampling strategies are supported using the {@link #STRATEGY} option:
 * <dl>
 * <dt>ConstantMemory</dt>
 * <dd>
 *     Downsamples a stream or file of SAMRecords using a hash-projection strategy such that it can run in constant memory.
 *     The downsampling is stochastic, and therefore the actual retained proportion will vary around the requested proportion. Due
 *     to working in fixed memory this strategy is good for large inputs, and due to the stochastic nature the accuracy of this strategy
 *     is highest with a high number of output records, and diminishes at low output volumes.
 * </dd>
 * <dt>HighAccuracy</dt>
 * <dd>
 *     Attempts (but does not guarantee) to provide accuracy up to a specified limit. Accuracy is defined as emitting
 *     a proportion of reads as close to the requested proportion as possible. In order to do so this strategy requires
 *     memory that is proportional to the number of template names in the incoming stream of reads, and will thus require
 *     large amounts of memory when running on large input files.
 * </dd>
 * <dt>Chained</dt>
 * <dd>
 *     Attempts to provide a compromise strategy that offers some of the advantages of both the ConstantMemory and HighAccuracy strategies.
 *     Uses a ConstantMemory strategy to downsample the incoming stream to approximately the desired proportion, and then a HighAccuracy
 *     strategy to finish. Works in a single pass, and will provide accuracy close to (but often not as good as) HighAccuracy while requiring
 *     memory proportional to the set of reads emitted from the ConstantMemory strategy to the HighAccuracy strategy. Works well when downsampling
 *     large inputs to small proportions (e.g. downsampling hundreds of millions of reads and retaining only 2%. Should be accurate 99.9% of the time
 *     when the input contains more than 50,000 templates (read names). For smaller inputs, HighAccuracy is recommended instead.
 * </dd>
 * </dl>
 *
 * The number of records written can be output to a {@link QualityYieldMetrics} metrics file via the {@link #METRICS_FILE}.
 *
 * <h3>Usage examples:</h3>
 * <h4>Downsample file, keeping about 10% of the reads</h4>
 * <pre>
 * java -jar picard.jar DownsampleSam \
 *       I=input.bam \
 *       O=downsampled.bam \
 *       P=0.1
 * </pre>
 *
 * <h4>Downsample file, keeping 2% of the reads </h4>
 * <pre>
 * java -jar picard.jar DownsampleSam \
 *       I=input.bam \
 *       O=downsampled.bam \
 *       STRATEGY=Chained \
 *       P=0.02 \
 *       ACCURACY=0.0001
 * </pre>
 *
 * <h4>Downsample file, keeping 0.001% of the reads (may require more memory)</h4>
 * <pre>
 * java -jar picard.jar DownsampleSam \
 *       I=input.bam \
 *       O=downsampled.bam \
 *       STRATEGY=HighAccuracy \
 *       P=0.00001 \
 *       ACCURACY=0.0000001
 * </pre>
 *
 * @author Tim Fennell
 */

@CommandLineProgramProperties(
        summary = DownsampleSam.USAGE_SUMMARY + DownsampleSam.USAGE_DETAILS,
        oneLineSummary = DownsampleSam.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature
public class DownsampleSam extends CommandLineProgram {

    final String PG_PROGRAM_NAME = getClass().getSimpleName();

    static final String USAGE_SUMMARY = "Downsample a SAM or BAM file.";
    static final String USAGE_DETAILS = "This tool applies a downsampling algorithm to a SAM or BAM file to retain " +
            "only a (deterministically random) subset of the reads. Reads from the same template (e.g. read-pairs, secondary " +
            "and supplementary reads) are all either kept or discarded as a unit, with the goal of retaining reads" +
            "from PROBABILITY * input templates. The results will contain approximately " +
            "PROBABILITY * input reads, however for very small PROBABILITIES this may not " +
            "be the case.\n" +
            "A number of different downsampling strategies are supported using the STRATEGY option:\n\n" +
            "ConstantMemory:\n " + DownsamplingIteratorFactory.CONSTANT_MEMORY_DESCRPTION + "\n" +
            "HighAccuracy:\n " + DownsamplingIteratorFactory.HIGH_ACCURACY_DESCRIPTION + "\n" +
            "Chained:\n " + DownsamplingIteratorFactory.CHAINED_DESCRIPTION + "\n" +
            "<h3>Usage examples:</h3>\n" +
            "<h4>Downsample file, keeping about 10% of the reads</h4>\n"+
            "\n"+
            "java -jar picard.jar DownsampleSam \\\n" +
            "      I=input.bam \\\n" +
            "      O=downsampled.bam \\\n" +
            "      P=0.2\n"+
            "\n" +
            "<h3>Downsample file, keeping about 2% of the reads </h3>\n"+
            "\n" +
            "java -jar picard.jar DownsampleSam \\\n" +
            "      I=input.bam \\\n" +
            "      O=downsampled.bam \\\n" +
            "      STRATEGY=Chained \\\n" +
            "      P=0.02 \\\n" +
            "      ACCURACY=0.0001\n" +
            "\n" +
            "<h3>Downsample file, keeping about 0.001% of the reads (may require more memory)</h3>\n"+
            "\n" +
            "java -jar picard.jar DownsampleSam \\\n" +
            "      I=input.bam \\\n" +
            "      O=downsampled.bam \\\n" +
            "      STRATEGY=HighAccuracy \\\n" +
            "      P=0.00001 \\\n" +
            "      ACCURACY=0.0000001\n";
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to downsample.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output, downsampled, SAM, BAM or CRAM file to write.")
    public File OUTPUT;

    @Argument(shortName="S", doc="The downsampling strategy to use. See usage for discussion.")
    public Strategy STRATEGY = Strategy.ConstantMemory;

    @Argument(shortName = "R", doc = "Random seed used for deterministic results. " +
            "Setting to null will cause multiple invocations to produce different results.  The header if the file will be checked for any previous runs " +
            "of DownsampleSam.  If DownsampleSam has been run before on this data with the same seed, the seed will be updated in a deterministic fashion " +
            "so the DownsampleSam will perform correctly, and still deterministically.")
    public Integer RANDOM_SEED = 1;

    @Argument(shortName = "P", doc = "The probability of keeping any individual read, between 0 and 1.")
    public double PROBABILITY = 1;

    @Argument(shortName = "A", doc = "The accuracy that the downsampler should try to achieve if the selected strategy supports it. " +
            "Note that accuracy is never guaranteed, but some strategies will attempt to provide accuracy within the requested bounds." +
            "Higher accuracy will generally require more memory.")
    public double ACCURACY = 0.0001;

    @Argument(shortName = "M", doc = "The metrics file (of type QualityYieldMetrics) which will contain information about the downsampled file.", optional=true)
    public File METRICS_FILE;

    private final Log log = Log.getInstance(DownsampleSam.class);

    public static final String RANDOM_SEED_TAG = "rs";

    @Override
    protected String[] customCommandLineValidation() {
        if (PROBABILITY < 0 || PROBABILITY > 1)
            return new String[]{"Downsampling requires 0<=PROBABILITY<=1. Found invalid value: " + PROBABILITY};

        return super.customCommandLineValidation();
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        // Warn the user if they are running with P=1 or P=0 (which are legal, but odd)
        if (PROBABILITY == 1) {
            log.warn("Running DownsampleSam with PROBABILITY=1! This will likely just recreate the input file.");
        }

        if (PROBABILITY == 0) {
            log.warn("Running DownsampleSam with PROBABILITY=0! This will create an empty file.");
        }

        if (RANDOM_SEED == null) {
            RANDOM_SEED = new Random().nextInt();
            log.warn(String.format(
                    "Drawing a random seed because RANDOM_SEED was not set. Set RANDOM_SEED to %s to reproduce these results in the future.", RANDOM_SEED));
        }

        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(SamInputResource.of(INPUT));
        final SAMFileHeader header = in.getFileHeader().clone();

        if (STRATEGY == Strategy.ConstantMemory || STRATEGY == Strategy.Chained) {
            //if running using ConstantMemory or Chained strategy, need to check if we have previously run with the same random seed
            //collect previously used seeds
            final Integer userSeed = RANDOM_SEED;
            final Set<Integer> previousSeeds = new HashSet<>();
            for (final SAMProgramRecord pg : header.getProgramRecords()) {
                if (pg.getProgramName() != null && pg.getProgramName().equals(PG_PROGRAM_NAME)) {
                    final String previousSeedString = pg.getAttribute(RANDOM_SEED_TAG);
                    if (previousSeedString == null) {
                        /* The previous seed was not recorded.  In this case, the current seed may be the same as the previous seed,
                        so we will change it to a randomly selected seed, which is very likely to be unique
                         */
                        RANDOM_SEED = new Random(pg.hashCode()).nextInt();
                        log.warn("DownsampleSam has been run before on this data, but the previous seed was not recorded.  The used seed will be changed to minimize the chance of using" +
                                " the same seed as in a previous run.");
                    } else {
                        previousSeeds.add(Integer.parseInt(previousSeedString));
                    }
                }
            }

            final Random rnd = new Random(RANDOM_SEED);
            while (previousSeeds.contains(RANDOM_SEED)) {
                RANDOM_SEED = rnd.nextInt();
                log.warn("DownsampleSam has been run before on this data with the seed " + RANDOM_SEED + ".  The random seed will be changed to avoid using the " +
                        "same seed as previously.");
            }
            if (!userSeed.equals(RANDOM_SEED)) {
                log.warn("RANDOM_SEED has been changed to " + RANDOM_SEED + ".");
            }
        }

        final SAMProgramRecord pgRecord = getPGRecord(header);
        pgRecord.setAttribute(RANDOM_SEED_TAG, RANDOM_SEED.toString());
        header.addProgramRecord(pgRecord);
        final SAMFileWriter out = new SAMFileWriterFactory().makeWriter(header, true, OUTPUT, referenceSequence.getReferenceFile());
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e7, "Wrote");
        final DownsamplingIterator iterator = DownsamplingIteratorFactory.make(in, STRATEGY, PROBABILITY, ACCURACY, RANDOM_SEED);
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
