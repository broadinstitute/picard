package picard.sam;

import htsjdk.samtools.DownsamplingIteratorFactory;
import htsjdk.samtools.DownsamplingIteratorFactory.Strategy;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.DownsamplingIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Random;

/**
 * Class to randomly downsample a BAM file while respecting that we should either get rid
 * of both ends of a pair or neither end of the pair!
 */
@CommandLineProgramProperties(
        usage = "Randomly down-sample a SAM or BAM file to retain only a subset of the reads in the file. " +
                "All reads for a templates are kept or discarded as a unit, with the goal or retaining reads" +
                "from PROBABILITY * input templates. While this will usually result in approximately " +
                "PROBABILITY * input reads being retained also, for very small PROBABILITIES this may not " +
                "be the case.\n" +
                "A number of different downsampling strategies are supported using the STRATEGY option:\n\n" +
                "ConstantMemory: " + DownsamplingIteratorFactory.CONSTANT_MEMORY_DESCRPTION + "\n\n" +
                "HighAccuracy: " + DownsamplingIteratorFactory.HIGH_ACCURACY_DESCRIPTION + "\n\n" +
                "Chained: " + DownsamplingIteratorFactory.CHAINED_DESCRIPTION + "\n",
        usageShort = "Down-sample a SAM or BAM file to retain a random subset of the reads",
        programGroup = SamOrBam.class
)
public class DownsampleSam extends CommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to downsample.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output, downsampled, SAM or BAM file to write.")
    public File OUTPUT;

    @Option(shortName="S", doc="The downsampling strategy to use. See usage for discussion.")
    public Strategy STRATEGY = Strategy.ConstantMemory;

    @Option(shortName = "R", doc = "Random seed to use if reproducibilty is desired.  " +
            "Setting to null will cause multiple invocations to produce different results.")
    public Integer RANDOM_SEED = 1;

    @Option(shortName = "P", doc = "The probability of keeping any individual read, between 0 and 1.")
    public double PROBABILITY = 1;

    @Option(shortName = "A", doc = "The accuracy that the downsampler should try to achieve if the selected strategy supports it. " +
            "Note that accuracy is never guaranteed, but some strategies will attempt to provide accuracy within the requested bounds." +
            "Higher accuracy will generally require more memory.")
    public double ACCURACY = 0.0001;

    private final Log log = Log.getInstance(DownsampleSam.class);

    public static void main(final String[] args) {
        new DownsampleSam().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final Random r = RANDOM_SEED == null ? new Random() : new Random(RANDOM_SEED);
        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(in.getFileHeader(), true, OUTPUT);
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e7, "Wrote");
        final DownsamplingIterator iterator = DownsamplingIteratorFactory.make(INPUT, STRATEGY, PROBABILITY, ACCURACY, RANDOM_SEED);

        while (iterator.hasNext()) {
            final SAMRecord rec = iterator.next();
            out.addAlignment(rec);
            progress.record(rec);
        }

        out.close();
        CloserUtil.close(in);
        final NumberFormat fmt = new DecimalFormat("0.00%");
        log.info("Finished downsampling.");
        log.info("Kept ", iterator.getAcceptedCount(), " out of ", iterator.getSeenCount(), " reads (", fmt.format(iterator.getAcceptedFraction()), ").");

        return 0;
    }
}
