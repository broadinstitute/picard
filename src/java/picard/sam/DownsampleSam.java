package picard.sam;

import htsjdk.samtools.DiscardLevel;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.HashDecisionProjector;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.Usage;

import java.io.File;
import java.util.Random;

/**
 * Class to randomly downsample a BAM file while respecting that we should either get rid
 * of both ends of a pair or neither end of the pair!
 */
public class DownsampleSam extends CommandLineProgram {
    @Usage public final String USAGE = getStandardUsagePreamble() + " Randomly down-sample a SAM or BAM file to retain " +
            "a random subset of the reads. Mate-pairs and reads marked as not primary are either all kept or all " +
            "discarded. Each read is given a probability P of being retained - results with the exact same input in " +
            "the same order and with the same value for RANDOM_SEED will produce the same results.";

    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The input SAM or BAM file to downsample.")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output, downsampled, SAM or BAM file to write.")
    public File OUTPUT;

    @Option(shortName="R", doc="Random seed to use if reproducibilty is desired.  " +
            "Setting to null will cause multiple invocations to produce different results.")
    public Integer RANDOM_SEED = 1;

    @Option(shortName="P", doc="The probability of keeping any individual read, between 0 and 1.")
    public double PROBABILITY = 1;

    @Option(shortName="D", doc="Whether or not to discard non-primary alignments.")
    public DiscardLevel DISCARD_LEVEL = DiscardLevel.NONE;

    private final Log log = Log.getInstance(DownsampleSam.class);

    public static void main(final String[] args) {
        new DownsampleSam().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final SAMFileReader in = new SAMFileReader(INPUT);
        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(in.getFileHeader(), true, OUTPUT);
        final HashDecisionProjector decision = new HashDecisionProjector(RANDOM_SEED == null ? new Random().nextInt() : RANDOM_SEED, PROBABILITY);

        long total = 0;
        long kept  = 0;
        
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e7, "Read");

        for (final SAMRecord rec : in) {
            if ( (DISCARD_LEVEL == DiscardLevel.SECONDARY && rec.getNotPrimaryAlignmentFlag()) || (DISCARD_LEVEL == DiscardLevel.ALL && rec.isSecondaryOrSupplementary())) continue;

            if (!rec.isSecondaryOrSupplementary()) ++total;

            final String key = rec.getReadName();

            if (decision.decide(key)) {
                out.addAlignment(rec);
                if (!rec.isSecondaryOrSupplementary()) ++kept;
            }

            progress.record(rec);
        }

        out.close();
        log.info("Finished! Kept " + kept + " out of " + total + " primary reads.");

        return 0;
    }
}
