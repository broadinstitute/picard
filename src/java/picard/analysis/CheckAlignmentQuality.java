package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;
import picard.cmdline.programgroups.SamOrBam;

import javax.xml.crypto.Data;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.List;

/**
 * A minor data quality check against the AlignmentSummaryMetrics file which writes out an assessment of "SUCCEEDED" or "DATA_QUALITY_FAILED"
 * to a specified output file.
 */
@CommandLineProgramProperties(
        usage = "Checks that alignment quality meets certain thresholds",
        usageShort = "Checks alignment quality",
        programGroup = Metrics.class
)
public class CheckAlignmentQuality extends CommandLineProgram {
    private static final Log log = Log.getInstance(CheckAlignmentQuality.class);

    static final String USAGE = "Reads in an alignment summary metrics file and writes an assessment to a file";

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Option(shortName = "T", optional = true)
    public Integer DATA_QUALITY_THRESHOLD = 1000;

    @Option(shortName = "P")
    public Double PCT_PF_READS_THRESHOLD = 0.1;

    public static void main(final String[] args) {new CheckAlignmentQuality().instanceMainWithExit(args);}

    public static final String DataQualitySuccess = "SUCCEEDED";
    public static final String DataQualityFailed = "DATA_QUALITY_FAILED";

    @Override
    public int doWork() {
        String dataQuality = DataQualitySuccess;

        try {
            final List<AlignmentSummaryMetrics> alignmentSummaryMetrics = (List<AlignmentSummaryMetrics>) MetricsFile.readBeans(INPUT);
            final AlignmentSummaryMetrics alignmentSummaryMetric = alignmentSummaryMetrics.get(alignmentSummaryMetrics.size() - 1);

            if (alignmentSummaryMetric.PF_READS < DATA_QUALITY_THRESHOLD) {
                log.error("DataQuality Failure: " + 100 * PCT_PF_READS_THRESHOLD + " or fewer PF reads in " + INPUT.getAbsolutePath());
                dataQuality = DataQualityFailed;
            }

            if (alignmentSummaryMetric.PCT_PF_READS < PCT_PF_READS_THRESHOLD) {
                log.error("Data Quality Failure: < " + DATA_QUALITY_THRESHOLD + " PF reads in " + INPUT.getAbsolutePath());
                dataQuality = DataQualityFailed;
            }
        } catch (final Exception e) {
            log.error("Failed to parse Alignment Summary Metrics file: " + INPUT.getAbsolutePath());
            dataQuality = DataQualityFailed;
        }

        try {
            final PrintWriter writer = new PrintWriter(OUTPUT);
            writer.write(dataQuality);
            writer.close();
        } catch (final FileNotFoundException fnfe) {
            throw new PicardException("Unable to open output file " + OUTPUT.getAbsolutePath() + " for writing", fnfe);
        }

        return 0;
    }
}
