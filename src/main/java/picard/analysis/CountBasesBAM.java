package picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import picard.util.RExecutor;


@CommandLineProgramProperties(
        summary = "Counts As in a SAM/BAM/CRAM file",
        oneLineSummary = "Count As in a SAM/BAM/CRAM file",
        programGroup = DiagnosticsAndQCProgramGroup.class
)

@DocumentedFeature
public class CountBasesBAM{

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM/BAM file to revert the state of.")
    public File INPUT;

    private long count = 0;

    public void apply(SAMRecord read) {
        String basesString = read.getReadString();
        for (int i = 0; i < basesString.length(); i++) {
            if (basesString.charAt(i) == 'A') {
                count += 1;
            }
        }
    }

    public Object onTraversalSuccess() {
        System.out.println(count);
        return count;
    }


}
