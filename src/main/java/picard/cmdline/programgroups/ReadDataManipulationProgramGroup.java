package picard.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import picard.util.help.HelpConstants;

/**
 * Tools that manipulate read data in SAM, BAM or CRAM format
 */
public class ReadDataManipulationProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_READ_DATA_MANIPULATION; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_READ_DATA_MANIPULATION_SUMMARY; }
}