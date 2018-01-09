package picard.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import picard.util.help.HelpConstants;

/**
 * Miscellaneous tools, e.g. those that aid in data streaming
 */
public class OtherProgramGroup  implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_OTHER; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_OTHER_SUMMARY; }
}
