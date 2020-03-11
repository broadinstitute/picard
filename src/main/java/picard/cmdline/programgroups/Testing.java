package picard.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import picard.util.help.HelpConstants;

/**
 * For internal test purposes only.
 *
 * @author nhomer
 */
public class Testing implements CommandLineProgramGroup {
    @Override
    public String getName() { return HelpConstants.DOC_CAT_TEST; }
    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_TEST_SUMMARY; }
}
