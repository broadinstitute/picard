package picard.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import picard.util.help.HelpConstants;

/**
 * Tools that process genomic intervals in various formats.
 *
 * @author nhomer
*/
public class IntervalsManipulationProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() { return HelpConstants.DOC_CAT_INTERVALS_MANIPULATION ; }
    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_INTERVALS_MANIPULATION_SUMMARY; }
}
