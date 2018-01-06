package picard.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import picard.util.help.HelpConstants;

/**
 * Tools that process sequencing machine data, e.g. Illumina base calls, and detect sequencing level attributes, e.g. adapters
 */
public class BaseCallingProgramGroup implements CommandLineProgramGroup{

    @Override
    public String getName() { return HelpConstants.DOC_CAT_BASE_CALLING; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_BASE_CALLING_SUMMARY; }
}
