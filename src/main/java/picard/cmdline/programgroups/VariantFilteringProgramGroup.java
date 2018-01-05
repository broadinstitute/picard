package picard.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import picard.util.help.HelpConstants;

/**
 * Tools that filter variants
 */
public class VariantFilteringProgramGroup  implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_VARIANT_FILTERING; }

    @Override
    public String getDescription() { return  HelpConstants.DOC_CAT_VARIANT_FILTERING_SUMMARY; }
}

