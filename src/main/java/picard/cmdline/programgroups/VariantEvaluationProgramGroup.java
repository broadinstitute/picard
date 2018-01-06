package picard.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import picard.util.help.HelpConstants;

/**
 * Tools that evaluate and refine variant calls, e.g. by adding annotations that are not calculated during variant calling
 */
public class VariantEvaluationProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_VARIANT_EVALUATION; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_VARIANT_EVALUATION_SUMMARY; }
}