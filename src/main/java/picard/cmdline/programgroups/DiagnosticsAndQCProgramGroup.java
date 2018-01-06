package picard.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import picard.util.help.HelpConstants;

/**
 * Tools that collect sequencing quality-related and comparative metrics
 */
public class DiagnosticsAndQCProgramGroup  implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_DIAGNOSTICS_AND_QC; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_DIAGNOSTICS_AND_QC_SUMMARY; }
}