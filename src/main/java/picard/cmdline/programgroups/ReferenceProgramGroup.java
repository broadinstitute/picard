package picard.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import picard.util.help.HelpConstants;

/**
 * Tools that analyze and manipulate FASTA format references
 */
public class ReferenceProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_REFERENCE; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_REFERENCE_SUMMARY; }
}
