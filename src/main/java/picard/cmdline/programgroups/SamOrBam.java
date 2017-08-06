package picard.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;

/**
* @author nhomer
*/
public class SamOrBam implements CommandLineProgramGroup {
    @Override
    public String getName() { return "SAM/BAM"; }
    @Override
    public String getDescription() { return "Tools for manipulating SAM, BAM, or related data."; }
}
