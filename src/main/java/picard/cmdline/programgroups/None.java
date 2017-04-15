package picard.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;

/**
* @author nhomer
*/
public class None implements CommandLineProgramGroup {
    @Override
    public String getName() { return "Miscellaneous Tools"; }
    @Override
    public String getDescription() { return "A set of miscellaneous tools."; }
}
