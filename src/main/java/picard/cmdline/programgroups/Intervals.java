package picard.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;

/**
* @author nhomer
*/
public class Intervals implements CommandLineProgramGroup {
    @Override
    public String getName() { return "Interval Tools"; }
    @Override
    public String getDescription() { return "Tools for manipulating Picard interval lists."; }
}
