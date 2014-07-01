package picard.cmdline.programgroups;

import picard.cmdline.CommandLineProgramGroup;

/**
* @author nhomer
*/
public class Crsp implements CommandLineProgramGroup {
    @Override
    public String getName() { return "Crsp"; }
    @Override
    public String getDescription() { return "Tools specific to the CRSP pipeline."; }
}
