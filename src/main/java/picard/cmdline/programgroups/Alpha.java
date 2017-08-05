package picard.cmdline.programgroups;

import picard.cmdline.CommandLineProgramGroup;

/**
 * @author ebanks
 */
public class Alpha implements CommandLineProgramGroup {
    @Override
    public String getName() { return "Alpha Tools"; }
    @Override
    public String getDescription() { return "Tools that are currently UNSUPPORTED until further testing and maturation."; }
}
