package picard.cmdline.programgroups;

import picard.cmdline.CommandLineProgramGroup;

public class Array implements CommandLineProgramGroup {
    @Override
    public String getName() { return "Array Tools"; }

    @Override
    public String getDescription() { return "Tools for manipulating data specific to genotyping arrays"; }
}

