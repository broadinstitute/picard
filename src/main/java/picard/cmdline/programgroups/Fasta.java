package picard.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;

/**
* @author nhomer
*/
public class Fasta implements CommandLineProgramGroup {
    @Override
    public String getName() { return "Fasta"; }
    @Override
    public String getDescription() { return "Tools for manipulating FASTA, or related data."; }
}
