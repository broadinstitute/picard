package picard.cmdline;

/**
 * Interface for groups of CommandLinePrograms.
 * @author Nils Homer
 */
public interface CommandLineProgramGroup {

    /** Gets the name of this program. **/
    public String getName();
    /** Gets the description of this program. **/
    public String getDescription();
}
