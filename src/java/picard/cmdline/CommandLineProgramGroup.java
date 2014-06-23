package picard.cmdline;

import java.util.Comparator;

/**
 * Interface for groups of CommandLinePrograms.
 * @author Nils Homer
 */
public interface CommandLineProgramGroup {

    /** Gets the name of this program. **/
    public String getName();
    /** Gets the description of this program. **/
    public String getDescription();
    /** Compares two program groups by name. **/
    static public Comparator<Class> comparator = new Comparator<Class>() {
        public int compare(Class a, Class b) {
            return a.getName().compareTo(b.getName());
        }
    };
}
