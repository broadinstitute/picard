package picard.cmdline;

import java.io.File;

/**
 * Embodies defaults for global values that affect how the Picard Command Line operates. Defaults are encoded in the class
 * and are also overridable using system properties.
 *
 * @author Nils Homer
 */
public class CommandLineDefaults {

    /** Implementation note, this is duplicate code stolen from HTSJDK's Default.java

    /**
     * Decides if we want to write colors to the terminal.
     */
    public static final boolean COLOR_STATUS;

    static {
        COLOR_STATUS = getBooleanProperty("color_status", true);
    }

    /** Gets a string system property, prefixed with "picard.cmdline." using the default if the property does not exist. */
    private static String getStringProperty(final String name, final String def) {
        return System.getProperty("picard.cmdline." + name, def);
    }

    /** Gets a boolean system property, prefixed with "picard.cmdline." using the default if the property does not exist. */
    private static boolean getBooleanProperty(final String name, final boolean def) {
        final String value = getStringProperty(name, new Boolean(def).toString());
        return Boolean.parseBoolean(value);
    }

    /** Gets an int system property, prefixed with "picard.cmdline." using the default if the property does not exist. */
    private static int getIntProperty(final String name, final int def) {
        final String value = getStringProperty(name, new Integer(def).toString());
        return Integer.parseInt(value);
    }

    /** Gets a File system property, prefixed with "picard.cmdline." using the default if the property does not exist. */
    private static File getFileProperty(final String name, final String def) {
        final String value = getStringProperty(name, def);
        // TODO: assert that it is readable
        return (null == value) ? null : new File(value);
    }
}

