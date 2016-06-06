package picard.sam.util;

/**
 * Small interface that provides access to the read name information.
 */
/*public interface ReadNameInterface {
    public static int NO_VALUE = -1;

    public String readname = null;

    public int read1ReferenceIndex = -1;


    public String getReadName();
}*/

public class ReadNameInterface {

    //public static int NO_VALUE = -1;

    public int read1IndexInFile = -1;

    public int setSize = -1;

    public String readname = null;

    public String getReadName() {return readname;};

    public void setReadName(final String name) { this.readname = name; }

}