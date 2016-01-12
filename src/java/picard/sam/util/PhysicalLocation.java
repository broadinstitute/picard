package picard.sam.util;

/**
 * Small interface that provides access to the physical location information about a cluster.
 * All values should be defaulted to -1 if unavailable.  ReadGroup and Tile should only allow
 * non-zero positive integers, x and y coordinates may be negative.
 */
public interface PhysicalLocation {
    public int NO_VALUE = -1;

    public short getReadGroup();

    public void setReadGroup(short rg);

    public short getTile();

    public void setTile(short tile);

    public int getX();

    public void setX(int x);

    public int getY();

    public void setY(int y);

    public short getLibraryId();

    public void setLibraryId(short libraryId);

    /** Default implementation of a method to check whether real location data has been set. */
    default public boolean hasLocation() {
        return getTile() != NO_VALUE;
    }
}
