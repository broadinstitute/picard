package picard.sam.util;

import picard.PicardException;

/**
 * Small class that provides access to the physical location information about a cluster.
 * All values should be defaulted to -1 if unavailable.  Tile should only allow
 * non-zero positive integers, x and y coordinates must be non-negative.
 * This is different from PhysicalLocationShort in that the x and y positions are ints, not shorts
 * thus, they do not overflow within a HiSeqX tile.
 */
public class PhysicalLocationInt implements PhysicalLocation {

    public short tile = -1;
    public int x = -1, y = -1;

    public short getReadGroup() { throw new PicardException("Not Implemented"); }

    public void setReadGroup(final short readGroup) { throw new PicardException("Not Implemented"); }

    public short getTile() { return tile; }

    public void setTile(final short tile) { this.tile = tile; }

    public int getX() { return x; }

    public void setX(final int x) { this.x = x; }

    public int getY() { return y; }

    public void setY(final int y) { this.y = y; }

    public short getLibraryId() { throw new PicardException("Not Implemented"); }

    public void setLibraryId(final short libraryId) { throw new PicardException("Not Implemented"); }

}
