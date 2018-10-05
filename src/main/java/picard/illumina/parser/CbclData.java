package picard.illumina.parser;

import picard.illumina.parser.readers.AbstractIlluminaPositionFileReader;

/**
 * This class provides that data structure for cbcls. This includes BCL data as well as PF (pass-filter) data and
 * positional information.
 */
public class CbclData extends BclData implements PfData, PositionalData {
    private final int tile;
    private AbstractIlluminaPositionFileReader.PositionInfo positionInfo;
    private boolean pfRead;


    public CbclData(int[] outputLengths, int tile) {
        super(outputLengths);
        this.tile = tile;
    }

    public int getTile() {
        return tile;
    }

    public void setPositionInfo(AbstractIlluminaPositionFileReader.PositionInfo positionInfo) {
        this.positionInfo = positionInfo;
    }

    public void setPfRead(boolean pfRead) {
        this.pfRead = pfRead;
    }


    @Override
    public int getXCoordinate() {
        return this.positionInfo.xQseqCoord;
    }

    @Override
    public int getYCoordinate() {
        return this.positionInfo.yQseqCoord;
    }

    @Override
    public boolean isPf() {
        return pfRead;
    }
}
