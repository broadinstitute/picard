/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package net.sf.picard.illumina.parser;

import net.sf.picard.PicardException;

/**
 * Store the data read from Illumina files for a single read, separated into two ends, and barcode, if appropriate.
 * Caller must call setFirstEnd, setSecondEnd, setBarcodeEnd as appropriate to initialize those properties.
 * 
 * @author alecw@broadinstitute.org
 */
public class IlluminaReadData {

    private int lane = -1;
    private int tile = -1;
    private int x = -1;
    private int y = -1;
    private IlluminaEndData firstEnd;
    private IlluminaEndData secondEnd;
    private IlluminaEndData barcodeRead;
    private Boolean pf;
    private String matchedBarcode;

    public String toString() {
        return "IlluminaReadData(lane: " + lane + "; tile: " + tile + "; x: " + x + "; y: " + y + "; pf: " + pf +
                "; matchedBarcode: " + matchedBarcode + ")";
    }

    public int getTile() {
        return tile;
    }

    public void setTile(final int tile) {
        this.tile = tile;
    }

    public boolean tileIsSet() {
        return tile != -1;
    }

    /**
     * Either set this value if not already set, or if already set, throw an exception if new value != current value.
     */
    public void setOrCheckTile(final int tile) {
        if (tileIsSet()) {
            if (this.tile != tile) {
                throw new PicardException("Tile number mismatch for " + this + " : " + this.tile + " != " + tile);
            }
        } else {
            this.tile = tile;
        }
    }

    public int getLane() {
        return lane;
    }

    public void setLane(final int lane) {
        this.lane = lane;
    }

    public boolean laneIsSet() {
        return lane != -1;
    }

    /**
     * Either set this value if not already set, or if already set, throw an exception if new value != current value.
     */
    public void setOrCheckLane(final int lane) {
        if (laneIsSet()) {
            if (this.lane != lane) {
                throw new PicardException("Lane number mismatch for " + this + " : " + this.lane + " != " + lane);
            }
        } else {
            this.lane = lane;
        }
    }

    public int getX() {
        return x;
    }

    public void setX(final int x) {
        this.x = x;
    }

    public boolean xIsSet() {
        return x != -1;
    }

    /**
     * Either set this value if not already set, or if already set, throw an exception if new value != current value.
     */
    public void setOrCheckX(final int x) {
        if (xIsSet()) {
            if (this.x != x) {
                throw new PicardException("X value mismatch for " + this + " : " + this.x + " != " + x);
            }
        } else {
            this.x = x;
        }
    }

    public int getY() {
        return y;
    }

    public void setY(final int y) {
        this.y = y;
    }

    public boolean yIsSet() {
        return y != -1;
    }

    /**
     * Either set this value if not already set, or if already set, throw an exception if new value != current value.
     */
    public void setOrCheckY(final int y) {
        if (yIsSet()) {
            if (this.y != y) {
                throw new PicardException("Y value mismatch for " + this + " : " + this.y + " != " + y);
            }
        } else {
            this.y = y;
        }
    }

    public Boolean isPf() {
        return pf;
    }

    public void setPf(final boolean pf) {
        this.pf = pf;
    }

    /**
     * Either set this value if not already set, or if already set, throw an exception if new value != current value.
     */
    public void setOrCheckPf(final boolean pf) {
        if (this.pf == null) {
            this.pf = pf;
        } else if (this.pf != pf) {
            throw new PicardException("PF value mismatch for " + this + " : ");
        }
    }

    public boolean isPairedEnd() {
        return secondEnd != null;
    }

    /**
     * @return 2 if paired-end, else 1.
     */
    public int getNumEnds() {
        return (isPairedEnd()? 2: 1);
    }
    
    /**
     * @return second end if oneBasedIndex == 2, else first end.
     */
    public IlluminaEndData getEnd(final int oneBasedEndIndex) {
        return (oneBasedEndIndex == 2? secondEnd: firstEnd);
    }

    public void setEnd(final IlluminaEndData end, final int oneBasedEndIndex) {
        if (oneBasedEndIndex == 1) {
            firstEnd = end;
        } else if (oneBasedEndIndex == 2) {
            secondEnd = end;
        } else {
            throw new IllegalArgumentException();
        }
    }

    public IlluminaEndData getFirstEnd() {
        return firstEnd;
    }

    public void setFirstEnd(final IlluminaEndData firstEnd) {
        this.firstEnd = firstEnd;
    }

    public IlluminaEndData getSecondEnd() {
        return secondEnd;
    }

    public void setSecondEnd(final IlluminaEndData secondEnd) {
        this.secondEnd = secondEnd;
    }

    public IlluminaEndData getBarcodeRead() {
        return barcodeRead;
    }

    public void setBarcodeRead(final IlluminaEndData barcodeRead) {
        this.barcodeRead = barcodeRead;
    }

    /**
     * @return The barcode matched (not the actual sequence from the read, which may not perfectly match
     * the barcode).
     */
    public String getMatchedBarcode() {
        return matchedBarcode;
    }

    public void setMatchedBarcode(final String matchedBarcode) {
        this.matchedBarcode = matchedBarcode;
    }

    public IlluminaEndData getEnd(final EndType whichEndType) {
        switch (whichEndType) {
            case BARCODE:
                return barcodeRead;
            case FIRST:
                return firstEnd;
            case SECOND:
                return secondEnd;
            default:
                throw new IllegalArgumentException("Null or strange value passed to getEnd");
        }
    }
}
