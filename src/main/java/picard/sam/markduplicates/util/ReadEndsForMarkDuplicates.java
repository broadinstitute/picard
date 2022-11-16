/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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

package picard.sam.markduplicates.util;

import htsjdk.samtools.SAMRecord;
import picard.sam.markduplicates.MarkDuplicatesForFlowHelper;

/**
 * Little struct-like class to hold read pair (and fragment) end data for MarkDuplicatesWithMateCigar
 *
 * @author Nils Homer
 */
public class ReadEndsForMarkDuplicates extends ReadEnds implements Cloneable {
    /*
    What do we need to store you ask?  Well, we need to store:
       - byte: orientation
       - short: libraryId, readGroup, tile, x, y, score
       - int: read1ReferenceIndex, read1Coordinate, read2ReferenceIndex, read2Coordinate, duplicateSetSize
       - long: read1IndexInFile, read2IndexInFile
     */
    protected static final int SIZE_OF = (1 * 1) + (5 * 2) + (5 * 4) + (8 * 2) + 1
            + 8 + // last 8 == reference overhead
            13; // This is determined experimentally with JProfiler

    public static int getSizeOf() {
        return SIZE_OF;
    }

    public short score = 0;
    public long read1IndexInFile = -1;
    public long read2IndexInFile = -1;
    public int duplicateSetSize = -1;

    public ReadEndsForMarkDuplicates() {}

    public ReadEndsForMarkDuplicates(final ReadEndsForMarkDuplicates read) {
        this.libraryId = read.getLibraryId();
        this.orientation = read.orientation;
        this.read1ReferenceIndex = read.read1ReferenceIndex;
        this.read1Coordinate = read.read1Coordinate;
        this.read2ReferenceIndex = read.read2ReferenceIndex;
        this.read2Coordinate = read.read2Coordinate;

        this.readGroup = read.getReadGroup();
        this.tile = read.getTile();
        this.x = read.x;
        this.y = read.y;

        this.orientationForOpticalDuplicates = read.orientationForOpticalDuplicates;

        this.score = read.score;

        this.read1IndexInFile = read.read1IndexInFile;
        this.read2IndexInFile = read.read2IndexInFile;
    }

    @Override
    public String toString() {
        return String.format("%d %d %d", read1IndexInFile, read1Coordinate, score);
    }

    @Override
    public ReadEndsForMarkDuplicates clone() {
        return new ReadEndsForMarkDuplicates(this);
    }

    /**
     * This method is used to generate the following two metrics:
     * UNPAIRED_DUPS_WITH_TLEN
     * UNPAIRED_DUPS_WITHOUT_TLEN
     *
     * It will return true if and only if the read is single ended and the exact fragment length is
     *  known (i.e. it was not quality trimmed)
     */
    public static boolean isSingleEndReadKnownFragment(final SAMRecord rec) {
        if ( rec.getReadUnmappedFlag() || rec.getReadPairedFlag() ) {
            return false;
        } else if ( MarkDuplicatesForFlowHelper.isAdapterClipped(rec) ) {
            return true;
        } else if ( !rec.getReadNegativeStrandFlag() ) {
            return rec.getEnd() != rec.getUnclippedEnd();
        } else {
            return rec.getStart() != rec.getUnclippedStart();
        }
    }

}