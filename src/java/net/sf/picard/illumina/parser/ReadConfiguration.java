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
import net.sf.samtools.util.CoordMath;

import java.util.List;
import java.util.ArrayList;

/**
 * Simple struct for holding the configuration of a run -- single end or paired end, length of the
 * ends, where the barcode is.  All indices are 1-based, inclusive.
 * 
 * @author alecw@broadinstitute.org
 */
public class ReadConfiguration {
    private boolean pairedEnd;
    private boolean barcoded;
    private ReadType barcodeRead;
    private int barcodeReads;

    private final InclusiveRange first =  new InclusiveRange("first");
    // These are not always used but it's easier to initialize them unconditionally.  There
    // are not many of these object so it shouldn't be a problem.
    private final InclusiveRange second = new InclusiveRange("second");
    private final InclusiveRange barcode = new InclusiveRange("barcode");


    /**
     * @throws PicardException if the read configuration doesn't make sense, e.g. ranges are not
     * contiguous and starting at 1.
     */
    public void assertValid() { assertValid(false); }


    /**
     * @throws PicardException if the read configuration doesn't make sense, e.g. ranges are not
     * contiguous and starting at 1.
     * @param allowZeroLengthFirstEnd If true, consider a read configuration with empty first end to be valid.
     * 
     */
    public void assertValid(boolean allowZeroLengthFirstEnd) {
        if (first.getLength() > 0 || ! allowZeroLengthFirstEnd) {
            first.assertValid();
        }
        if (pairedEnd) {
            second.assertValid();
            if (first.compareTo(second) >= 0) {
                throw new PicardException("Second end precedes first end in ReadConfiguration.");
            }
        }
        if (barcoded) {
            barcode.assertValid();
        }
        // For checking contiguity
        final List<InclusiveRange> ranges = new ArrayList<InclusiveRange>();
        ranges.add(first);
        if (pairedEnd) {
            ranges.add(second);
        }
        if (barcoded) {
            boolean added = false;
            for (int i = 0; i < ranges.size(); ++i) {
                if (barcode.compareTo(ranges.get(i)) < 0) {
                    ranges.add(i, barcode);
                    added = true;
                    break;
                }
            }
            if (!added) {
                ranges.add(barcode);
            }
        }
        // Test for contiguity, and starting at cycle 1, and second > first
        if (ranges.get(0).start != 1) {
            throw new PicardException("ReadConfiguration does not start at cycle 1.");
        }
        for (int i = 0; i < ranges.size() - 1; ++i) {
            if (!ranges.get(i).abuts(ranges.get(i+1))) {
                throw new PicardException("ReadConfiguration has cycle number gap.");
            }
            if (ranges.get(i).compareTo(ranges.get(i+1)) >= 0) {
                throw new PicardException("That's unpossible");
            }
        }
    }

    /**
     * @return 1-based highest cycle number in this ReadConfiguration.
     */
    public int getMaxCycleNumber() {
        final int ret;
        if (pairedEnd) {
            ret = second.end;
        } else {
            ret = first.end;
        }
        if (barcoded) {
            return Math.max(ret, barcode.end);
        }
        return ret;
    }

    public boolean isBarcoded() {
        return barcoded;
    }

    public void setBarcoded(final boolean barcoded) {
        this.barcoded = barcoded;
    }

    public int getBarcodeReads() { return barcodeReads; }

    public void setBarcodeReads(final int barcodeReads) { this.barcodeReads = barcodeReads; }

    /**
     * Which end contains the barcode.
     */
    public ReadType getBarcodeRead() {
        return barcodeRead;
    }

    /**
     * Set the end that contains the barcode.
     */
    public void setBarcodeRead(final ReadType barcodeRead) {
        this.barcodeRead = barcodeRead;
    }

    /**
     * @return last one-based cycle number (inclusive) of the barcode.
     */
    public int getBarcodeEnd() {
        return barcode.end;
    }

    public void setBarcodeEnd(final int barcodeEnd) {
        this.barcode.end = barcodeEnd;
    }

    /**
     * @return first one-based cycle number (inclusive) of the barcode.
     */
    public int getBarcodeStart() {
        return barcode.start;
    }

    public void setBarcodeStart(final int barcodeStart) {
        this.barcode.start = barcodeStart;
    }

    /**
     * @return last one-based cycle number (inclusive) of the first end.
     */
    public int getFirstEnd() {
        return first.end;
    }

    public void setFirstEnd(final int firstEnd) {
        this.first.end = firstEnd;
    }

    /**
     * @return first one-based cycle number (inclusive) of the first end.
     */
    public int getFirstStart() {
        return first.start;
    }

    public void setFirstStart(final int firstStart) {
        this.first.start = firstStart;
    }

    public boolean isPairedEnd() {
        return pairedEnd;
    }

    public void setPairedEnd(final boolean pairedEnd) {
        this.pairedEnd = pairedEnd;
    }

    /**
     * @return last one-based cycle number (inclusive) of the second end.
     */
    public int getSecondEnd() {
        return second.end;
    }

    public void setSecondEnd(final int secondEnd) {
        this.second.end = secondEnd;
    }

    /**
     * @return first one-based cycle number (inclusive) of the first end.
     */
    public int getSecondStart() {
        return second.start;
    }

    public void setSecondStart(final int secondStart) {
        this.second.start = secondStart;
    }

    public int getFirstLength() {
        return first.getLength();
    }

    /**
     * @return length of second end, or 0 if not paired end.
     */
    public int getSecondLength() {
        if (!pairedEnd) {
            return 0;
        }
        return second.getLength();
    }

    /**
     * @return length of second end, or 0 if not barcoded.
     */
    public int getBarcodeLength() {
        if (!barcoded) {
            return 0;
        }
        return barcode.getLength();
    }

    public InclusiveRange getFirstRange() {
        return first;
    }

    public InclusiveRange getSecondRange() {
        return second;
    }

    public InclusiveRange getBarcodeRange() {
        return barcode;
    }

    /**
     * Get 0-based offset of start of barcode for splitting bases & quals out of a read.
     */
    public int getOffsetOfBarcodeInRead() {
        if (!barcoded) {
            throw new IllegalStateException();
        }
        final InclusiveRange range = (barcodeRead == ReadType.FIRST? first: second);
        if (getBarcodeEnd() + 1 == range.getStart()) {
            // Barcode is at start of read
            return 0;
        } else {
            // Barcode is at end of read
            return range.getLength();
        }
    }

    /**
     * Get 0-based offset of start of non-barcode bases for splitting bases & quals out of a read.
     */
    public int getOffsetOfNonBarcodeInRead() {
        if (!barcoded) {
            throw new IllegalStateException();
        }
        final InclusiveRange range = (barcodeRead == ReadType.FIRST? first: second);
        if (getBarcodeEnd() + 1 == range.getStart()) {
            // Barcode is at start of read
            return getBarcodeLength();
        } else {
            // Barcode is at end of read
            return 0;
        }
    }

    /**
     * Determine the end that a cycle number between 1 and total read length corresponds to
     * @param cycle one-based cycle number.
     * @return Which end the cycle number corresponds to.
     */
    public ReadType getEndTypeForCycle(final int cycle) {
        if (cycle >= first.getStart() && cycle <= first.getEnd()) {
            return ReadType.FIRST;
        }
        if (cycle >= second.getStart() && cycle <= second.getEnd()) {
            return ReadType.SECOND;
        }
        if (cycle >= barcode.getStart() && cycle <= barcode.getEnd()) {
            return ReadType.BARCODE;
        }
        throw new IllegalArgumentException("Invalid cycle number for read configuration: " + cycle);
    }

    /**
     * Get zero-based offset into an IlluminaEndData for a cycle number.
     * @param cycle one-based cycle number.
     * @return Zero-based offset into first end, second end or barcode corresponding to the cycle number.
     */
    public int getOffsetForCycle(final int cycle) {
        if (cycle >= first.getStart() && cycle <= first.getEnd()) {
            return cycle - first.getStart();
        }
        if (cycle >= second.getStart() && cycle <= second.getEnd()) {
            return cycle - second.getStart();
        }
        if (cycle >= barcode.getStart() && cycle <= barcode.getEnd()) {
            return cycle - barcode.getStart();
        }
        throw new IllegalArgumentException("Invalid cycle number for read configuration: " + cycle);
    }

    public static class InclusiveRange {
        String name;
        int start = 0;
        int end = 0;

        InclusiveRange(final String name) {
            this.name = name;
        }

        public int getEnd() {
            return end;
        }

        public void setEnd(final int end) {
            this.end = end;
        }

        public int getStart() {
            return start;
        }

        public void setStart(final int start) {
            this.start = start;
        }

        /**
         * @return length, or 0 if range is not valid
         */
        public int getLength() {
            if (end == 0) {
                return 0;
            }
            return CoordMath.getLength(start, end);
        }

        public void assertValid() {
            if (start < 1) {
                throw new PicardException("ReadConfiguration.Range " + name + " has start < 1");
            }
            if (end < start) {
                throw new PicardException("ReadConfiguration.Range " + name + " has end < start");
            }
        }

        public boolean overlaps(final InclusiveRange other) {
            return CoordMath.overlaps(start, end, other.start, other.end);
        }

        public int compareTo(final InclusiveRange that) {
            if (this.start != that.start) {
                return this.start - that.start;
            }
            return this.end - that.end;
        }

        public boolean abuts(final InclusiveRange that) {
            return this.end + 1 == that.start || that.end + 1 == this.start;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final InclusiveRange that = (InclusiveRange) o;

            if (end != that.end) return false;
            if (start != that.start) return false;
            if (name != null ? !name.equals(that.name) : that.name != null) return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result = name != null ? name.hashCode() : 0;
            result = 31 * result + start;
            result = 31 * result + end;
            return result;
        }
    }
}
