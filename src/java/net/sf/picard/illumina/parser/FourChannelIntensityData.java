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

import java.util.Arrays;

/**
 * Holds a 4 short values for each cycle of a read.  This is used, to store raw intensities,
 * processed intensities, or noise.  Note that for Illumina 1.1 and 1.3, these are floating point values,
 * but are truncated to shorts to store here.
 *
 * Indices into the channel arrays are zero-based, i.e. the first cycle is 0.
 *
 * @author jburke@broadinstitute.org
 */
public class FourChannelIntensityData {
    /**
     * Major index: channel number; minor index: cycle number (zero based)
     */
    private short [] a;
    private short [] c;
    private short [] g;
    private short [] t;

    public FourChannelIntensityData(final int numberOfCycles) {
        a = new short[numberOfCycles];
        c = new short[numberOfCycles];
        g = new short[numberOfCycles];
        t = new short[numberOfCycles];
    }

    public short[] getChannel(final IntensityChannel channel) {
        switch(channel) {
            case A: return a;
            case C: return c;
            case G: return g;
            case T: return t;
        }

        throw new PicardException("Unexpected intensity channel " + channel);
    }

    public short[] getA() {
        return a;
    }

    public short[] getC() {
        return c;
    }

    public short[] getG() {
        return g;
    }

    public short[] getT() {
        return t;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final FourChannelIntensityData that = (FourChannelIntensityData)o;
        return Arrays.equals(this.a, that.a) &&
               Arrays.equals(this.c, that.c) &&
               Arrays.equals(this.g, that.g) &&
               Arrays.equals(this.t, that.t);
    }

    @Override
    public int hashCode() {
        int ret = 0; //TODO: Check this
        ret = ret * 31 + Arrays.hashCode(a);
        ret += ret * 31 + Arrays.hashCode(c);
        ret += ret * 31 + Arrays.hashCode(g);
        ret += ret * 31 + Arrays.hashCode(t);
        return ret;
    }
}
