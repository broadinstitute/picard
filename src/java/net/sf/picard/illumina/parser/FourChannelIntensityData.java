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

import java.util.Arrays;

/**
 * Holds a 4 short values for each cycle of a read.  This is used, e.g. to store raw intensities,
 * processed intensities, or noise.  Note that for Illumina 1.1 and 1.3, these are floating point values,
 * but are truncated to shorts to store here.
 *
 * Indices into the channel arrays are zero-based, i.e. the first cycle is 0.
 *
 * @author alecw@broadinstitute.org
 */
public class FourChannelIntensityData {
    /**
     * Major index: channel number; minor index: cycle number (zero based)
     */
    private short[][] channels = new short[IntensityChannel.NUM_CHANNELS][];

    /**
     * This ctor does not allocate the arrays for each channel.  They must be set with
     * the setters below.
     */
    public FourChannelIntensityData() {
    }

    /**
     * This ctor allocates the array for each channel.
     * @param length how big to make the array for each channel.
     */
    public FourChannelIntensityData(final int length) {
        for (int i = 0; i < channels.length; ++i) {
            channels[i] = new short[length];
        }
    }

    public short[] getChannel(final IntensityChannel channel) {
        return channels[channel.ordinal()];
    }

    public void setChannel(final IntensityChannel channel, final short[] array) {
        channels[channel.ordinal()] = array;
    }

    public short[] getA() {
        return getChannel(IntensityChannel.A_CHANNEL);
    }

    public void setA(final short[] a) {
        setChannel(IntensityChannel.A_CHANNEL, a);
    }

    public short[] getC() {
        return getChannel(IntensityChannel.C_CHANNEL);
    }

    public void setC(final short[] c) {
        setChannel(IntensityChannel.C_CHANNEL, c);
    }

    public short[] getG() {
        return getChannel(IntensityChannel.G_CHANNEL);
    }

    public void setG(final short[] g) {
        setChannel(IntensityChannel.G_CHANNEL, g);
    }

    public short[] getT() {
        return getChannel(IntensityChannel.T_CHANNEL);
    }

    public void setT(final short[] t) {
        setChannel(IntensityChannel.T_CHANNEL, t);
    }

    /**
     * @return Major index: channel number; minor index: cycle number (zero based)
     */
    public short[][]getChannels() {
        return channels;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final FourChannelIntensityData other = (FourChannelIntensityData)o;

        for (int i = 0; i < channels.length; ++i) {
            if (!Arrays.equals(channels[i], other.channels[i])) {
                return false;
            }
        }
        return true;
    }

    @Override
    public int hashCode() {
        int ret = 0;
        for (final short[] channel : channels) {
            ret = ret * 31 + Arrays.hashCode(channel);
        }
        return ret;
    }
}
