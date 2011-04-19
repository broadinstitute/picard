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

/**
 * Data for a single end of a paired-end read, a barcode read, or for the entire read if not paired end.
 *
 * @author alecw@broadinstitute.org
 */
public class IlluminaEndData {
    private byte[] bases;
    private byte[] qualities;
    private FourChannelIntensityData rawIntensities;
    private FourChannelIntensityData noise;
    private FourChannelIntensityData processedIntensities;

    /**
     * @return ASCII byte representation of bases.
     */
    public byte[] getBases() {
        return bases;
    }

    public void setBases(final byte[] bases) {
        this.bases = bases;
    }

    /**
     * @return Noise values as produced by Illumina software, converted to shorts.
     */
    public FourChannelIntensityData getNoise() {
        return noise;
    }

    public void setNoise(final FourChannelIntensityData noise) {
        this.noise = noise;
    }

    /**
     * @return Processed intensity values as produced by Illumina software, converted to shorts.
     */
    public FourChannelIntensityData getProcessedIntensities() {
        return processedIntensities;
    }

    public void setProcessedIntensities(final FourChannelIntensityData processedIntensities) {
        this.processedIntensities = processedIntensities;
    }

    /**
     * @return Phred-binary scaled qualities.  E.g. Q20 is the byte with value==20.
     */
    public byte[] getQualities() {
        return qualities;
    }

    public void setQualities(final byte[] qualities) {
        this.qualities = qualities;
    }

    /**
     * @return Raw intensity values as produced by Illumina software, converted to shorts.
     */
    public FourChannelIntensityData getRawIntensities() {
        return rawIntensities;
    }

    public void setRawIntensities(final FourChannelIntensityData rawIntensities) {
        this.rawIntensities = rawIntensities;
    }

}
