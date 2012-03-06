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
 * There is one IlluminaData sub-interface for each IlluminaDataType enum value.
 * IlluminaParsers must return objects implementing at least one of the interfaces below.
 * IlluminaDataProvider will take IlluminaData objects created by IlluminaParsers and cast them to the types they
 * implement and these objects will then be used to populate the ClusterData object.
 *
 * @author jburke@broadinstitute.org
 */
interface IlluminaData {
}

// Note: PositionalData was spun out this round but since every parser has means of retrieving lane/tile from the
// file name, we are going to move lane/tile to be queryable from parsers in future revisions and therefore if you
// want lane/tile info you will NOT have to parse one of the Positional Data formats (pos, locs, clocs, qseqs)
interface PositionalData extends IlluminaData {
    public int getXCoordinate();
    public int getYCoordinate();
}

interface BaseData extends IlluminaData {
    public byte [][] getBases();
}

interface QualityData extends IlluminaData {
    public byte [][] getQualities();
}

interface NoiseData extends IlluminaData {
    public FourChannelIntensityData [] getNoise();
}

interface RawIntensityData extends IlluminaData{
    public FourChannelIntensityData [] getRawIntensities();
}

interface PfData extends IlluminaData {
    public boolean isPf();
}

interface BarcodeData extends IlluminaData {
    public String getBarcode();
}


