/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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

import java.io.File;
import java.util.List;

/**
 * Parser for Illumina RTA CNF noise files.  These files are produced per cycle, with all the clusters
 * in a tile in a single file, so in order to get the values for a read, multiple files must be accessed.
 * 
 * @author alecw@broadinstitute.org
 */
public class CnfParser extends IlluminaCycleFileSetParser {

    public CnfParser(final ReadConfiguration readConfiguration, final File directory, final int lane,
                     final List<Integer> tiles) {
        super(readConfiguration, directory, lane, ClusterIntensityFileReader.FileType.cnf, tiles);
    }

    /**
     * Tell abstract base class which slot is being parsed.
     * @param readData Object into which to stuff the FourChannelIntensityData.
     * @param fourChannelIntensityData thing to stuff into the right slot of end.
     */
    protected void setFCID(final ReadData readData, final FourChannelIntensityData fourChannelIntensityData) {
        readData.setNoise(fourChannelIntensityData);
    }

    /**
     * @param readData Object from which to get the FourChannelIntensityData.
     * @return FourChannelIntensityData from the slot being parsed by this concrete class.
     */
    protected FourChannelIntensityData getFCID(final ReadData readData) {
        return readData.getNoise();
    }
}
