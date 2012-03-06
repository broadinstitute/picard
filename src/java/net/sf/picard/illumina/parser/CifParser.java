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

import net.sf.picard.util.CollectionUtil;

import java.io.File;
import java.util.Collections;
import java.util.Set;

/**
 * CifParser takes a directory, lane, a map of tiles to Cycled file iterators, and a list of desired lengths for the
 * output FourChannelIntensityData and allows iteration over the clusters of all the provided tiles in the given lane.
 *
 * Note: Files passed by CycledIlluminaFileMap are not checked for proper extension (e.g. cif or cnf) so a CifParser
 * can read a map to cnf files and put it in a RawIntensityData, you've been warned!
 *
 * @author jburke@broadinstitute.org
 */
class CifParser extends IlluminaIntensityParser<RawIntensityData> {
    private static final Set<IlluminaDataType> SupportedTypes = Collections.unmodifiableSet(CollectionUtil.makeSet(IlluminaDataType.RawIntensities));
    /**
     * @param directory BaseCallsDir or analogue directory containing folders like L001/C1.1/s_1_1.cif
     * @param lane The lane being analyzed
     * @param tilesToCycleFiles A map of tile -> CycledFilesIterator for each tile to consider
     * by RawIntensityData
     */
    public CifParser(final File directory, final int lane, final CycleIlluminaFileMap tilesToCycleFiles, final OutputMapping outputMapping) {
        super(directory, lane, tilesToCycleFiles, outputMapping);
    }

    /**
     * Populate the RawIntensityData from a short value and input indices for a cycle
     * @param rawIntensityData RawIntensityData to populate
     * @param index And index to the correct FCID and the correct index in that FCID
     * @param channel A,C,G, or T channel to populate
     * @param intensity
     */
    @Override
    protected void addIntensityToIlluminaData(final RawIntensityData rawIntensityData, final OutputMapping.TwoDIndex index, final IntensityChannel channel, final short intensity) {
        rawIntensityData.getRawIntensities()[index.majorIndex].getChannel(channel)[index.minorIndex] = intensity;
    }

    @Override
    protected RawIntensityData intensityToIlluminaData(final FourChannelIntensityData[] fcids) {
        return new RawIntensityData() {
            @Override
            public FourChannelIntensityData[] getRawIntensities() {
                return fcids;
            }
        };
    }

    @Override
    public Set<IlluminaDataType> supportedTypes() {
        return SupportedTypes;
    }
}
