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
 * CnfParser takes a directory, lane, a map of tiles to Cycled file iterators, and a list of desired lengths for the
 * output FourChannelIntensityData and allows iteration over the clusters of all the provided tiles in the given lane.
 *
 * Note: Files passed by CycledIlluminaFileMap are not checked for proper extension (e.g. cif or cnf) so a CnfParser
 * can read a map to cif files and put it in the a NoiseData, you've been warned!
 * 
 * @author Jonathan Burke
 */
class CnfParser extends IlluminaIntensityParser<NoiseData> {
    private static final Set<IlluminaDataType> SupportedTypes = Collections.unmodifiableSet(CollectionUtil.makeSet(IlluminaDataType.Noise));

    public CnfParser(final File directory, final int lane, final CycleIlluminaFileMap tilesToCycleFiles, int[] outputLengths) {
        super(directory, lane, tilesToCycleFiles, outputLengths);
    }

    @Override
    protected void addIntensityToIlluminaData(NoiseData illData, final CompositeIndex index, IntensityChannel channel, short intensity) {
        illData.getNoise()[index.arrayIndex].getChannel(channel)[index.elementIndex ] = intensity;
    }

    @Override
    protected NoiseData intensityToIlluminaData(final FourChannelIntensityData[] fcids) {
        return new NoiseData() {
            @Override
            public FourChannelIntensityData[] getNoise() {
                return fcids;
            }
        };
    }

    @Override
    public Set<IlluminaDataType> supportedTypes() {
        return SupportedTypes;
    }
}
