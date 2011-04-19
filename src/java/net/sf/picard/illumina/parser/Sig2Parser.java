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

import java.io.File;
import java.util.List;

/**
 * Parser for s_<lane>_<tile>_sig2.txt files, produced by Illumina 1.1 and 1.3
 * 
 * @author alecw@broadinstitute.org
 */
public class Sig2Parser extends IntensitiesOrNoiseParser{
    public Sig2Parser(final ReadConfiguration readConfiguration, final File directory, final int lane,
                      final List<Integer> tiles) {
        super(readConfiguration, directory, lane, "sig2", tiles);
    }

    /**
     * Stuff the FourChannelIntensityData into the processedIntensities property.
     *
     * @param end Object on which to store the data.
     * @param values Data to store.
     */
    @Override
    protected void setValues(final IlluminaEndData end, final FourChannelIntensityData values) {
        end.setProcessedIntensities(values);
    }
}
