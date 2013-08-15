/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
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
package org.broad.tribble.gelitext;

import net.sf.samtools.util.CollectionUtil;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.exception.CodecLineParsingException;
import org.broad.tribble.readers.LineIterator;

import java.util.Arrays;


/**
 * <p/>
 * A codec for parsing geli text files, which is the text version of the geli binary format.
 * <p/>
 * <p/>
 * GELI text has the following tab-seperated fields:
 * contig             the contig (string)
 * position           the position on the contig (long)
 * refBase            the reference base (char)
 * depthOfCoverage    the depth of coverage at this position (int)
 * maximumMappingQual the maximum mapping quality of a read at this position (int)
 * genotype           the called genotype (string)
 * LODBestToReference the LOD score of the best to the reference (double)
 * LODBestToNext      the LOD score of the best to the next best genotype (double)
 * likelihoods        the array of all genotype likelihoods, in ordinal ordering (array of 10 doubles, in ordinal order)
 *
 * @author aaron
 */
public class GeliTextCodec extends AsciiFeatureCodec<GeliTextFeature> {
    public GeliTextCodec() {
        super(GeliTextFeature.class);
    }

    public Feature decodeLoc(final String line) {
        return decode(line);
    }

    @Override
    public GeliTextFeature decode(final String line) {
        // clean out header lines and comments
        if (line.startsWith("#") || line.startsWith("@"))
            return null;

        // parse into tokens
        final String[] parts = line.trim().split("\\s+");
        return decode(parts);
    }

    @Override
    public Object readActualHeader(LineIterator reader) {
        return null;
    }

    public GeliTextFeature decode(final String[] tokens) {
        try {
            // check that we got the correct number of tokens in the split
            if (tokens.length != 18)
                throw new CodecLineParsingException("Invalid GeliTextFeature row found -- incorrect element count.  Expected 18, got " + tokens.length + " line = " + CollectionUtil.join(Arrays.asList(tokens), " "));

            // UPPER case and sort
            final char[] x = tokens[5].toUpperCase().toCharArray();
            Arrays.sort(x);
            final String bestGenotype = new String(x);

            final double[] genotypeLikelihoods = new double[10];
            for (int pieceIndex = 8, offset = 0; pieceIndex < 18; pieceIndex++, offset++) {
                genotypeLikelihoods[offset] = Double.valueOf(tokens[pieceIndex]);
            }
            return new GeliTextFeature(tokens[0],
                    Long.valueOf(tokens[1]),
                    Character.toUpperCase(tokens[2].charAt(0)),
                    Integer.valueOf(tokens[3]),
                    Integer.valueOf(tokens[4]),
                    DiploidGenotype.toDiploidGenotype(bestGenotype),
                    Double.valueOf(tokens[6]),
                    Double.valueOf(tokens[7]),
                    genotypeLikelihoods);
        } catch (CodecLineParsingException e) {
            e.printStackTrace();
            throw new RuntimeException("Unable to parse line " + CollectionUtil.join(Arrays.asList(tokens), " "), e);
        } catch (NumberFormatException e) {
            e.printStackTrace();
            throw new RuntimeException("Unable to parse line " + CollectionUtil.join(Arrays.asList(tokens), " "), e);
        }
    }
}
