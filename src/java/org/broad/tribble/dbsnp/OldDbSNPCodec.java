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
package org.broad.tribble.dbsnp;

import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.annotation.Strand;


/**
 * @author aaron
 *
 * Example format:
 * 585 chr1 433 433 rs56289060  0  +  - - -/C  genomic  insertion unknown 0  0  unknown  between  1
 * 585 chr1 491 492 rs55998931  0  +  C C C/T  genomic  single   unknown 0 0 unknown exact 1
 */
public class OldDbSNPCodec extends AsciiFeatureCodec<OldDbSNPFeature> {

    // the number of tokens we expect to parse from a dbSNP line
    static final int expectedTokenCount = 18;

    public OldDbSNPCodec() {
        super(OldDbSNPFeature.class);
    }

    public Feature decodeLoc(String line) {
        return decode(line);
    }
    
    /**
     * Decode a line as a db SNP Feature.
     *
     * @param line the line to decode
     *
     * @return Return the Feature encoded by the line, or null if the line does not represent a feature (e.g. is
     *         a comment)
     *
     * The ordering of db SNP fields from the UCSC track browser:
     * 1 bin
     * 2 chromosome
     * 3 chromosome Start
     * 4 chromosome End
     * 5 name
     * 6 score
     * 7 strand
     * 8 reference base NCBI
     * 9 reference base UCSC
     * 10 observed base
     * 11 mol. Type
     * 12 class
     * 13 valid
     * 14 avHet
     * 15 avHetSE
     * 16 functions
     * 17 locType
     * 18 weight
     */
    public OldDbSNPFeature decode(String line) {

        // we may be asked to process a header line; ignore it
        if (line.startsWith("#")) return null;

        // split the line
        String[] tokens = line.split("\\t+");
        return decode(tokens);
    }

    @Override
    public OldDbSNPFeature decode(String[] tokens){
        // check to see if we've parsed the string into the right number of tokens (expectedTokenCount)
        if (tokens.length != expectedTokenCount)
            return null;
        //    throw new CodecLineParsingException("the dbSNP line didn't have the expected number of tokens " +
        //                                        "(expected = " + expectedTokenCount + ", saw = " + tokens.length + " on " +
        //
            //                                         "line = " + line + ")");
        // create a new feature from the line
        int start = Integer.valueOf(tokens[2])+1;
        int stop = Integer.valueOf(tokens[3]);
        stop = (stop < start) ? start : stop; // Indels can be of length zero in dbSNP, we make them length one
        OldDbSNPFeature feature = new OldDbSNPFeature(tokens[1],
						      start,
						      stop);

        feature.setRsID(tokens[4]);
        feature.setScore(Integer.valueOf(tokens[5]));
        feature.setStrand(tokens[6].equals("+") ? Strand.POSITIVE : Strand.NEGATIVE);
        feature.setNCBIRefBase(tokens[7]);
        feature.setUCSCRefBase(tokens[8]);
        // split the observed bases
        feature.setObserved(tokens[9].split("/"));        
        feature.setMolType(tokens[10]);
        feature.setVariantType(tokens[11]);
        feature.setValidationStatus(tokens[12]);
        feature.setAvHet(Double.valueOf(tokens[13]));
        feature.setAvHetSE(Double.valueOf(tokens[14]));
        feature.setFunction(tokens[15]);
        feature.setLocationType(tokens[16]);
        feature.setWeight(Integer.valueOf(tokens[17]));

        // return the setup feature                                                                                
        return feature;
    }
}
