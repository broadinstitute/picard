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
package org.broad.tribble.bed;

import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.annotation.Strand;
import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.util.ParsingUtils;

import java.util.regex.Pattern;

/**
 * Codec for parsing BED file, as described by UCSC
 * See https://genome.ucsc.edu/FAQ/FAQformat.html#format1
 *
 * @author jrobinso
 *         Date: Dec 20, 2009
 */
public class BEDCodec extends AsciiFeatureCodec<BEDFeature> {

    private static final Pattern SPLIT_PATTERN = Pattern.compile("\\t|( +)");
    private final int startOffsetValue;

    /**
     * Calls {@link #BEDCodec(StartOffset)} with an argument
     * of {@code StartOffset.ONE}
     */
    public BEDCodec() {
        this(StartOffset.ONE);
    }


    /**
     * BED format is 0-based, but Tribble is 1-based.
     * Set desired start position at either ZERO or ONE
     */
    public BEDCodec(final StartOffset startOffset) {
        super(BEDFeature.class);
        this.startOffsetValue = startOffset.value();
    }


    public BEDFeature decodeLoc(String line) {
        return decode(line);
    }

    @Override
    public BEDFeature decode(String line) {

        if (line.trim().length() == 0) {
            return null;
        }

        if (line.startsWith("#") || line.startsWith("track") || line.startsWith("browser")) {
            this.readHeaderLine(line);
            return null;
        }

        String[] tokens = SPLIT_PATTERN.split(line, -1);
        return decode(tokens);
    }

    @Override
    public Object readActualHeader(LineIterator reader) {
        return null;
    }

    public BEDFeature decode(String[] tokens) {
        int tokenCount = tokens.length;

        // The first 3 columns are non optional for BED.  We will relax this
        // and only require 2.

        if (tokenCount < 2) {
            return null;
        }

        String chr = tokens[0];

        // The BED format uses a first-base-is-zero convention,  Tribble features use 1 => add 1.
        int start = Integer.parseInt(tokens[1]) + startOffsetValue;

        int end = start;
        if (tokenCount > 2) {
            end = Integer.parseInt(tokens[2]);
        }

        FullBEDFeature feature = new FullBEDFeature(chr, start, end);

        // The rest of the columns are optional.  Stop parsing upon encountering
        // a non-expected value

        // Name
        if (tokenCount > 3) {
            String name = tokens[3].replaceAll("\"", "");
            feature.setName(name);
        }

        // Score
        if (tokenCount > 4) {
            try {
                float score = Float.parseFloat(tokens[4]);
                feature.setScore(score);
            } catch (NumberFormatException numberFormatException) {

                // Unexpected, but does not invalidate the previous values.
                // Stop parsing the line here but keep the feature
                // Don't log, would just slow parsing down.
                return feature;
            }
        }

        // Strand
        if (tokenCount > 5) {
            String strandString = tokens[5].trim();
            char strand = (strandString.length() == 0)
                    ? ' ' : strandString.charAt(0);

            if (strand == '-') {
                feature.setStrand(Strand.NEGATIVE);
            } else if (strand == '+') {
                feature.setStrand(Strand.POSITIVE);
            } else {
                feature.setStrand(Strand.NONE);
            }
        }

        //Color
        if (tokenCount > 8) {
            String colorString = tokens[8];
            feature.setColor(ParsingUtils.parseColor(colorString));
        }

        // Coding information is optional
        if (tokenCount > 11) {
            createExons(start, tokens, feature, feature.getStrand());
        }

        return feature;
    }

    protected boolean readHeaderLine(String line) {
        //We don't parse BED header
        return false;
    }

    private void createExons(int start, String[] tokens, FullBEDFeature gene,
                             Strand strand) throws NumberFormatException {

        int cdStart = Integer.parseInt(tokens[6]) + startOffsetValue;
        int cdEnd = Integer.parseInt(tokens[7]);

        int exonCount = Integer.parseInt(tokens[9]);
        String[] exonSizes = new String[exonCount];
        String[] startsBuffer = new String[exonCount];
        ParsingUtils.split(tokens[10], exonSizes, ',');
        ParsingUtils.split(tokens[11], startsBuffer, ',');

        int exonNumber = (strand == Strand.NEGATIVE ? exonCount : 1);

        if (startsBuffer.length == exonSizes.length) {
            for (int i = 0; i < startsBuffer.length; i++) {
                int exonStart = start + Integer.parseInt(startsBuffer[i]);
                int exonEnd = exonStart + Integer.parseInt(exonSizes[i]) - 1;
                gene.addExon(exonStart, exonEnd, cdStart, cdEnd, exonNumber);

                if (strand == Strand.NEGATIVE) {
                    exonNumber--;
                } else {
                    exonNumber++;
                }
            }
        }
    }

    @Override
    public boolean canDecode(final String path) {
        return path.toLowerCase().endsWith(".bed");
    }

    public int getStartOffset() {
        return this.startOffsetValue;
    }

    /**
     * Indicate whether co-ordinates or 0-based or 1-based.
     * <p/>
     * Tribble uses 1-based, BED files use 0.
     * e.g.:
     * start_position = bedline_start_position - startIndex.value()
     */
    public enum StartOffset {
        ZERO(0),
        ONE(1);
        private int start;

        private StartOffset(int start) {
            this.start = start;
        }

        public int value() {
            return this.start;
        }
    }

}
