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
package net.sf.samtools;

import net.sf.samtools.util.StringUtil;

/**
 * Convert between String and Cigar class representations of CIGAR.
 */
public class TextCigarCodec
{
    private static final byte ZERO_BYTE = "0".getBytes()[0];
    private static final byte NINE_BYTE = "9".getBytes()[0];
    
    private static final TextCigarCodec singleton = new TextCigarCodec();

    /**
     * It is not necessary to get the singleton but it is preferable to use the same one
     * over and over vs. creating a new object for each BAMRecord.  There is no state in this
     * class so this is thread-safe.
     * @return A singleton TextCigarCodec useful for converting Cigar classes to and from strings
     */
    public static TextCigarCodec getSingleton() {
        return singleton;
    }


    /**
     * Convert from Cigar class representation to String.
     * @param cigar in Cigar class format
     * @return CIGAR in String form ala SAM text file.  "*" means empty CIGAR.
     */
    public String encode(final Cigar cigar) {
        if (cigar.isEmpty()) {
            return SAMRecord.NO_ALIGNMENT_CIGAR;
        }
        final StringBuilder ret = new StringBuilder();
        for (final CigarElement cigarElement : cigar.getCigarElements()) {
            ret.append(cigarElement.getLength());
            ret.append(cigarElement.getOperator());
        }
        return ret.toString();
    }

    /**
     * Convert from String CIGAR representation to Cigar class representation.  Does not
     * do validation beyond the most basic CIGAR string well-formedness, i.e. each operator is
     * valid, and preceded by a decimal length.
     * @param textCigar CIGAR in String form ala SAM text file.  "*" means empty CIGAR.
     * @throws RuntimeException if textCigar is invalid at the most basic level.
     * @return cigar in Cigar class format
     */
    public Cigar decode(final String textCigar) {
        if (SAMRecord.NO_ALIGNMENT_CIGAR.equals(textCigar)) {
            return new Cigar();
        }
        final Cigar ret = new Cigar();
        final byte[] cigarBytes = StringUtil.stringToBytes(textCigar);
        for (int i = 0; i < cigarBytes.length; ++i) {
            if (!isDigit(cigarBytes[i])) {
                throw new IllegalArgumentException("Malformed CIGAR string: " + textCigar);
            }
            int length = (cigarBytes[i] - ZERO_BYTE);
            for (++i; isDigit(cigarBytes[i]); ++i) {
                length = (length * 10) + cigarBytes[i] - ZERO_BYTE;
            }
            final CigarOperator operator = CigarOperator.characterToEnum(cigarBytes[i]);
            ret.add(new CigarElement(length, operator));
        }
        return ret;
    }
    
    private boolean isDigit(final byte c) {
        return c >= ZERO_BYTE && c <= NINE_BYTE;
    }

    
        
}

/******************************************************************/
/**************************[END OF TextCigarCodec.java]*************************/
/******************************************************************/
