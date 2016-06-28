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
package picard.fingerprint;

import java.util.HashMap;
import java.util.Map;

import picard.PicardException;

import htsjdk.samtools.util.StringUtil;

/**
 * A genotype produced by one of the concrete implementations of AbstractAlleleCaller.
 * DO NOT ADD TO OR REORDER THIS ENUM AS THAT WOULD BREAK THE GELI FILE FORMAT.
 */
public enum DiploidGenotype {
    AA('A','A'),
    AC('A','C'),
    AG('A','G'),
    AT('A','T'),
    CC('C','C'),
    CG('C','G'),
    CT('C','T'),
    GG('G','G'),
    GT('G','T'),
    TT('T','T');
    
    private static final Map<Integer, DiploidGenotype> genotypes = new HashMap<Integer, DiploidGenotype>();
    
    static {
        for (final DiploidGenotype genotype : values()) {
            // this relies on the fact that the integer sum of allele1 and allele2 is unique
            if (genotypes.put(genotype.allele1 + genotype.allele2, genotype) != null) {
                // this check is just for safety, this should never happen
                throw new PicardException("sum of allele values are not unique!!!");
            }
        }
    }
    
    /** Converts a pair of bases into a DiploidGenotype regardless of base order or case */
    public static DiploidGenotype fromBases(final byte[] bases) {
        if (bases.length != 2) {
            throw new IllegalArgumentException("bases must contain 2 and only 2 bases, it actually contained " + bases.length);
        }
        return fromBases(bases[0], bases[1]);
    }

    /** Converts a pair of bases into a DiploidGenotype regardless of base order or case */
    public static DiploidGenotype fromBases(final byte base1, final byte base2) {
        final byte first = StringUtil.toUpperCase(base1);
        final byte second = StringUtil.toUpperCase(base2);
        final DiploidGenotype genotype = genotypes.get(first + second);
        if (genotype == null) {
            throw new IllegalArgumentException("Unknown genotype string [" + 
                    StringUtil.bytesToString(new byte[] {base1, base2}) +
                    "], any pair of ACTG case insensitive is acceptable");
        }
        return genotype;
    }

    /**
     * @return true if this is a valid base, i.e. one of [ACGTacgt]
     */
    public static boolean isValidBase(final byte base) {
        switch(StringUtil.toUpperCase(base)) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
                return true;
            default:
                return false;
        }
    }

    private final byte allele1;
    private final byte allele2;

    private DiploidGenotype(final char allele1, final char allele2) {
        this.allele1 = (byte)(allele1 & 0xff);
        this.allele2 = (byte)(allele2 & 0xff);
    }
    
    public byte getAllele1() { return allele1; }
    public byte getAllele2() { return allele2; }
    public boolean isHeterozygous() { return this.allele1 != this.allele2; }
    public boolean isHomomozygous() { return this.allele1 == this.allele2; }
}
