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


/**
 * Class DiploidGenotype
 *
 * Enum describing all possible combinations of diploid genotype variations;
 * AA, AC, etc.
 *
 * @author aaron
 */
public enum DiploidGenotype {
    AA, AC, AG, AT, CC, CG, CT, GG, GT, TT;

    public static DiploidGenotype toDiploidGenotype(String genotype) {
        if (genotype.length() != 2)
            throw new DiploidGenotypeException("Genotype string for conversion should be of length 2, we were passed = " + genotype);
        genotype = genotype.toUpperCase();
        for (DiploidGenotype g: DiploidGenotype.values())
            if (g.toString().equals(genotype)) return g;
        throw new DiploidGenotypeException("Unable to find genotype matching " + genotype);
    }

    public boolean isHet() {
        return toString().toCharArray()[0] != toString().toCharArray()[1];
    }

    public boolean containsBase(char base) {
        return (toString().charAt(0) == base || toString().charAt(1) == base);
    }
}

class DiploidGenotypeException extends RuntimeException {
    DiploidGenotypeException(String s) {
        super(s);
    }

    DiploidGenotypeException(String s, Throwable throwable) {
        super(s, throwable);
    }
}