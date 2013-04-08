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

import org.broad.tribble.Feature;

import java.util.Arrays;


/**
 *         <p/>
 *         Class GeliTextFeature
 *         <p/>
 *         This is a feature for the Geli text object, which is the text version of the Geli binary genotyping format.
 *
 * @author aaron
 */
public class GeliTextFeature implements Feature {

    private final String contig;                // the contig name
    private final long position;                // the position on the contig
    private final char refBase;                 // the reference base
    private final int depthOfCoverage;          // the depth of coverage at this position
    private final int maximumMappingQual;       // the maximum mapping quality of a read at this position
    private final DiploidGenotype genotype;     // the called genotype
    private final double LODBestToReference;    // the LOD score of the best to the reference
    private final double LODBestToNext;         // the LOD score of the best to the next best genotype
    private final double likelihoods[];         // the array of all genotype likelihoods, in ordinal order

    /**
     * Create a geli text feature, given:
     *
     * @param contig             the contig
     * @param position           the position on the contig
     * @param refBase            the reference base
     * @param depthOfCoverage    the depth of coverage at this position
     * @param maximumMappingQual the maximum mapping quality of a read at this position
     * @param genotype           the called genotype
     * @param LODBestToReference the LOD score of the best to the reference
     * @param LODBestToNext      the LOD score of the best to the next best genotype
     * @param likelihoods        the array of all genotype likelihoods, in ordinal ordering
     */
    public GeliTextFeature(String contig,
                           long position,
                           char refBase,
                           int depthOfCoverage,
                           int maximumMappingQual,
                           DiploidGenotype genotype,
                           double LODBestToReference,
                           double LODBestToNext,
                           double[] likelihoods) {
        this.contig = contig;
        this.position = position;
        this.refBase = refBase;
        this.depthOfCoverage = depthOfCoverage;
        this.maximumMappingQual = maximumMappingQual;
        this.genotype = genotype;
        this.LODBestToReference = LODBestToReference;
        this.LODBestToNext = LODBestToNext;
        this.likelihoods = likelihoods;
    }

    /** Return the features reference sequence name, e.g chromosome or contig */
    public String getChr() {
        return this.contig;
    }

    /** Return the start position in 1-based coordinates (first base is 1) */
    public int getStart() {
        return (int) this.position;
    }

    /**
     * Return the end position following 1-based fully closed conventions.  The length of a feature is
     * end - start + 1;
     */
    public int getEnd() {
        return (int) this.position;
    }

    public char getRefBase() {
        return refBase;
    }

    public int getDepthOfCoverage() {
        return depthOfCoverage;
    }

    public int getMaximumMappingQual() {
        return maximumMappingQual;
    }

    public DiploidGenotype getGenotype() {
        return genotype;
    }

    public double getLODBestToNext() {
        return LODBestToNext;
    }

    public double getLODBestToReference() {
        return LODBestToReference;
    }

    public double[] getLikelihoods() {
        return likelihoods;
    }

    private static double Epsilon = 0.0001;
    public boolean equals(Object o) {
        if (!(o instanceof GeliTextFeature)) return false;
        GeliTextFeature other = (GeliTextFeature)o;
        if (!Arrays.equals(likelihoods,other.likelihoods)) return false;
        if (!contig.equals(other.contig)) return false;
        if (!(position == other.position)) return false;
        if (!(refBase == other.refBase)) return false;
        if (!(depthOfCoverage == other.depthOfCoverage)) return false;
        if (!(maximumMappingQual == other.maximumMappingQual)) return false;
        if (!(genotype == other.genotype)) return false;
        if (!(Math.abs(LODBestToReference - other.LODBestToReference) < Epsilon)) return false;
        if (!(Math.abs(LODBestToNext - other.LODBestToNext) < Epsilon)) return false;
        return true;
    }
}
