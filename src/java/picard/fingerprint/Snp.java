/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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

import htsjdk.samtools.util.StringUtil;

import java.util.ArrayList;
import java.util.List;

/**
 * Class to represent a SNP in context of a haplotype block that is used in fingerprinting.
 *
 * @author Tim Fennell
 */
public class Snp implements Comparable<Snp> {
    private final String name;
    private final String chrom;
    private final int pos;
    private final byte allele1;
    private final byte allele2;
    private final double maf; // technically the allele frequency of allele2
    private final List<String> fingerprintPanels;

    private final DiploidGenotype[] genotypes = new DiploidGenotype[3];

    public Snp(final String name, final String chrom, final int pos, final byte allele1, final byte allele2,
               final double maf, final List<String> fingerprintPanels) {
        this.name = name;
        this.chrom = chrom;
        this.pos = pos;
        this.allele1 = StringUtil.toUpperCase(allele1);
        this.allele2 = StringUtil.toUpperCase(allele2);
        this.maf = maf;
        this.fingerprintPanels = fingerprintPanels == null ? new ArrayList<String>() : fingerprintPanels;

        // Construct the genotypes for ease of comparison
        this.genotypes[0] = DiploidGenotype.fromBases(allele1, allele1);
        this.genotypes[1] = DiploidGenotype.fromBases(allele1, allele2);
        this.genotypes[2] = DiploidGenotype.fromBases(allele2, allele2);
    }

    /** Returns a new SNP object with the alleles swapped and MAF corrected. */
    public Snp flip() {
        return new Snp(name, chrom, pos, allele2, allele1, 1-maf, fingerprintPanels);
    }

    public String getName() { return name; }
    public String getChrom() { return chrom; }
    public int getPos() { return pos; }
    public byte getAllele1() { return allele1; }
    public byte getAllele2() { return allele2; }
    public double getMaf() { return maf; }
    public List<String> getFingerprintPanels() { return this.fingerprintPanels; }

    public DiploidGenotype getHomozygousAllele1Genotype() { return this.genotypes[0]; }
    public DiploidGenotype getHeterogyzousGenotype() { return this.genotypes[1]; }
    public DiploidGenotype getHomozygousAllele2Genotype() { return this.genotypes[2]; }

    /** Gets the genotype with the given index. */
    DiploidGenotype getGenotype(final DiploidHaplotype haplotype) { return this.genotypes[haplotype.ordinal()]; }

    /** Gets the index of the supplied genotype within the genotypes for this SNP. */
    int indexOf(final DiploidGenotype gt) {
        for (int i=0; i<this.genotypes.length; ++i) {
            if (gt == this.genotypes[i]) return i;
        }

        throw new IllegalArgumentException("Genotype " + gt + " is not valid for this SNP.");
    }

    public String getAlleleString() {
        return StringUtil.bytesToString(new byte[] {allele1, StringUtil.toLowerCase(allele2)});
    }    

    @Override
    public int compareTo(final Snp that) {
        int retval = this.chrom.compareTo(that.chrom);
        if (retval == 0) retval = this.pos - that.pos;
        return retval;
    }

    @Override
    public boolean equals(final Object o) {
        return (this == o) || ((o instanceof Snp) && compareTo((Snp) o) == 0);
    }

    @Override
    public int hashCode() {
        int result = chrom.hashCode();
        result = 31 * result + pos;
        return result;
    }

    @Override
    public String toString() {
        return this.chrom + ":" + this.pos;
    }
}
