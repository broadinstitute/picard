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

import picard.PicardException;

import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * Represents information about a group of SNPs that form a haplotype in perfect LD
 * with one another.
 *
 * @author Tim Fennell
 */
public class HaplotypeBlock implements Comparable<HaplotypeBlock> {
    private final double maf;
    private final Map<String,Snp> snpsByName     = new HashMap<String,Snp>();
    private final double[] haplotypeFrequencies  = new double[3];

    private Snp firstSnp;
    private String chrom;
    private int start;
    private int end;

    /** Constructs a haplotype block with the provided minor allele frequency. */
    public HaplotypeBlock(final double maf) {
        this.maf = maf;

        // Set the haplotype frequencies assuming hardy-weinberg
        final double majorAf = (1 - maf);
        this.haplotypeFrequencies[0] = majorAf * majorAf;
        this.haplotypeFrequencies[1] = majorAf * maf * 2;
        this.haplotypeFrequencies[2] = maf * maf;
    }

    /** Gets the set of haplotype frequencies. */
    public double[] getHaplotypeFrequencies() { return this.haplotypeFrequencies; }

    /** Adds a SNP to the haplotype.  Will throw an exception if the SNP is on the wrong chromosome. */
    public void addSnp(final Snp snp) {
        if (this.snpsByName.isEmpty()) {
            this.chrom = snp.getChrom();
            this.start = snp.getPos();
            this.end   = snp.getPos();
            this.firstSnp = snp;
        }
        else if (!this.chrom.equals(snp.getChrom())) {
            throw new PicardException("Snp chromosome " + snp.getChrom() +
                " does not agree with chromosome of existing snp(s): " + this.chrom);
        }
        else {
            if (snp.getPos() < this.start) {
                this.start = snp.getPos();
                this.firstSnp = snp;
            }
            if (snp.getPos() > this.end) {
                this.end = snp.getPos();
            }
        }

        this.snpsByName.put(snp.getName(), snp);
    }

    /** Gets a SNP by name if it belongs to this haplotype. */
    public Snp getSnp(final String name) {
        return this.snpsByName.get(name);
    }

    /** Gets the arbitrarily first SNP in the haplotype. */
    public Snp getFirstSnp() {
        return this.firstSnp;
    }

    /** Returns true if the SNP is contained within the haplotype block, false otherwise. */
    public boolean contains(final Snp snp) {
        // Check is performed on SNP name and position because of the fact that some SNPs
        // have multiple mappings in the genome and we're paranoid!
        final Snp contained = this.snpsByName.get(snp.getName());
        return contained != null && contained.getChrom().equals(snp.getChrom()) &&
                contained.getPos() == snp.getPos();
    }

    /** Returns the number of SNPs within the haplotype block. */
    public int size() {
        return snpsByName.size();
    }    

    /** Returns an unmodifiable, unordered, collection of all SNPs in this haplotype block. */
    public Collection<Snp> getSnps() {
        return Collections.unmodifiableCollection(this.snpsByName.values());
    }

    /**
     * Gets the frequency of the i'th diploid haplotype where haplotypes are ordered accorinding
     * to DiploidHaplotype.
     */
    public double getHaplotypeFrequency(final int i) {
        if (i < 0 || i > 2) throw new IllegalArgumentException("Illegal haplotype index " + i);
        else return this.haplotypeFrequencies[i];
    }

    /** Returns the minor allele frequency of this haplotype. */
    public double getMaf() { return this.maf; }

    /**
     * Gets the expected genotype of the provided SNP given the provided haplotype of this
     * haplotype block.
     */
    public DiploidGenotype getSnpGenotype(final Snp snp, final DiploidHaplotype haplotype) {
        if (!contains(snp)) throw new IllegalArgumentException("Snp is not part of haplotype " + snp);
        return snp.getGenotype(haplotype);
    }

    /**
     * Gets the diploid haplotype for this haplotype block given the provided SNP and SNP
     * genotype.
     */
    public DiploidHaplotype getDiploidHaplotype(final Snp snp, final DiploidGenotype gt) {
        if (!contains(snp)) throw new IllegalArgumentException("Snp is not part of haplotype " + snp);
        return DiploidHaplotype.values()[snp.indexOf(gt)];
    }

    @Override
    public int compareTo(final HaplotypeBlock that) {
        int retval = this.chrom.compareTo(that.chrom);
        if (retval == 0) retval = this.start - that.start;
        if (retval == 0) retval = this.end   - that.end;
        return retval;
    }

    @Override public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        else return this.compareTo((HaplotypeBlock) o) == 0;
    }

    @Override public int hashCode() {
        return this.start;
    }

    @Override public String toString() {
        return this.chrom + "[" + this.start + "-" + this.end + "]";
    }
}
