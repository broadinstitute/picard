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

import org.broad.tribble.Feature;
import org.broad.tribble.annotation.Strand;


/**
 * @author aaron
 *
 * This class represents a DBSNP track, as retrieved from the UCSC browser.
 *
 * Example format:
 * bin	chrom	chromStart	chromEnd	rsID	score	strand	refNCBI	refUCSC	observed	molType	class	valid	avHet	avHetSE	func	locType	weight
 * 585 chr1 433 433 rs56289060  0  +  - - -/C  genomic  insertion unknown 0  0  unknown  between  1
 * 585 chr1 491 492 rs55998931  0  +  C C C/T  genomic  single   unknown 0 0 unknown exact 1
 */
public class OldDbSNPFeature implements Feature {

    private String contig;                      // our contig location
    private int start;                          // our starting location, zero based
    private int stop;                           // our stopping location

    private String rsID = "";                   // the snp indentifier
    private int score = 0;
    private Strand strand = Strand.NONE;        // Which DNA strand contains the observed alleles

    private String ncbiRefBase = "N";             // the reference base according to NCBI
    private String ucscRefBase = "N";             // the reference base according to UCSC

    private String[] mObserved = null;          // The sequences of the observed alleles

    private String molType = "genomic";         // molecule type
    private String classType = "unknown";       // The class of variant (simple, insertion, deletion, range, etc.)
    private String validationStatus;            // The validation status of the SNP

    private double avHet = 0.0;                 // The average heterozygosity from all observations
    private double avHetSE = 0.0;               // The Standard Error for the average heterozygosity

    private String function = "unknown";        // The functional category of the SNP (coding-synon, coding-nonsynon, intron, etc.)
    private String locationType = "unknown";    // How the variant affects the reference sequence
    private int weight = 0;                     // The quality of the alignment

    /**
     * create the dbSNP feature, given the following information:
     *
     * @param contig the contig rsID
     * @param start  the start position, one based
     * @param stop   the stop position, one based
     */
    OldDbSNPFeature(String contig,
		    int start,
		    int stop) {
        this.contig = contig;
        this.start = start;
        this.stop = stop;
    }

    /*
     * the required getting and setter methods
     */

    public String getChr() {
        return contig;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return stop;
    }

    /*
     * getter and setter methods for the rest of the dbSNP data
     */

    public int getScore() {
        return score;
    }

    public void setScore(int score) {
        this.score = score;
    }

    public Strand getStrand() {
        return strand;
    }

    public void setStrand(Strand strand) {
        this.strand = strand;
    }

    public String getNCBIRefBase() {
        return ncbiRefBase;
    }

    public void setNCBIRefBase(String ncbiRefBase) {
        this.ncbiRefBase = ncbiRefBase;
    }

    public String getUCSCRefBase() {
        return ucscRefBase;
    }

    public void setUCSCRefBase(String ucscRefBase) {
        this.ucscRefBase = ucscRefBase;
    }

    public String[] getObserved() {
        return mObserved;
    }

    public void setObserved(String[] mObserved) {
        this.mObserved = mObserved;
    }

    public String getMolType() {
        return molType;
    }

    public void setMolType(String molType) {
        this.molType = molType;
    }

    public String getVariantType() {
        return classType;
    }

    public void setVariantType(String classType) {
        this.classType = classType;
    }

    public String getValidationStatus() {
        return validationStatus;
    }

    public void setValidationStatus(String validationStatus) {
        this.validationStatus = validationStatus;
    }

    public double getAvHet() {
        return avHet;
    }

    public void setAvHet(double avHet) {
        this.avHet = avHet;
    }

    public double getAvHetSE() {
        return avHetSE;
    }

    public void setAvHetSE(double avHetSE) {
        this.avHetSE = avHetSE;
    }

    public String getFunction() {
        return function;
    }

    public void setFunction(String function) {
        this.function = function;
    }

    public String getLocationType() {
        return locationType;
    }

    public void setLocationType(String locationType) {
        this.locationType = locationType;
    }

    public int getWeight() {
        return weight;
    }

    public void setWeight(int weight) {
        this.weight = weight;
    }

    public String getRsID() {
        return rsID;
    }

    public void setRsID(String rsID) {
        this.rsID = rsID;
    }
}
