/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
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

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;

/**
 * Program to create a fingerprint for the <b>contaminating</b> sample when the level of contamination is both known and
 * uniform in the genome.
 *
 * @author Yossi Farjoun
 */
@CommandLineProgramProperties(
        summary = "Computes the fingerprint genotype likelihoods from the supplied SAM/BAM file and a contamination estimate." +
                "NOTA BENE: the fingerprint is provided for the contamination (by default) for the main sample. " +
                "It is given as a list of PLs at the fingerprinting sites.",
        oneLineSummary = "Computes a fingerprint from the supplied SAM/BAM file, given a contamination estimate.",
        programGroup = DiagnosticsAndQCProgramGroup.class)

public class IdentifyContaminant extends CommandLineProgram {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM or BAM file.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output fingerprint file (VCF).")
    public File OUTPUT;

    @Argument(shortName = "H", doc = "A file of haplotype information. The file lists a set of SNPs, optionally arranged in high-LD blocks, to be used for fingerprinting. See " +
            "https://software.broadinstitute.org/gatk/documentation/article?id=9526 for details.")
    public File HAPLOTYPE_MAP;

    @Argument(shortName = "C", doc = "A value of estimated contamination in the input. ", minValue = 0D, maxValue = 1D)
    public double CONTAMINATION;

    @Argument(doc = "The sample alias to associate with the resulting fingerprint. When null, <SAMPLE> is extracted from the input file and \"<SAMPLE>-contamination\" is used.", optional = true)
    public String SAMPLE_ALIAS = null;

    @Argument(doc = "The maximum number of reads to use as evidence for any given locus. This is provided as a way to limit the " +
            "effect that any given locus may have.")
    public int LOCUS_MAX_READS = 200;

    @Argument(doc = "Extract a fingerprint for the contaminated sample (instead of the contaminant). Setting to true changes the effect of SAMPLE_ALIAS when null. " +
            "It names the sample in the VCF <SAMPLE>, using the SM value from the SAM header.")
    public boolean EXTRACT_CONTAMINATED = false;

    @Override
    protected boolean requiresReference() {
        return true;
    }

    @Override
    protected int doWork() {

        final ExtractFingerprint extractFingerprint = new ExtractFingerprint();

        extractFingerprint.INPUT = INPUT;
        extractFingerprint.OUTPUT = OUTPUT;
        extractFingerprint.EXTRACT_CONTAMINATION = !EXTRACT_CONTAMINATED;
        extractFingerprint.HAPLOTYPE_MAP = HAPLOTYPE_MAP;
        extractFingerprint.CONTAMINATION = CONTAMINATION;
        extractFingerprint.LOCUS_MAX_READS = LOCUS_MAX_READS;
        extractFingerprint.SAMPLE_ALIAS = SAMPLE_ALIAS;
        extractFingerprint.VALIDATION_STRINGENCY = VALIDATION_STRINGENCY;
        extractFingerprint.VERBOSITY = VERBOSITY;
        extractFingerprint.referenceSequence = referenceSequence;

        extractFingerprint.doWork();
        return 0;
    }
}
