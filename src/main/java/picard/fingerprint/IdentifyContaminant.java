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

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.argumentcollections.ReferenceArgumentCollection;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.util.Map;

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

    @Argument(shortName = "C", doc = "A value of estimated contamination (must be between 0 and 1) in the file. ", minValue = 0D, maxValue = 1D)
    public double CONTAMINATION;

    @Argument(doc = "The sample alias to associate with the resulting fingerprint. When null, uses \"<SAMPLE>-contamination\", where <SAMPLE> is extracted from the input SAM/BAM.", optional = true)
    public String SAMPLE_ALIAS = null;

    @Argument(doc = "The maximum number of reads to use as evidence for any given locus. This is provided as a way to limit the " +
            "effect that any given locus may have.")
    public int LOCUS_MAX_READS = 200;

    @Argument(doc = "Extract a fingerprint for the contaminated sample (instead of the contaminant). Setting to true changes the effect of SAMPLE_ALIAS when null. It names the sample in the VCF <SAMPLE>-contaminated, using the SM value from the SAM header.")
    public boolean EXTRACT_CONTAMINATED = false;

    private final Log log = Log.getInstance(IdentifyContaminant.class);

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(HAPLOTYPE_MAP);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);

        final FingerprintChecker checker = new FingerprintChecker(HAPLOTYPE_MAP);

        //if we want the contaminated fingerprint instead, we need to change the value of CONTAMINATION:
        if (EXTRACT_CONTAMINATED) CONTAMINATION = 1 - CONTAMINATION;

        final Map<String, Fingerprint> fingerprintMap = checker.identifyContaminant(INPUT, CONTAMINATION, LOCUS_MAX_READS);
        if (fingerprintMap.size() != 1) {
            log.error("Expected exactly 1 fingerprint, found " + fingerprintMap.size());
            throw new IllegalArgumentException("Expected exactly 1 fingerprint in Input file, found " + fingerprintMap.size());
        }

        final Map.Entry<String, Fingerprint> sampleFingerprintEntry = fingerprintMap.entrySet().iterator().next();
        final String sampleToUse;
        if (SAMPLE_ALIAS == null) {
            // there must be exactly one element in fingerprintMap, this is checked above.
            final String sample = sampleFingerprintEntry.getKey() + "-contamination";
            if (EXTRACT_CONTAMINATED) {
                sampleToUse = String.format("%s-contaminated", sample);
            } else {
                sampleToUse = String.format("%s-contamination", sample);
            }

        } else {
            sampleToUse = SAMPLE_ALIAS;
        }

        try {
            FingerprintUtils.writeFingerPrint(sampleFingerprintEntry.getValue(), OUTPUT, REFERENCE_SEQUENCE, sampleToUse, "PLs derived from " + INPUT + " using an assumed contamination of " + this.CONTAMINATION);
        } catch (Exception e) {
            log.error(e);
        }

        return 0;
    }

    // Return a custom (required) Reference Argument Collection
    @Override
    protected ReferenceArgumentCollection makeReferenceArgumentCollection() {
        return new ReferenceArgumentCollection() {

            @Argument(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc = "The reference sequence to map the genotypes to.")
            public File REFERENCE_SEQUENCE;

            @Override
            public File getReferenceFile() {
                return REFERENCE_SEQUENCE;
            }
        };
    }

    @Override
    protected boolean requiresReference() {
        return true;
    }
}
