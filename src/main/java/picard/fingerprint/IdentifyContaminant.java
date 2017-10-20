/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.argumentcollections.ReferenceArgumentCollection;
import picard.cmdline.programgroups.Fingerprinting;

import java.io.File;
import java.io.IOException;
import java.util.Map;

/**
 * Program to create a fingerprint for the <b>contaminating</b> sample when the level of contamination is both known and
 * uniform in the genome.
 *
 * @author Yossi Farjoun
 */
@CommandLineProgramProperties(
        summary = "Computes the fingerprint genotype likelihoods from the supplied SAM/BAM file and a contamination estimate." +
                "NOTA BENE: the fingerprint is provided for the contamination (by default) or the main sample. " +
                "It is given as a list of PLs at the fingerprinting sites.",
        oneLineSummary = "Computes a fingerprint (of the contaminating or contaminated sample) from the supplied SAM/BAM file.",
        programGroup = Fingerprinting.class)

public class IdentifyContaminant extends CommandLineProgram {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM or BAM file.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Name of output fingerprint file (VCF)")
    public File OUTPUT;

    @Argument(shortName = "H", doc = "A file of haplotype information.")
    public File HAPLOTYPE_MAP;

    @Argument(shortName = "C", doc = "A value of estimated contamination (must be between 0 and 1) in the file. ", optional = false)
    public double CONTAMINATION;

    @Argument(doc = "The sample alias to associate with the resulting fingerprint. Default value is <SAMPLE>-contamination, where <SAMPLE> is extracted from the header.", optional = true)
    public String SAMPLE_ALIAS;

    @Argument(doc = "The maximum number of reads to use as evidence for any given locus", optional = true)
    public int LOCUS_MAX_READS = 200;

    @Argument(doc = "Extract the contaminated fingerprint instead of the contaminant. If true changes the default value of SAMPLE_ALIAS to be the actual sample-id found in the header.", optional = true)
    public boolean EXTRACT_CONTAMINATED = false;

    private final Log log = Log.getInstance(IdentifyContaminant.class);

    // Stock main method
    public static void main(final String[] args) {
        new IdentifyContaminant().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(HAPLOTYPE_MAP);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);

        if (CONTAMINATION < 0.0 || CONTAMINATION > 1.0) {
            final PicardException e = new PicardException("Value for contamination must be in the range 0-1, value received=" + CONTAMINATION);

            log.error(e);
            throw e;
        }

        final SAMFileHeader header = getHeader();

        checkSingleSample(header);

        final FingerprintChecker checker = new FingerprintChecker(HAPLOTYPE_MAP);

        //if we want the contaminated fingerprint instead, we need to change the value of CONTAMINATION:
        if (EXTRACT_CONTAMINATED) CONTAMINATION = 1 - CONTAMINATION;

        final Map<String, Fingerprint> fingerprintMap = checker.identifyContaminant(INPUT, CONTAMINATION, LOCUS_MAX_READS);
        if (fingerprintMap.size() != 1) {
            log.error("Expected exactly 1 fingerprint, found " + fingerprintMap.size());
            throw new IllegalArgumentException("Expected exactly 1 fingerprint in Input file, found " + fingerprintMap.size());
        }

        final Fingerprint results = fingerprintMap.values().iterator().next();

        try {
            FingerprintUtils.writeFingerPrint(results, OUTPUT, REFERENCE_SEQUENCE, SAMPLE_ALIAS, "PLs derived from " + INPUT + " using an assumed contamination of " + this.CONTAMINATION);
        } catch(Exception e){
            log.error(e);
        }

        return 0;
    }

    // Return a custom (required)  argument collection
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

    // take the the sample in the first RG and check that the others are the same.
    // create a <SAMPLE>-contamination sample id from it.
    private String getSampleAlias(final SAMFileHeader header) {

        if (header.getReadGroups() == null || header.getReadGroups().isEmpty() || header.getReadGroups().get(0).getSample() == null) {
            log.error("Cannot find a read-group in bam header. This program requires that input files that include a read-group with an @SM tag.");
            throw new IllegalArgumentException("Input file doesn't have a readgroup with an @SM tag: Not supported.");
        }
        final String bamSample = header.getReadGroups().get(0).getSample();
        for (int i = 1; i < header.getReadGroups().size(); i++) {
            if (!bamSample.equals(header.getReadGroups().get(i).getSample())) {
                log.error("File has more than one sample in it (or a null sample): '%s' and '%s'. Not supported.", bamSample, header.getReadGroups().get(i).getSample());
                throw new IllegalArgumentException("File has more than one sample in it: Not supported.");
            }
        }

        final String sampleToUse;
        if (EXTRACT_CONTAMINATED) {
            sampleToUse = bamSample;
        } else {
            sampleToUse = String.format("%s-contamination", bamSample);
        }
        log.info("No sample name provided, Using '%s'", sampleToUse);

        return sampleToUse;
    }

    private SAMFileHeader getHeader() {
        final SamReader in = SamReaderFactory.makeDefault().open(INPUT);
        final SAMFileHeader header = in.getFileHeader();
        try {
            in.close();
        } catch (final IOException e) {
            log.error(e);
        }
        return header;
    }

    private void checkSingleSample(final SAMFileHeader header) {
        String sampleInFile = null;
        for (final SAMReadGroupRecord rec : header.getReadGroups()) {
            if (sampleInFile == null) {
                sampleInFile = rec.getSample();
            } else if (!sampleInFile.equals(rec.getSample())) {
                final PicardException e = new PicardException("SAM File must not contain data from multiple samples. Found: " + sampleInFile + " and " + rec.getSample());
                log.error(e);
                throw e;
            }
        }
    }
}
