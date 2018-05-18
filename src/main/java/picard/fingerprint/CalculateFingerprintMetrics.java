/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.analysis.FingerprintingDetailMetrics;
import picard.analysis.FingerprintingSummaryMetrics;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.text.Normalizer;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Checks the sample identity
 *
 * @author Yossi Farjoun
 */

@CommandLineProgramProperties(
        summary = CalculateFingerprintMetrics.USAGE_DETAILS,
        oneLineSummary = "blah",
        programGroup = DiagnosticsAndQCProgramGroup.class
)

@DocumentedFeature
public class CalculateFingerprintMetrics extends CommandLineProgram {

    static final String USAGE_DETAILS =
            "blah";

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input file SAM/BAM or VCF.  If a VCF is used, " +
            "it must have at least one sample.  If there are more than one samples in the VCF, the parameter OBSERVED_SAMPLE_ALIAS must " +
            "be provided in order to indicate which sample's data to use.  If there are no samples in the VCF, an exception will be thrown.")
    public String INPUT;

      @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output files to write.")
    public File OUTPUT;

    @Argument(shortName = "H", doc = "The file lists a set of SNPs, optionally arranged in high-LD blocks, to be used for fingerprinting. See " +
            "https://software.broadinstitute.org/gatk/documentation/article?id=9526 for details.")
    public File HAPLOTYPE_MAP;

    @Argument(optional = true, shortName = "IGNORE_RG", doc = "If the input is a SAM/BAM, and this parameter is true, treat the " +
            "entire input BAM as one single read group in the calculation, " +
            "ignoring RG annotations, and producing a single fingerprint metric for the entire BAM.")
    public boolean IGNORE_READ_GROUPS = false;

    private final Log log = Log.getInstance(CalculateFingerprintMetrics.class);

    private Path inputPath;

    @Override
    protected int doWork() {

        try {
            inputPath = IOUtil.getPath(INPUT);
        } catch (IOException e) {
            throw new IllegalArgumentException(e);
        }
        IOUtil.assertFileIsReadable(inputPath);
        IOUtil.assertFileIsReadable(HAPLOTYPE_MAP);
        IOUtil.assertFileIsWritable(OUTPUT);

        final FingerprintChecker checker = new FingerprintChecker(HAPLOTYPE_MAP);

        String observedSampleAlias = null;
        final boolean isBamOrSamFile = isBamOrSam(inputPath);
        final MetricsFile<Fingerprint.FingerprintMetrics, ?> metricsFile = getMetricsFile();
        if (isBamOrSamFile) {
            log.info("Comparing Dictionaries");
            SequenceUtil.assertSequenceDictionariesEqual(SAMSequenceDictionaryExtractor.extractDictionary(inputPath), checker.getHeader().getSequenceDictionary(), true);

            // Verify that there's only one sample in the SAM/BAM.
            final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(inputPath);
            for (final SAMReadGroupRecord rec : in.getFileHeader().getReadGroups()) {
                if (observedSampleAlias == null) {
                    observedSampleAlias = rec.getSample();
                } else if (!observedSampleAlias.equals(rec.getSample())) {
                    throw new PicardException("inputPath SAM/BAM file must not contain data from multiple samples.");
                }
            }
            CloserUtil.close(in);
            log.info("Reading fingerprint");

            final Map<FingerprintIdDetails, Fingerprint> fingerprintIdDetailsFingerprintMap = checker.fingerprintSamFile(inputPath, new HaplotypeMap(HAPLOTYPE_MAP).getIntervalList());
            fingerprintIdDetailsFingerprintMap.forEach((fd,f)-> metricsFile.addMetric(f.getFingerprintMetrics()));

            final Fingerprint combinedFp = new Fingerprint(observedSampleAlias, inputPath, "<MERGED>");
            fingerprintIdDetailsFingerprintMap.values().forEach(combinedFp::merge);

            metricsFile.addMetric(combinedFp.getFingerprintMetrics());
        } else {            // Input is a VCF
            // Note that FingerprintChecker.loadFingerprints() verifies that the VCF's Sequence Dictionaries agree with that of the Haplotye Map File

            // Verify that there is only one sample in the VCF
            final VCFFileReader fileReader = new VCFFileReader(inputPath, false);
            final VCFHeader fileHeader = fileReader.getFileHeader();
            if (fileHeader.getNGenotypeSamples() < 1) {
                throw new PicardException("inputPath VCF file must contain at least one sample.");
            }
            log.info("Loading fingerprints");
            final Map<String, Fingerprint> fingerprintMap = checker.loadFingerprints(inputPath, null);

            fingerprintMap.forEach((fd, f) -> metricsFile.addMetric(f.getFingerprintMetrics()));
        }
        metricsFile.write(OUTPUT);

        return 0;
    }

    @Override
    protected String[] customCommandLineValidation() {
        try {
            if (IGNORE_READ_GROUPS && !isBamOrSam(IOUtil.getPath(INPUT))) {
                return new String[]{"The parameter IGNORE_READ_GROUPS can only be used with BAM/SAM inputs."};
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        return super.customCommandLineValidation();
    }

    static boolean isBamOrSam(final Path p) {
        return (p.toUri().getRawPath().endsWith(BamFileIoUtils.BAM_FILE_EXTENSION) || p.toUri().getRawPath().endsWith(IOUtil.SAM_FILE_EXTENSION));
    }
}
