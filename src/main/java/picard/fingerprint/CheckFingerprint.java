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
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.analysis.FingerprintingDetailMetrics;
import picard.analysis.FingerprintingSummaryMetrics;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Fingerprinting;

import java.io.File;
import java.util.Collections;
import java.util.List;

/**
 * Attempts to check the sample identity of the sequence/genotype data in the provided file (SAM/BAM or VCF)
 * against a set of known genotypes in the supplied genotype file (in either GELI or VCF format).
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = CheckFingerprint.USAGE_DETAILS,
        oneLineSummary = "Computes a fingerprint from the supplied input (SAM/BAM or VCF) file and compares it to the provided genotypes",
        programGroup = Fingerprinting.class
)
public class CheckFingerprint extends CommandLineProgram {

    static final String USAGE_DETAILS = "Computes a fingerprint from the supplied input file (SAM/BAM or VCF) file and " +
            "compares it to the expected fingerprint genotypes provided. The key output is a LOD score " +
            "which represents the relative likelihood of the sequence data originating from the same " +
            "sample as the genotypes vs. from a random sample.  Two outputs are produced: (1) a summary " +
            "metrics file that gives metrics at the single sample level (if the input was a VCF) or at the read " +
            "level (lane or index within a lane) (if the input was a SAM/BAM) " +
            "versus a set of known genotypes for the expected sample, and (2) a detail metrics file that " +
            "contains an individual SNP/Haplotype comparison within a fingerprint comparison.  The two " +
            "files may be specified individually using the SUMMARY_OUTPUT and DETAIL_OUTPUT options.  " +
            "Alternatively the OUTPUT option may be used instead to give the base of the two output " +
            "files, with the summary metrics having a file extension '" + CheckFingerprint.FINGERPRINT_SUMMARY_FILE_SUFFIX + "' " +
            "and the detail metrics having a file extension '" + CheckFingerprint.FINGERPRINT_DETAIL_FILE_SUFFIX + "'.";

    @Argument(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input file SAM/BAM or VCF.  If a VCF is used, " +
            "it must have at least one sample.  If there are more than one samples in the VCF, the parameter OBSERVED_SAMPLE_ALIAS must " +
            "be provided in order to indicate which sample's data to use.  If there are no samples in the VCF, an exception will be thrown.")
    public File INPUT;

    @Argument(optional = true, doc = "If the input is a VCF, this parameters used to select which sample's data in the VCF to use.")
    public String OBSERVED_SAMPLE_ALIAS;

    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The base prefix of output files to write.  The summary metrics " +
            "will have the file extension '" + CheckFingerprint.FINGERPRINT_SUMMARY_FILE_SUFFIX + "' and the detail metrics will have " +
            "the extension '" + CheckFingerprint.FINGERPRINT_DETAIL_FILE_SUFFIX + "'.",  mutex = {"SUMMARY_OUTPUT", "DETAIL_OUTPUT"})
    public String OUTPUT;

    @Argument(shortName = "S", doc = "The text file to which to write summary metrics.", mutex = {"OUTPUT"})
    public File SUMMARY_OUTPUT;

    @Argument(shortName = "D", doc = "The text file to which to write detail metrics.",  mutex = {"OUTPUT"})
    public File DETAIL_OUTPUT;

    @Argument(shortName="G", doc = "File of genotypes (VCF or GELI) to be used in comparison. May contain " +
            "any number of genotypes; CheckFingerprint will use only those that are usable for fingerprinting.")
    public File GENOTYPES;

    @Argument(shortName = "SAMPLE_ALIAS", optional=true, doc = "This parameter can be used to specify which sample's genotypes to use from the " +
            "expected VCF file (the GENOTYPES file).  If it is not supplied, the sample name from the input " +
            "(VCF or BAM read group header) will be used.")
    public String EXPECTED_SAMPLE_ALIAS;

    @Argument(shortName="H", doc = "The file lists a set of SNPs, optionally arranged in high-LD blocks, to be used for fingerprinting. See " +
            "https://software.broadinstitute.org/gatk/documentation/article?id=9526 for details.")
    public File HAPLOTYPE_MAP;

    @Argument(shortName="LOD", doc = "When counting haplotypes checked and matching, count only haplotypes " +
            "where the most likely haplotype achieves at least this LOD.")
    public double GENOTYPE_LOD_THRESHOLD = 5;

    @Argument(optional=true, shortName="IGNORE_RG", doc = "If the input is a SAM/BAM, and this parameter is true, treat the " +
            "entire input BAM as one single read group in the calculation, " +
            "ignoring RG annotations, and producing a single fingerprint metric for the entire BAM.")
    public boolean IGNORE_READ_GROUPS = false;

    private final Log log = Log.getInstance(CheckFingerprint.class);

    public static final String FINGERPRINT_SUMMARY_FILE_SUFFIX = "fingerprinting_summary_metrics";
    public static final String FINGERPRINT_DETAIL_FILE_SUFFIX = "fingerprinting_detail_metrics";


    // Stock main method
    public static void main(final String[] args) {
        new CheckFingerprint().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        final File outputDetailMetricsFile, outputSummaryMetricsFile;
        if (OUTPUT == null) {
            outputDetailMetricsFile = DETAIL_OUTPUT;
            outputSummaryMetricsFile = SUMMARY_OUTPUT;
        }
        else {
            if (!OUTPUT.endsWith(".")) OUTPUT = OUTPUT + ".";
            outputDetailMetricsFile = new File(OUTPUT + FINGERPRINT_DETAIL_FILE_SUFFIX);
            outputSummaryMetricsFile = new File(OUTPUT + FINGERPRINT_SUMMARY_FILE_SUFFIX);
        }

        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(HAPLOTYPE_MAP);
        IOUtil.assertFileIsReadable(GENOTYPES);
        IOUtil.assertFileIsWritable(outputDetailMetricsFile);
        IOUtil.assertFileIsWritable(outputSummaryMetricsFile);

        final FingerprintChecker checker = new FingerprintChecker(HAPLOTYPE_MAP);
        List<FingerprintResults> results;

        String observedSampleAlias = null;
        final boolean isBamOrSamFile = isBamOrSamFile(INPUT);
        if (isBamOrSamFile) {
            SequenceUtil.assertSequenceDictionariesEqual(SAMSequenceDictionaryExtractor.extractDictionary(INPUT), SAMSequenceDictionaryExtractor.extractDictionary(GENOTYPES), true);
            SequenceUtil.assertSequenceDictionariesEqual(SAMSequenceDictionaryExtractor.extractDictionary(INPUT), checker.getHeader().getSequenceDictionary(), true);

            // Verify that there's only one sample in the SAM/BAM.
            final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
            for (final SAMReadGroupRecord rec : in.getFileHeader().getReadGroups()) {
                if (observedSampleAlias == null) {
                    observedSampleAlias = rec.getSample();
                }
                else if (!observedSampleAlias.equals(rec.getSample())) {
                    throw new PicardException("INPUT SAM/BAM file must not contain data from multiple samples.");
                }
            }
            CloserUtil.close(in);

            // If expected sample alias isn't supplied, assume it's the one from the INPUT file's RGs
            if (EXPECTED_SAMPLE_ALIAS == null) {
                EXPECTED_SAMPLE_ALIAS = observedSampleAlias;
            }

            results = checker.checkFingerprints(
                    Collections.singletonList(INPUT),
                    Collections.singletonList(GENOTYPES),
                    EXPECTED_SAMPLE_ALIAS,
                    IGNORE_READ_GROUPS);
        } else {            // Input is a VCF
            // Note that FingerprintChecker.loadFingerprints() verifies that the VCF's Sequence Dictionaries agree with that of the Haplotye Map File

            // Verify that there is only one sample in the VCF
            final VCFFileReader fileReader = new VCFFileReader(INPUT, false);
            final VCFHeader fileHeader = fileReader.getFileHeader();
            if (fileHeader.getNGenotypeSamples() < 1) {
                throw new PicardException("INPUT VCF file must contain at least one sample.");
            }
            if ((fileHeader.getNGenotypeSamples() > 1) && (OBSERVED_SAMPLE_ALIAS == null)) {
                throw new PicardException("INPUT VCF file contains multiple samples and yet the OBSERVED_SAMPLE_ALIAS parameter is not set.");
            }
            // set observedSampleAlias to the parameter, if set.  Otherwise, if here, this must be a single sample VCF, get it's sample
            observedSampleAlias = (OBSERVED_SAMPLE_ALIAS != null) ? OBSERVED_SAMPLE_ALIAS : fileHeader.getGenotypeSamples().get(0);

            // Now verify that observedSampleAlias is, in fact, in the VCF
            if (!fileHeader.getGenotypeSamples().contains(observedSampleAlias)) {
                throw new PicardException("INPUT VCF file does not contain OBSERVED_SAMPLE_ALIAS: " + observedSampleAlias);
            }

            if (OBSERVED_SAMPLE_ALIAS == null) {
                observedSampleAlias = fileHeader.getGenotypeSamples().get(0);
            }
            fileReader.close();

            // If expected sample alias isn't supplied, assume it's the one from the INPUT file
            if (EXPECTED_SAMPLE_ALIAS == null) {
                EXPECTED_SAMPLE_ALIAS = observedSampleAlias;
            }

            results = checker.checkFingerprints(
                    Collections.singletonList(INPUT),
                    Collections.singletonList(GENOTYPES),
                    observedSampleAlias,
                    EXPECTED_SAMPLE_ALIAS);
        }

        final MetricsFile<FingerprintingSummaryMetrics,?> summaryFile = getMetricsFile();
        final MetricsFile<FingerprintingDetailMetrics,?> detailsFile = getMetricsFile();

        for (final FingerprintResults fpr : results) {
            final MatchResults mr = fpr.getMatchResults().first();

            final FingerprintingSummaryMetrics metrics = new FingerprintingSummaryMetrics();
            metrics.READ_GROUP                = fpr.getReadGroup();
            metrics.SAMPLE                    = EXPECTED_SAMPLE_ALIAS;
            metrics.LL_EXPECTED_SAMPLE        = mr.getSampleLikelihood();
            metrics.LL_RANDOM_SAMPLE          = mr.getPopulationLikelihood();
            metrics.LOD_EXPECTED_SAMPLE = mr.getLOD();

            for (final LocusResult lr : mr.getLocusResults()) {
                final DiploidGenotype expectedGenotype = lr.getExpectedGenotype();
                final DiploidGenotype observedGenotype = lr.getMostLikelyGenotype();
                // Update the summary metrics
                metrics.HAPLOTYPES_WITH_GENOTYPES++;
                if (lr.getLodGenotype() >= GENOTYPE_LOD_THRESHOLD) {
                    metrics.HAPLOTYPES_CONFIDENTLY_CHECKED++;

                    if (lr.getExpectedGenotype() == lr.getMostLikelyGenotype()) {
                        metrics.HAPLOTYPES_CONFIDENTLY_MATCHING++;
                    }

                    if (expectedGenotype.isHeterozygous() && observedGenotype.isHomomozygous()) {
                        metrics.HET_AS_HOM++;
                    }

                    if (expectedGenotype.isHomomozygous() && observedGenotype.isHeterozygous()) {
                        metrics.HOM_AS_HET++;
                    }

                    if (expectedGenotype.isHomomozygous() && observedGenotype.isHomomozygous()
                            && expectedGenotype.compareTo(observedGenotype) != 0) {
                        metrics.HOM_AS_OTHER_HOM++;
                    }
                }

                // Build the detail metrics
                final FingerprintingDetailMetrics details = new FingerprintingDetailMetrics();
                details.READ_GROUP         = fpr.getReadGroup();
                details.SAMPLE             = EXPECTED_SAMPLE_ALIAS;
                details.SNP                = lr.getSnp().getName();
                details.SNP_ALLELES        = lr.getSnp().getAlleleString();
                details.CHROM              = lr.getSnp().getChrom();
                details.POSITION           = lr.getSnp().getPos();
                details.EXPECTED_GENOTYPE  = expectedGenotype.toString();
                details.OBSERVED_GENOTYPE  = observedGenotype.toString();
                details.LOD                = lr.getLodGenotype();
                details.OBS_A              = lr.getAllele1Count();
                details.OBS_B              = lr.getAllele2Count();
                detailsFile.addMetric(details);
            }

            summaryFile.addMetric(metrics);
            log.info("Read Group: " + metrics.READ_GROUP + " / " + observedSampleAlias + " vs. " + metrics.SAMPLE + ": LOD = " + metrics.LOD_EXPECTED_SAMPLE);
        }

        summaryFile.write(outputSummaryMetricsFile);
        detailsFile.write(outputDetailMetricsFile);

        return 0;
    }

    protected String[] customCommandLineValidation() {
        IOUtil.assertFileIsReadable(INPUT);

        boolean isBamOrSamFile = isBamOrSamFile(INPUT);
        if (!isBamOrSamFile && IGNORE_READ_GROUPS) {
            return new String[]{"The parameter IGNORE_READ_GROUPS can only be used with BAM/SAM inputs."};
        }
        if (isBamOrSamFile && OBSERVED_SAMPLE_ALIAS != null) {
            return new String[]{"The parameter OBSERVED_SAMPLE_ALIAS can only be used with a VCF input."};
        }
        return super.customCommandLineValidation();
    }

    static boolean isBamOrSamFile(final File f) {
        return (BamFileIoUtils.isBamFile(f) || f.getName().endsWith(IOUtil.SAM_FILE_EXTENSION));
    }
}
