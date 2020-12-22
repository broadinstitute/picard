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
import java.util.Collections;
import java.util.List;

/**
 * Checks the sample identity of the sequence/genotype data in the provided file (SAM/BAM or VCF)
 * against a set of known genotypes in the supplied genotype file (in VCF format).
 * <p>
 * <h3> Summary </h3>
 * Computes a fingerprint (essentially, genotype information from different parts of the genome)
 * from the supplied input file (SAM/BAM or VCF) file and
 * compares it to the expected fingerprint genotypes provided. The key output is a LOD score
 * which represents the relative likelihood of the sequence data originating from the same
 * sample as the genotypes vs. from a random sample.
 * <br/>
 * Two outputs are produced:
 * <ol>
 * <li>A summary metrics file that gives metrics of the fingerprint matches when comparing the input to a set of
 * genotypes for the expected sample.  At the single sample level (if the input was a VCF) or at the read
 * level (lane or index within a lane) (if the input was a SAM/BAM) </li>
 * <li>A detail metrics file that contains an individual SNP/Haplotype comparison within a fingerprint comparison.
 * </li>
 * </ol>
 * The metrics files fill the fields of the classes {@link FingerprintingSummaryMetrics} and
 * {@link FingerprintingDetailMetrics}.
 * The output files may be specified individually using the SUMMARY_OUTPUT and DETAIL_OUTPUT options.
 * Alternatively the OUTPUT option may be used instead to give the base of the two output
 * files, with the summary metrics having a file extension {@value #FINGERPRINT_SUMMARY_FILE_SUFFIX},
 * and the detail metrics having a file extension {@value #FINGERPRINT_DETAIL_FILE_SUFFIX}.
 * <br/>
 * <h3>Example comparing a bam against known genotypes:</h3>
 * <pre>
 *     java -jar picard.jar CheckFingerprint \
 *          INPUT=sample.bam \
 *          GENOTYPES=sample_genotypes.vcf \
 *          HAPLOTYPE_MAP=fingerprinting_haplotype_database.txt \
 *          OUTPUT=sample_fingerprinting
 * </pre>
 * <br/>
 * <h3>Detailed Explanation</h3>
 * <p>
 * This tool calculates a single number that reports the LOD score for identity check between the {@link #INPUT}
 * and the {@link #GENOTYPES}. A positive value indicates that the data seems to have come from the same individual
 * or, in other words the identity checks out. The scale is logarithmic (base 10), so a LOD of 6 indicates
 * that it is 1,000,000 more likely that the data matches the genotypes than not. A negative value indicates
 * that the data do not match. A score that is near zero is inconclusive and can result from low coverage
 * or non-informative genotypes.
 * <p>
 * The identity check makes use of haplotype blocks defined in the {@link #HAPLOTYPE_MAP} file to enable it to have higher
 * statistical power for detecting identity or swap by aggregating data from several SNPs in the haplotype block. This
 * enables an identity check of samples with very low coverage (e.g. ~1x mean coverage).
 * <p>
 * When provided a VCF, the identity check looks at the PL, GL and GT fields (in that order) and uses the first one that
 * it finds.
 *
 * @author Tim Fennell
 */

@CommandLineProgramProperties(
        summary = CheckFingerprint.USAGE_DETAILS,
        oneLineSummary = "Computes a fingerprint from the supplied input (SAM/BAM/CRAM or VCF) file and compares it to the provided genotypes",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@DocumentedFeature
public class CheckFingerprint extends CommandLineProgram {

    static final String USAGE_DETAILS =
            "Checks the sample identity of the sequence/genotype data in the provided file (SAM/BAM/CRAM or VCF) " +
                    "against a set of known genotypes in the supplied genotype file (in VCF format).\n " +
                    "\n " +
                    "<h3>Summary</h3> " +
                    "Computes a fingerprint (essentially, genotype information from different parts of the genome) " +
                    "from the supplied input file (SAM/BAM/CRAM or VCF) file and " +
                    "compares it to the expected fingerprint genotypes provided. " +
                    "The key output is a LOD score which represents the relative likelihood of " +
                    "the sequence data originating from the same sample as the genotypes vs. from a random sample. " +
                    "<br/> " +
                    "Two outputs are produced: " +
                    "<ol> " +
                    "<li>" +
                    "A summary metrics file that gives metrics of the fingerprint matches when comparing the input to a set of " +
                    "genotypes for the expected sample.  " +
                    "At the single sample level (if the input was a VCF) or at the read level (lane or index within a lane) " +
                    "(if the input was a SAM/BAM) " +
                    "</li> " +
                    "<li>" +
                    "A detail metrics file that contains an individual SNP/Haplotype comparison within a fingerprint comparison." +
                    "</li> " +
                    "</ol> " +
                    "The metrics files fill the fields of the classes FingerprintingSummaryMetrics and FingerprintingDetailMetrics. " +
                    "The output files may be specified individually using the SUMMARY_OUTPUT and DETAIL_OUTPUT options. " +
                    "Alternatively the OUTPUT option may be used instead to give the base of the two output files, " +
                    "with the summary metrics having a file extension \".fingerprinting_summary_metrics\", " +
                    "and the detail metrics having a file extension \".fingerprinting_detail_metrics\". " +
                    "<br/> " +
                    "<h3>Example comparing a bam against known genotypes:</h3> " +
                    "<pre> " +
                    "    java -jar picard.jar CheckFingerprint \\\n " +
                    "         INPUT=sample.bam \\\n " +
                    "         GENOTYPES=sample_genotypes.vcf \\\n " +
                    "         HAPLOTYPE_MAP=fingerprinting_haplotype_database.txt \\\n " +
                    "         OUTPUT=sample_fingerprinting " +
                    "</pre> " +
                    "<br/> " +
                    "<h3>Detailed Explanation</h3>" +
                    "This tool calculates a single number that reports the LOD score for identity check between the INPUT and the GENOTYPES. " +
                    "A positive value indicates that the data seems to have come from the same individual or, " +
                    "in other words the identity checks out. " +
                    "The scale is logarithmic (base 10), " +
                    "so a LOD of 6 indicates that it is 1,000,000 more likely that the data matches the genotypes than not. " +
                    "A negative value indicates that the data do not match. " +
                    "A score that is near zero is inconclusive and can result from low coverage or non-informative genotypes. " +
                    "\n\n " +
                    "The identity check makes use of haplotype blocks defined in the HAPLOTYPE_MAP file to enable it to have higher " +
                    "statistical power for detecting identity or swap by aggregating data from several SNPs in the haplotype block. " +
                    "This enables an identity check of samples with very low coverage (e.g. ~1x mean coverage). " +
                    "\n\n " +
                    "When provided a VCF, the identity check looks at the PL, GL and GT fields (in that order) " +
                    "and uses the first one that it finds. ";

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input file SAM/BAM/CRAM or VCF.  If a VCF is used, " +
            "it must have at least one sample.  If there are more than one samples in the VCF, the parameter OBSERVED_SAMPLE_ALIAS must " +
            "be provided in order to indicate which sample's data to use.  If there are no samples in the VCF, an exception will be thrown.")
    public String INPUT;

    @Argument(optional = true, doc = "If the input is a VCF, this parameters used to select which sample's data in the VCF to use.")
    public String OBSERVED_SAMPLE_ALIAS;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The base prefix of output files to write.  The summary metrics " +
            "will have the file extension '" + CheckFingerprint.FINGERPRINT_SUMMARY_FILE_SUFFIX + "' and the detail metrics will have " +
            "the extension '" + CheckFingerprint.FINGERPRINT_DETAIL_FILE_SUFFIX + "'.", mutex = {"SUMMARY_OUTPUT", "DETAIL_OUTPUT"})
    public String OUTPUT;

    @Argument(shortName = "S", doc = "The text file to which to write summary metrics.", mutex = {"OUTPUT"})
    public File SUMMARY_OUTPUT;

    @Argument(shortName = "D", doc = "The text file to which to write detail metrics.", mutex = {"OUTPUT"})
    public File DETAIL_OUTPUT;

    @Argument(shortName = "G", doc = "File of genotypes (VCF) to be used in comparison. May contain " +
            "any number of genotypes; CheckFingerprint will use only those that are usable for fingerprinting.")
    public String GENOTYPES;

    @Argument(shortName = "SAMPLE_ALIAS", optional = true, doc = "This parameter can be used to specify which sample's genotypes to use from the " +
            "expected VCF file (the GENOTYPES file).  If it is not supplied, the sample name from the input " +
            "(VCF or BAM read group header) will be used.")
    public String EXPECTED_SAMPLE_ALIAS;

    @Argument(shortName = "H", doc = "The file lists a set of SNPs, optionally arranged in high-LD blocks, to be used for fingerprinting. See " +
            "https://software.broadinstitute.org/gatk/documentation/article?id=9526 for details.")
    public File HAPLOTYPE_MAP;

    @Argument(shortName = "LOD", doc = "When counting haplotypes checked and matching, count only haplotypes " +
            "where the most likely haplotype achieves at least this LOD.")
    public double GENOTYPE_LOD_THRESHOLD = 5;

    @Argument(optional = true, shortName = "IGNORE_RG", doc = "If the input is a SAM/BAM/CRAM, and this parameter is true, treat the " +
            "entire input BAM as one single read group in the calculation, " +
            "ignoring RG annotations, and producing a single fingerprint metric for the entire BAM.")
    public boolean IGNORE_READ_GROUPS = false;

    private final Log log = Log.getInstance(CheckFingerprint.class);

    public static final String FINGERPRINT_SUMMARY_FILE_SUFFIX = "fingerprinting_summary_metrics";
    public static final String FINGERPRINT_DETAIL_FILE_SUFFIX = "fingerprinting_detail_metrics";

    private Path inputPath;
    private Path genotypesPath;

    @Override
    protected int doWork() {
        final File outputDetailMetricsFile, outputSummaryMetricsFile;
        if (OUTPUT == null) {
            outputDetailMetricsFile = DETAIL_OUTPUT;
            outputSummaryMetricsFile = SUMMARY_OUTPUT;
        } else {
            OUTPUT += ".";
            outputDetailMetricsFile = new File(OUTPUT + FINGERPRINT_DETAIL_FILE_SUFFIX);
            outputSummaryMetricsFile = new File(OUTPUT + FINGERPRINT_SUMMARY_FILE_SUFFIX);
        }

        try {
            inputPath = IOUtil.getPath(INPUT);
            genotypesPath = IOUtil.getPath(GENOTYPES);
        } catch (IOException e) {
            throw new IllegalArgumentException(e);
        }
        IOUtil.assertFileIsReadable(inputPath);
        IOUtil.assertFileIsReadable(HAPLOTYPE_MAP);
        IOUtil.assertFileIsReadable(genotypesPath);
        IOUtil.assertFileIsWritable(outputDetailMetricsFile);
        IOUtil.assertFileIsWritable(outputSummaryMetricsFile);

        final FingerprintChecker checker = new FingerprintChecker(HAPLOTYPE_MAP);
        checker.setReferenceFasta(REFERENCE_SEQUENCE);
        List<FingerprintResults> results;

        String observedSampleAlias = null;
        if (fileContainsReads(inputPath)) {
            SequenceUtil.assertSequenceDictionariesEqual(SAMSequenceDictionaryExtractor.extractDictionary(inputPath), SAMSequenceDictionaryExtractor.extractDictionary(genotypesPath), true);
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

            // If expected sample alias isn't supplied, assume it's the one from the INPUT file's RGs
            if (EXPECTED_SAMPLE_ALIAS == null) {
                EXPECTED_SAMPLE_ALIAS = observedSampleAlias;
            }

            results = checker.checkFingerprints(
                    Collections.singletonList(inputPath),
                    Collections.singletonList(genotypesPath),
                    EXPECTED_SAMPLE_ALIAS,
                    IGNORE_READ_GROUPS);
        } else {            // Input is a VCF
            // Note that FingerprintChecker.loadFingerprints() verifies that the VCF's Sequence Dictionaries agree with that of the Haplotye Map File

            // Verify that there is only one sample in the VCF
            final VCFFileReader fileReader = new VCFFileReader(inputPath, false);
            final VCFHeader fileHeader = fileReader.getFileHeader();
            if (fileHeader.getNGenotypeSamples() < 1) {
                throw new PicardException("inputPath VCF file must contain at least one sample.");
            }
            if ((fileHeader.getNGenotypeSamples() > 1) && (OBSERVED_SAMPLE_ALIAS == null)) {
                throw new PicardException("inputPath VCF file contains multiple samples and yet the OBSERVED_SAMPLE_ALIAS parameter is not set.");
            }
            // set observedSampleAlias to the parameter, if set.  Otherwise, if here, this must be a single sample VCF, get it's sample
            observedSampleAlias = (OBSERVED_SAMPLE_ALIAS != null) ? OBSERVED_SAMPLE_ALIAS : fileHeader.getGenotypeSamples().get(0);

            // Now verify that observedSampleAlias is, in fact, in the VCF
            if (!fileHeader.getGenotypeSamples().contains(observedSampleAlias)) {
                throw new PicardException("inputPath VCF file does not contain OBSERVED_SAMPLE_ALIAS: " + observedSampleAlias);
            }

            if (OBSERVED_SAMPLE_ALIAS == null) {
                observedSampleAlias = fileHeader.getGenotypeSamples().get(0);
            }
            fileReader.close();

            // If expected sample alias isn't supplied, assume it's the one from the INPUT file
            if (EXPECTED_SAMPLE_ALIAS == null) {
                EXPECTED_SAMPLE_ALIAS = observedSampleAlias;
            }

            results = checker.checkFingerprintsFromPaths(
                    Collections.singletonList(inputPath),
                    Collections.singletonList(genotypesPath),
                    observedSampleAlias,
                    EXPECTED_SAMPLE_ALIAS);
        }

        final MetricsFile<FingerprintingSummaryMetrics, ?> summaryFile = getMetricsFile();
        final MetricsFile<FingerprintingDetailMetrics, ?> detailsFile = getMetricsFile();

        boolean allZeroLod = true;
        for (final FingerprintResults fpr : results) {
            final MatchResults mr = fpr.getMatchResults().first();

            final FingerprintingSummaryMetrics metrics = new FingerprintingSummaryMetrics();
            metrics.READ_GROUP = fpr.getReadGroup();
            metrics.SAMPLE = EXPECTED_SAMPLE_ALIAS;
            metrics.LL_EXPECTED_SAMPLE = mr.getSampleLikelihood();
            metrics.LL_RANDOM_SAMPLE = mr.getPopulationLikelihood();
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
                details.READ_GROUP = fpr.getReadGroup();
                details.SAMPLE = EXPECTED_SAMPLE_ALIAS;
                details.SNP = lr.getSnp().getName();
                details.SNP_ALLELES = lr.getSnp().getAlleleString();
                details.CHROM = lr.getSnp().getChrom();
                details.POSITION = lr.getSnp().getPos();
                details.EXPECTED_GENOTYPE = expectedGenotype.toString();
                details.OBSERVED_GENOTYPE = observedGenotype.toString();
                details.LOD = lr.getLodGenotype();
                details.OBS_A = lr.getAllele1Count();
                details.OBS_B = lr.getAllele2Count();
                detailsFile.addMetric(details);
            }

            summaryFile.addMetric(metrics);
            log.info("Read Group: " + metrics.READ_GROUP + " / " + observedSampleAlias + " vs. " + metrics.SAMPLE + ": LOD = " + metrics.LOD_EXPECTED_SAMPLE);

            allZeroLod &= metrics.LOD_EXPECTED_SAMPLE == 0;
        }

        summaryFile.write(outputSummaryMetricsFile);
        detailsFile.write(outputDetailMetricsFile);

        if (allZeroLod) {
            log.error("No non-zero results found. This is likely an error. " +
                    "Probable cause: EXPECTED_SAMPLE (if provided) or the sample name from INPUT (if EXPECTED_SAMPLE isn't provided)" +
                    "isn't a sample in GENOTYPES file.");
            return 1;
        }

        return 0;
    }

    protected String[] customCommandLineValidation() {

        try {
            final boolean fileContainsReads = fileContainsReads(IOUtil.getPath(INPUT));
            if (!fileContainsReads && IGNORE_READ_GROUPS) {
                return new String[]{"The parameter IGNORE_READ_GROUPS can only be used with BAM/SAM/CRAM inputs."};
            }
            if (fileContainsReads && OBSERVED_SAMPLE_ALIAS != null) {
                return new String[]{"The parameter OBSERVED_SAMPLE_ALIAS can only be used with a VCF input."};
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        if (REFERENCE_SEQUENCE == null && INPUT.endsWith(SamReader.Type.CRAM_TYPE.fileExtension())) {
            return new String[]{"REFERENCE must be provided when using CRAM as input."};
        }

        return super.customCommandLineValidation();
    }

    static boolean fileContainsReads(final Path p) {
        return (p.toUri().getRawPath().endsWith(SamReader.Type.BAM_TYPE.fileExtension()) ||
                p.toUri().getRawPath().endsWith(SamReader.Type.SAM_TYPE.fileExtension()) ||
                p.toUri().getRawPath().endsWith(SamReader.Type.CRAM_TYPE.fileExtension()));
    }
}
