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

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import picard.PicardException;
import picard.analysis.FingerprintingDetailMetrics;
import picard.analysis.FingerprintingSummaryMetrics;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Alpha;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * Attempts to check the sample identity of the sequence data in the provided SAM/BAM file
 * against a set of known genotypes in the supplied genotype file (in either GELI or VCF format).
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage = CheckFingerprint.USAGE_DETAILS,
        usageShort = "Computes a fingerprint from the supplied SAM/BAM file and compares it to the provided genotypes",
        programGroup = Alpha.class // TODO -- when mature please move to a to-be-created Fingerprinting.class
)
public class CheckFingerprint extends CommandLineProgram {

    static final String USAGE_DETAILS = "Computes a fingerprint from the supplied SAM/BAM file and " +
            "compares it to the expected fingerprint genotypes provided. The key output is a LOD score " +
            "which represents the relative likelihood of the sequence data originating from the same " +
            "sample as the genotypes vs. from a random sample.  Two outputs are produced: (1) a summary " +
            "metrics file that gives metrics related single read group (lane or index within a lane) " +
            "versus a set of known genotypes for the expected sample, and (2) a detail metrics file that " +
            "contains an individual SNP/Haplotype comparison within a fingerprint comparison.  The two " +
            "files may be specified individually using the SUMMARY_OUTPUT and DETAIL_OUTPUT options.  " +
            "Alternatively the OUTPUT option may be used instead to give the base of the two output " +
            "files, with the summary metrics having a file extension '" + CheckFingerprint.FINGERPRINT_SUMMARY_FILE_SUFFIX + "' " +
            "and the detail metrics having a file extension '" + CheckFingerprint.FINGERPRINT_DETAIL_FILE_SUFFIX + "'.";

    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM or BAM file.")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The base of output files to write.  The summary metrics " +
            "will have the file extension '" + CheckFingerprint.FINGERPRINT_SUMMARY_FILE_SUFFIX + "' and the detail metrics will have " +
            "the extension '" + CheckFingerprint.FINGERPRINT_DETAIL_FILE_SUFFIX + "'.",  mutex = {"SUMMARY_OUTPUT", "DETAIL_OUTPUT"})
    public String OUTPUT;

    @Option(shortName = "S", doc = "The text file to which to write summary metrics.", mutex = {"OUTPUT"})
    public File SUMMARY_OUTPUT;

    @Option(shortName = "D", doc = "The text file to which to write detail metrics.",  mutex = {"OUTPUT"})
    public File DETAIL_OUTPUT;

    @Option(shortName="G", doc = "File of genotypes (VCF or GELI) to be used in comparison. May contain " +
            "any number of genotypes; CheckFingerprint will use only those that are usable for fingerprinting.")
    public File GENOTYPES;

    @Option(optional=true, doc = "If using VCF format genotypes, this parameter can be used to specify which sample's " +
            "genotypes to use from the VCF file.  If not supplied the sample name from the BAM read group header " +
            "is used instead.")
    public String SAMPLE_ALIAS;

    @Option(shortName="H", doc = "A file of haplotype information produced by the CheckFingerprint program.")
    public File HAPLOTYPE_MAP;

    @Option(shortName="LOD", doc = "When counting haplotypes checked and matching, count only haplotypes " +
            "where the most likely haplotype achieves at least this LOD.")
    public double GENOTYPE_LOD_THRESHOLD = 5;

    @Option(shortName="IGNORE_RG", doc = "If true, treat the entire input BAM as one single read group in the calculation, " +
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

        SequenceUtil.assertSequenceDictionariesEqual(SAMSequenceDictionaryExtractor.extractDictionary(INPUT), SAMSequenceDictionaryExtractor.extractDictionary(GENOTYPES), true);
        SequenceUtil.assertSequenceDictionariesEqual(SAMSequenceDictionaryExtractor.extractDictionary(INPUT), checker.getHeader().getSequenceDictionary(), true);

        // If sample alias isn't supplied, assume it's the one from the INPUT file's RGs
        if (SAMPLE_ALIAS == null) {
            final SamReader in = SamReaderFactory.makeDefault().open(INPUT);
            for (final SAMReadGroupRecord rec : in.getFileHeader().getReadGroups()) {
                if (SAMPLE_ALIAS == null) {
                    SAMPLE_ALIAS = rec.getSample();
                }
                else if (!SAMPLE_ALIAS.equals(rec.getSample())) {
                    throw new PicardException("SAM File must not contain data from multiple samples.");
                }
            }
            CloserUtil.close(in);
        }


        final List<FingerprintResults> results = checker.checkFingerprints(
                Arrays.asList(INPUT),
                Arrays.asList(GENOTYPES),
                SAMPLE_ALIAS,
                IGNORE_READ_GROUPS);

        final MetricsFile<FingerprintingSummaryMetrics,?> summaryFile = getMetricsFile();
        final MetricsFile<FingerprintingDetailMetrics,?> detailsFile = getMetricsFile();

        for (final FingerprintResults fpr : results) {
            final MatchResults mr = fpr.getMatchResults().first();

            final FingerprintingSummaryMetrics metrics = new FingerprintingSummaryMetrics();
            metrics.READ_GROUP                = fpr.getReadGroup();
            metrics.SAMPLE                    = SAMPLE_ALIAS;
            metrics.LL_EXPECTED_SAMPLE = mr.getSampleLikelihood();
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
                details.READ_GROUP = fpr.getReadGroup();
                details.SAMPLE     = SAMPLE_ALIAS;
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
            log.info(metrics.READ_GROUP + " vs. " + metrics.SAMPLE + ": LOD = " + metrics.LOD_EXPECTED_SAMPLE);
        }

        summaryFile.write(outputSummaryMetricsFile);
        detailsFile.write(outputDetailMetricsFile);

        return 0;
    }
}
