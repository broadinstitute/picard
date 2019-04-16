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

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.util.MathUtils;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import picard.util.MathUtil;

import java.io.File;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static picard.fingerprint.HaplotypeProbabilities.Genotype.HET_ALLELE12;
import static picard.fingerprint.HaplotypeProbabilities.Genotype.HOM_ALLELE1;
import static picard.fingerprint.HaplotypeProbabilities.Genotype.HOM_ALLELE2;

/**
 * Calculates various metrics on a sample fingerprint, indicating whether the fingerprint satisfies the assumptions we have.
 * For example, if too many sites are heterozygous, that would get flagged.
 *
 * @author Yossi Farjoun
 */

@CommandLineProgramProperties(
        summary = CalculateFingerprintMetrics.USAGE_SUMMARY + CalculateFingerprintMetrics.USAGE_DETAILS,
        oneLineSummary = CalculateFingerprintMetrics.USAGE_SUMMARY,
        programGroup = DiagnosticsAndQCProgramGroup.class
)

@DocumentedFeature
public class CalculateFingerprintMetrics extends CommandLineProgram {

    static final String USAGE_DETAILS =
            "This tools collects various statistics that pertain to a single fingerprint (<b>not</b> the comparison, or " +
                    "\'fingerprinting\' of two distinct samples) and reports the results in a metrics file. " +
                    "<p>" +
                    "The statistics collected are p-values, where the null-hypothesis is that the fingerprint is collected from " +
                    "a non-contaminated, diploid human, whose genotypes are modelled by the probabilities given in the " +
                    "HAPLOTYPE_MAP file." +
                    "<p>" +
                    "<h3>Example</h3>\n" +
                    "<pre>\" +\n" +
                    "java -jar picard.jar CalculateFingerprintMetrics \\\n" +
                    "      INPUT=sample.bam \\\n" +
                    "      HAPLOTYPE_DATABASE=fingerprinting_haplotype_database.txt \\\n" +
                    "      OUTPUT=sample.fingerprint_metrics\n" +
                    " </pre>\n";

    static final String USAGE_SUMMARY="Calculate statistics on fingerprints, checking their viability";

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "One or more input files (SAM/BAM/CRAM or VCF).")
    public List<String> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output files to write.")
    public File OUTPUT;

    @Argument(shortName = "H", doc = "The file lists a set of SNPs, optionally arranged in high-LD blocks, to be used for fingerprinting. See " +
            "https://software.broadinstitute.org/gatk/documentation/article?id=9526 for details.")
    public File HAPLOTYPE_MAP;

    @Argument(doc = "Specificies which data-type should be used as the basic unit. Fingerprints from readgroups can " +
            "be \"rolled-up\" to the LIBRARY, SAMPLE, or FILE level before being used." +
            " Fingerprints from VCF can be be examined by SAMPLE or FILE.")
    public CrosscheckMetric.DataType CALCULATE_BY = CrosscheckMetric.DataType.READGROUP;

    @Argument(doc="LOD score threshold for considering a genotype to be definitive.")
    public final double GENOTYPE_LOD_THRESHOLD = 3;

    @Argument(doc="Number of randomization trials for calculating the DISCRIMINATORY_POWER metric.")
    public final int NUMBER_OF_SAMPLING = 100;

    // a fixed random seed for reproducibility;
    private static final int RANDOM_SEED = 42;
    private static final ChiSquareTest chiSquareTest = new ChiSquareTest();

    @Override
    protected int doWork() {

        final List<Path> inputPaths = IOUtil.getPaths(INPUT);
        IOUtil.assertPathsAreReadable(inputPaths);
        IOUtil.assertFileIsReadable(HAPLOTYPE_MAP);
        IOUtil.assertFileIsWritable(OUTPUT);

        final FingerprintChecker checker = new FingerprintChecker(HAPLOTYPE_MAP);

        final MetricsFile<FingerprintMetrics, ?> metricsFile = getMetricsFile();
        final Map<FingerprintIdDetails, Fingerprint> fpMap = checker.fingerprintFiles(inputPaths, 1, 1, TimeUnit.DAYS);
        final Map<FingerprintIdDetails, Fingerprint> mergedFpMap = Fingerprint.mergeFingerprintsBy(fpMap,Fingerprint.getFingerprintIdDetailsStringFunction(CALCULATE_BY));

        metricsFile.addAllMetrics(mergedFpMap.values().stream().map(this::getFingerprintMetrics).collect(Collectors.toList()));
        metricsFile.write(OUTPUT);

        return 0;
    }

    public FingerprintMetrics getFingerprintMetrics(final Fingerprint fingerprint) {

        final RandomGenerator rg = new MersenneTwister(RANDOM_SEED);
        //get expectation of counts and expected fractions
        double[] genotypeCounts = new double[] {0, 0, 0};
        double[] expectedRatios = new double[] {0, 0, 0};

        for (final HaplotypeProbabilities haplotypeProbabilities : fingerprint.values()) {
            genotypeCounts = MathUtil.sum(genotypeCounts, haplotypeProbabilities.getPosteriorProbabilities());
            expectedRatios = MathUtil.sum(expectedRatios, haplotypeProbabilities.getPriorProbablities());
        }

        final double[] hetVsHomExpect = new double[]{expectedRatios[HOM_ALLELE1.v] + expectedRatios[HOM_ALLELE2.v], expectedRatios[HET_ALLELE12.v]};
        final double[] hetVsHomCounts = new double[]{genotypeCounts[HOM_ALLELE1.v] + genotypeCounts[HOM_ALLELE2.v], genotypeCounts[HET_ALLELE12.v]};

        final double[] homAllele1VsAllele2Expect = new double[]{expectedRatios[HOM_ALLELE1.v], expectedRatios[HOM_ALLELE2.v]};
        final double[] homAllele1VsAllele2Counts = new double[]{genotypeCounts[HOM_ALLELE1.v], genotypeCounts[HOM_ALLELE2.v]};

        final long[] roundedHomRefVsVarCounts = MathUtil.round(homAllele1VsAllele2Counts);
        final long[] roundedHetVsHomCounts = MathUtil.round(hetVsHomCounts);
        final long[] roundedGenotypeCounts = MathUtil.round(genotypeCounts);

        final FingerprintMetrics fingerprintMetrics = new FingerprintMetrics();

        fingerprintMetrics.SAMPLE_ALIAS = fingerprint.getSample();
        fingerprintMetrics.SOURCE = Optional.ofNullable(fingerprint.getSource()).map(p->p.toUri().toString()).orElse("");
        fingerprintMetrics.INFO = fingerprint.getInfo();
        fingerprintMetrics.HAPLOTYPES = fingerprint.values().size();
        fingerprintMetrics.HAPLOTYPES_WITH_EVIDENCE = fingerprint.values().stream().filter(HaplotypeProbabilities::hasEvidence).count();
        fingerprintMetrics.DEFINITE_GENOTYPES = fingerprint.values().stream().filter(h -> h.getLodMostProbableGenotype() > GENOTYPE_LOD_THRESHOLD).count();
        fingerprintMetrics.NUM_HOM_ALLELE1 = roundedGenotypeCounts[HOM_ALLELE1.v];
        fingerprintMetrics.NUM_HET = roundedGenotypeCounts[HET_ALLELE12.v];
        fingerprintMetrics.NUM_HOM_ALLELE2 = roundedGenotypeCounts[HOM_ALLELE2.v];

        // calculate p-value
        final double chiSquaredTest = chiSquareTest.chiSquareTest(expectedRatios, roundedGenotypeCounts);
        fingerprintMetrics.CHI_SQUARED_PVALUE = chiSquaredTest;
        fingerprintMetrics.LOG10_CHI_SQUARED_PVALUE = Math.log10(chiSquaredTest);

        // calculate LOD (cross-entropy)
        fingerprintMetrics.CROSS_ENTROPY_LOD = MathUtil.klDivergance(genotypeCounts, expectedRatios);

        // calculate p-value
        final double hetsChiSquaredTest = chiSquareTest.chiSquareTest(hetVsHomExpect, roundedHetVsHomCounts);
        fingerprintMetrics.HET_CHI_SQUARED_PVALUE = hetsChiSquaredTest;
        fingerprintMetrics.HET_LOG10_CHI_SQUARED_PVALUE = Math.log10(hetsChiSquaredTest);

        // calculate LOD (cross-entropy)
        fingerprintMetrics.HET_CROSS_ENTROPY_LOD = MathUtil.klDivergance(hetVsHomCounts, hetVsHomExpect);

        // calculate p-value
        final double homsChiSquaredTest = chiSquareTest.chiSquareTest(homAllele1VsAllele2Expect, roundedHomRefVsVarCounts);
        fingerprintMetrics.HOM_CHI_SQUARED_PVALUE = homsChiSquaredTest;
        fingerprintMetrics.HOM_LOG10_CHI_SQUARED_PVALUE = Math.log10(homsChiSquaredTest);

        // calculate LOD (cross-entropy)
        fingerprintMetrics.HOM_CROSS_ENTROPY_LOD = MathUtil.klDivergance(homAllele1VsAllele2Counts, homAllele1VsAllele2Expect);;

        final MathUtils.RunningStat randomTrials = new MathUtils.RunningStat();

        // get a bunch of random permutations of the fingerprint and compare to self
        IntStream.range(0, NUMBER_OF_SAMPLING).forEach(i ->
                randomTrials.push(FingerprintChecker.calculateMatchResults(fingerprint, randomizeFingerprint(fingerprint, rg)).getLOD()));

        fingerprintMetrics.LOD_SELF_CHECK = FingerprintChecker.calculateMatchResults(fingerprint, fingerprint).getLOD();
        fingerprintMetrics.DISCRIMINATORY_POWER = fingerprintMetrics.LOD_SELF_CHECK - randomTrials.mean();

        return fingerprintMetrics;
    }

    /** Creates a new fingerprint from the current one by randomizing the probabilities within each haplotype
     *
     */
    private static Fingerprint randomizeFingerprint(final Fingerprint fingerprint, final RandomGenerator rg) {
        final Fingerprint retVal = new Fingerprint(null, null, null);

        final RandomDataGenerator rng = new RandomDataGenerator(rg);

        fingerprint.forEach((key, hp) -> {
            final HaplotypeProbabilitiesFromGenotypeLikelihoods permutedHaplotypeProbabilities = new HaplotypeProbabilitiesFromGenotypeLikelihoods(key);
            permutedHaplotypeProbabilities.addToLogLikelihoods(
                    hp.getRepresentativeSnp(),
                    hp.getRepresentativeSnp().getAlleles(),
                    MathUtil.permute(hp.getLogLikelihoods(), rng));
            retVal.add(permutedHaplotypeProbabilities);
        });

        return retVal;
    }
}
