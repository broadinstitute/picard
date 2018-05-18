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

import htsjdk.samtools.metrics.MetricBase;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import picard.util.MathUtil;

import java.nio.file.Path;
import java.util.*;

/**
 * Small class to represent a genetic fingerprint as a set of HaplotypeProbabilities
 * objects that give the relative probabilities of each of the possible haplotypes
 * at a locus.
 *
 * @author Tim Fennell
 */
public class Fingerprint extends TreeMap<HaplotypeBlock, HaplotypeProbabilities> {
    private static final double GENOTYPE_LOD_THRESHOLD = 3;
    private static final int NUMBER_OF_SAMPLING = 100;
    private static final int RANDOM_SEED = 42;
    private final String sample;
    private final Path source;
    private final String info;

    public Fingerprint(final String sample, final Path source, final String info) {
        this.sample = sample;
        this.source = source;
        this.info = info;
    }

    public String getSample() { return sample; }

    public Path getSource() { return source; }

    public String getInfo() { return info; }

    public String getPrintableId() {
        return getSample() + "@" + (source == null ? "" : source.toUri().toString()) + (info == null ? "" : (":" + info));
    }

    public void add(final HaplotypeProbabilities h) {
        put(h.getHaplotype(), h);
    }

    /**
     * Merges the likelihoods from the supplied Fingerprint into the likelihoods for this fingerprint.
     */
    public void merge(final Fingerprint other) {
        final Set<HaplotypeBlock> haps = new HashSet<>();
        haps.addAll(keySet());
        haps.addAll(other.keySet());

        for (final HaplotypeBlock haplotype : haps) {
            HaplotypeProbabilities probabilities = get(haplotype);
            final HaplotypeProbabilities otherProbabilities = other.get(haplotype);
            if (probabilities == null) {
                probabilities = otherProbabilities;
                put(haplotype, probabilities);
            } else if (otherProbabilities != null) {
                probabilities.merge(otherProbabilities);
            }
        }
    }

    public FingerprintMetrics getFingerprintMetrics() {

        //get expectation of counts and expected fractions
        double[] genotypeCounts = new double[] {0, 0, 0};
        double[] expectedRatios = new double[] {0, 0, 0};


        for (HaplotypeProbabilities haplotypeProbabilities : this.values()) {
            genotypeCounts = MathUtil.sum(genotypeCounts, haplotypeProbabilities.getPosteriorProbabilities());
            expectedRatios = MathUtil.sum(expectedRatios, haplotypeProbabilities.getPriorProbablities());
        }
        final long[] actualGenotypeCounts = new long[3];
        for (int i = 0; i < HaplotypeProbabilities.Genotype.values().length; i++) {
            // TODO: it would be more accurate to not round the genotype counts, but commons doesn't have a
            // TODO: chi-squared statistic method with double[], double[] signature....(can be remedied)
            actualGenotypeCounts[i] = Math.round(genotypeCounts[i]);
        }
        double[] hetVsHomexpectedRatios = new double[]{expectedRatios[0] + expectedRatios[2], expectedRatios[1]};
        double[] homRefVsVarexpectedRatios = new double[]{expectedRatios[0], expectedRatios[2]};

        double[] hetVsHomCounts = new double[]{actualGenotypeCounts[0] + actualGenotypeCounts[2], actualGenotypeCounts[1]};
        double[] homRefVsVarCounts = new double[]{genotypeCounts[0], genotypeCounts[2]};

        long[] actualHomRefVsVarCounts = new long[]{Math.round(homRefVsVarCounts[0]), Math.round(homRefVsVarCounts[1])};
        long[] actualHetVsHomCounts = new long[]{Math.round(hetVsHomCounts[0]), Math.round(hetVsHomCounts[1])};

        // calculate p-value
        ChiSquareTest chiSquareTest = new ChiSquareTest();
        final double chiSquaredTest = chiSquareTest.chiSquareTest(expectedRatios, actualGenotypeCounts);
        // calculate LOD (cross-entropy)
        final double crossEntropy = -MathUtil.klDivergance(genotypeCounts, expectedRatios);

        // calculate p-value
        final double hetsChiSquaredTest = chiSquareTest.chiSquareTest(hetVsHomexpectedRatios, actualHetVsHomCounts);
        // calculate LOD (cross-entropy)
        final double hetsCrossEntropy = -MathUtil.klDivergance(hetVsHomexpectedRatios, hetVsHomCounts);

        // calculate p-value
        final double homsChiSquaredTest = chiSquareTest.chiSquareTest(homRefVsVarexpectedRatios, actualHomRefVsVarCounts);
        // calculate LOD (cross-entropy)
        final double homsCrossEntropy = -MathUtil.klDivergance(homRefVsVarexpectedRatios, homRefVsVarCounts);



        final double lodSelfCheck = FingerprintChecker.calculateMatchResults(this, this).getLOD();

        final double[] randomizationTrials = new double[NUMBER_OF_SAMPLING];
        RandomGenerator rg = new MersenneTwister(RANDOM_SEED);
        for (int i = 0; i < NUMBER_OF_SAMPLING; i++) {
            randomizationTrials[i] = FingerprintChecker.calculateMatchResults(this, randomizeFingerprint(rg)).getLOD();
        }

        FingerprintMetrics fingerprintMetrics = new FingerprintMetrics();

        fingerprintMetrics.SAMPLE_NAME = sample;
        fingerprintMetrics.SOURCE = source.toUri().toString();
        fingerprintMetrics.INFO = info;
        fingerprintMetrics.HAPLOTYPE = values().size();
        fingerprintMetrics.HAPLOTYPES_WITH_EVIDENCE = values().stream().filter(HaplotypeProbabilities::hasEvidence).count();
        fingerprintMetrics.DEFINITE_GENOTYPES = values().stream().filter(h -> h.getLodMostProbableGenotype() > GENOTYPE_LOD_THRESHOLD).count();
        fingerprintMetrics.NUM_HOM_REF = actualGenotypeCounts[0];
        fingerprintMetrics.NUM_HET = actualGenotypeCounts[1];
        fingerprintMetrics.NUM_HOM_VAR = actualGenotypeCounts[2];
        fingerprintMetrics.CHI_SQUARED_PVALUE = chiSquaredTest;
        fingerprintMetrics.LOG10_CHI_SQUARED_PVALUE = Math.log10(chiSquaredTest);
        fingerprintMetrics.CROSS_ENTROPY_LOD = crossEntropy;
        fingerprintMetrics.HET_CHI_SQUARED_PVALUE = hetsChiSquaredTest;
        fingerprintMetrics.HET_LOG10_CHI_SQUARED_PVALUE = Math.log10(hetsChiSquaredTest);
        fingerprintMetrics.HET_CROSS_ENTROPY_LOD = hetsCrossEntropy;
        fingerprintMetrics.HOM_CHI_SQUARED_PVALUE = homsChiSquaredTest;
        fingerprintMetrics.HOM_LOG10_CHI_SQUARED_PVALUE = Math.log10(homsChiSquaredTest);
        fingerprintMetrics.HOM_CROSS_ENTROPY_LOD = homsCrossEntropy;
        fingerprintMetrics.LOD_SELF_CHECK = lodSelfCheck;
        fingerprintMetrics.DISCRIMINATORY_POWER = lodSelfCheck - MathUtil.mean(randomizationTrials);

        return fingerprintMetrics;
    }

    /** Creates a new fingerprint from the current one by randomizing the genotypes
     *
     * @return
     */

    private Fingerprint randomizeFingerprint(final RandomGenerator rg) {
        final Fingerprint retVal = new Fingerprint(null, null, null);

        RandomDataGenerator rng = new RandomDataGenerator(rg);
        for (Map.Entry<HaplotypeBlock, HaplotypeProbabilities> entry : entrySet()) {
            final HaplotypeProbabilities haplotypeProbabilities = entry.getValue();
            final HaplotypeProbabilitiesFromGenotypeLikelihoods permutedHaplotypeProbabilities = new HaplotypeProbabilitiesFromGenotypeLikelihoods(entry.getKey());
            permutedHaplotypeProbabilities.addToLogLikelihoods(
                    haplotypeProbabilities.getRepresentativeSnp(),
                    haplotypeProbabilities.getRepresentativeSnp().getAlleles(),
                    MathUtil.permute(haplotypeProbabilities.getLogLikelihoods(), rng));
            retVal.add(permutedHaplotypeProbabilities);
        }
        return retVal;
    }


    /**
     * Attempts to filter out haplotypes that may have suspect genotyping by removing haplotypes that reach
     * a minimum confidence score yet have a significant fraction of observations from a third or fourth allele.
     */
    public void filterSuspectSites() {
        final Iterator<Map.Entry<HaplotypeBlock, HaplotypeProbabilities>> iterator = entrySet().iterator();
        while (iterator.hasNext()) {
            final Map.Entry<HaplotypeBlock, HaplotypeProbabilities> entry = iterator.next();
            final HaplotypeProbabilities p = entry.getValue();
            if (p instanceof HaplotypeProbabilitiesFromSequence) {
                final HaplotypeProbabilitiesFromSequence probs = (HaplotypeProbabilitiesFromSequence) p;

                if (probs.getLodMostProbableGenotype() >= 3 && probs.getFractionUnexpectedAlleleObs() > 0.1) {
                    iterator.remove();
                }
            }
        }
    }

    public static class FingerprintMetrics extends MetricBase {
        public String SAMPLE_NAME;
        public String SOURCE;
        public String INFO;
        public long HAPLOTYPE;
        public long HAPLOTYPES_WITH_EVIDENCE;
        public long DEFINITE_GENOTYPES;
        public long NUM_HOM_REF;
        public long NUM_HET;
        public long NUM_HOM_VAR;
        public double CHI_SQUARED_PVALUE;
        public double LOG10_CHI_SQUARED_PVALUE;
        public double CROSS_ENTROPY_LOD;
        public double HET_CHI_SQUARED_PVALUE;
        public double HET_LOG10_CHI_SQUARED_PVALUE;
        public double HET_CROSS_ENTROPY_LOD;
        public double HOM_CHI_SQUARED_PVALUE;
        public double HOM_LOG10_CHI_SQUARED_PVALUE;
        public double HOM_CROSS_ENTROPY_LOD;
        public double DISCRIMINATORY_POWER;
        public double LOD_SELF_CHECK;

    }
}
