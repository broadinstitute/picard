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

import htsjdk.tribble.util.MathUtils;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import picard.PicardException;
import picard.util.MathUtil;

import java.nio.file.Path;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.stream.Collectors;

import static picard.fingerprint.HaplotypeProbabilities.Genotype.HET_ALLELE12;
import static picard.fingerprint.HaplotypeProbabilities.Genotype.HOM_ALLELE1;
import static picard.fingerprint.HaplotypeProbabilities.Genotype.HOM_ALLELE2;

/**
 * class to represent a genetic fingerprint as a set of HaplotypeProbabilities
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


        for (final HaplotypeProbabilities haplotypeProbabilities : this.values()) {
            genotypeCounts = MathUtil.sum(genotypeCounts, haplotypeProbabilities.getPosteriorProbabilities());
            expectedRatios = MathUtil.sum(expectedRatios, haplotypeProbabilities.getPriorProbablities());
        }
        final long[] roundedGenotypeCounts = MathUtil.round(genotypeCounts);

        final double[] hetVsHomExpect = new double[]{expectedRatios[HOM_ALLELE1.v] + expectedRatios[HOM_ALLELE2.v], expectedRatios[HET_ALLELE12.v]};
        final double[] hetVsHomCounts = new double[]{genotypeCounts[HOM_ALLELE1.v] + genotypeCounts[HOM_ALLELE2.v], genotypeCounts[HET_ALLELE12.v]};

        final double[] homAllele1VsAllele2Expect = new double[]{expectedRatios[HOM_ALLELE1.v], expectedRatios[HOM_ALLELE2.v]};
        final double[] homAllele1VsAllele2Counts = new double[]{genotypeCounts[HOM_ALLELE1.v], genotypeCounts[HOM_ALLELE2.v]};

        final long[] roundedHomRefVsVarCounts = MathUtil.round(homAllele1VsAllele2Counts);
        final long[] roundedHetVsHomCounts = MathUtil.round(hetVsHomCounts);

        // calculate p-value
        final ChiSquareTest chiSquareTest = new ChiSquareTest();
        final double chiSquaredTest = chiSquareTest.chiSquareTest(expectedRatios, roundedGenotypeCounts);
        // calculate LOD (cross-entropy)
        final double crossEntropy = -MathUtil.klDivergance(genotypeCounts, expectedRatios);

        // calculate p-value
        final double hetsChiSquaredTest = chiSquareTest.chiSquareTest(hetVsHomExpect, roundedHetVsHomCounts);

        // calculate LOD (cross-entropy)
        final double hetsCrossEntropy = -MathUtil.klDivergance(hetVsHomCounts, hetVsHomExpect);

        // calculate p-value
        final double homsChiSquaredTest = chiSquareTest.chiSquareTest(homAllele1VsAllele2Expect, roundedHomRefVsVarCounts);

        // calculate LOD (cross-entropy)
        final double homsCrossEntropy = -MathUtil.klDivergance(homAllele1VsAllele2Counts, homAllele1VsAllele2Expect);

        final double lodSelfCheck = FingerprintChecker.calculateMatchResults(this, this).getLOD();

        final double[] randomizationTrials = new double[NUMBER_OF_SAMPLING];
        final RandomGenerator rg = new MersenneTwister(RANDOM_SEED);

        MathUtils.RunningStat runningStat= new MathUtils.RunningStat();

        for (int i = 0; i < NUMBER_OF_SAMPLING; i++) {
            runningStat.push(FingerprintChecker.calculateMatchResults(this, randomizeFingerprint(rg)).getLOD());
        }

        FingerprintMetrics fingerprintMetrics = new FingerprintMetrics();

        fingerprintMetrics.SAMPLE_ALIAS = sample;
        fingerprintMetrics.SOURCE = source.toUri().toString();
        fingerprintMetrics.INFO = info;
        fingerprintMetrics.HAPLOTYPES = values().size();
        fingerprintMetrics.HAPLOTYPES_WITH_EVIDENCE = values().stream().filter(HaplotypeProbabilities::hasEvidence).count();
        fingerprintMetrics.DEFINITE_GENOTYPES = values().stream().filter(h -> h.getLodMostProbableGenotype() > GENOTYPE_LOD_THRESHOLD).count();
        fingerprintMetrics.NUM_HOM_ALLELE1 = roundedGenotypeCounts[0];
        fingerprintMetrics.NUM_HET = roundedGenotypeCounts[1];
        fingerprintMetrics.NUM_HOM_ALLELE2 = roundedGenotypeCounts[2];
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

    /** Creates a new fingerprint from the current one by randomizing the probabilities within each haplotype
     *
     * @return
     */
    private Fingerprint randomizeFingerprint(final RandomGenerator rg) {
        final Fingerprint retVal = new Fingerprint(null, null, null);

        final RandomDataGenerator rng = new RandomDataGenerator(rg);
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

    public static Function<FingerprintIdDetails, String> getFingerprintIdDetailsStringFunction(CrosscheckMetric.DataType CROSSCHECK_BY) {
        final Function<FingerprintIdDetails, String> groupByTemp;
        switch (CROSSCHECK_BY) {
            case READGROUP:
                groupByTemp = details -> details.platformUnit;
                break;
            case LIBRARY:
                groupByTemp = details -> details.sample + "::" + details.library;
                break;
            case FILE:
                groupByTemp = details -> details.file + "::" + details.sample;
                break;
            case SAMPLE:
                groupByTemp = details -> details.sample;
                break;
            default:
                throw new PicardException("unpossible");
        }

        // if the groupBy string is null (e.g. a vcf file has no read group info) then the hashcode is
        // used intending to be unique per object (ignoring possible collisions)
        return key -> {
            final String temp = groupByTemp.apply(key);
            return temp == null ? Integer.toString(key.hashCode()) : temp;
        };
    }

    public static Map<FingerprintIdDetails, Fingerprint> mergeFingerprintsBy(
            final Map<FingerprintIdDetails, Fingerprint> fingerprints,
            final Function<FingerprintIdDetails, String> by) {

        // collect the various entries according to the grouping "by"

        final Map<String, List<Map.Entry<FingerprintIdDetails, Fingerprint>>> collection =
                fingerprints.entrySet()
                        .stream()
                        .collect(Collectors.groupingBy(entry -> by.apply(entry.getKey())));

        return collection.entrySet().stream()
                .collect(Collectors.toMap(
                        entry -> {
                            // merge the keys (unequal values are eliminated by merge).

                            final FingerprintIdDetails finalId = new FingerprintIdDetails();
                            entry.getValue().forEach(id -> finalId.merge(id.getKey()));
                            finalId.group = entry.getKey();
                            return finalId;

                        }, entry -> {
                            // merge the values by merging the fingerprints.

                            final FingerprintIdDetails firstDetail = entry.getValue().get(0).getKey();
                            //use the "by" function to determine the "info" part of the fingerprint
                            final Fingerprint sampleFp = new Fingerprint(firstDetail.sample, null, by.apply(firstDetail));
                            entry.getValue().stream().map(Map.Entry::getValue).collect(Collectors.toSet()).forEach(sampleFp::merge);
                            return sampleFp;

                        }));
    }

    enum CrosscheckMode implements CommandLineParser.ClpEnum {
        CHECK_SAME_SAMPLE {
            @Override
            public String getHelpDoc() {
                return "In this mode, each sample in INPUT will only be checked against a single corresponding sample in SECOND_INPUT. " +
                        "If a corresponding sample cannot be found, the program will proceed, but report the missing samples" +
                        " and return the value specified in EXIT_CODE_WHEN_MISMATCH. The corresponding samples are those that equal each other, after possible renaming " +
                        "via INPUT_SAMPLE_MAP and SECOND_INPUT_SAMPLE_MAP. In this mode CROSSCHECK_BY must be SAMPLE.";
            }
        },
        CHECK_ALL_OTHERS {
            @Override
            public String getHelpDoc() {
                return "In this mode, each sample in INPUT will be checked against all the samples in SECOND_INPUT.";
            }
        }
    }
}
