/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

package picard.analysis;

import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.QualityUtil;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import picard.PicardException;
import picard.util.MathUtil;

import java.util.*;
import java.util.stream.IntStream;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.distribution.BinomialDistribution;

/**
 * Created by David Benjamin on 5/13/15.
 */
public class TheoreticalSensitivity {

    private static final Log log = Log.getInstance(TheoreticalSensitivity.class);
    private static final int SAMPLING_MAX = 600; //prevent 'infinite' loops
    private static final int MAX_CONSIDERED_DEPTH_HET_SENS = 1000; // No point in looking any deeper than this, otherwise GC overhead is too high.  Only used for HET sensitivity.
    private static final double DEPTH_BIN_WIDTH = 0.01; // Minimal fraction of depth histogram to use when integrating theoretical sensitivity.  This ensures we don't calculate theoretical sensitivity at every depth, which would be computationally expensive.
    private static final int RANDOM_SEED = 51;

    /**
     * @param depthDistribution   the probability of depth n is depthDistribution[n] for n = 0, 1. . . N - 1
     * @param qualityDistribution the probability of quality q is qualityDistribution[q] for q = 0, 1. . . Q
     * @param sampleSize          sample size is the number of random sums of quality scores for each m
     * @param logOddsThreshold    is the log_10 of the likelihood ratio required to call a SNP,
     *                            for example 5 if the variant likelihood must be 10^5 times greater
     */
    public static double hetSNPSensitivity(final double[] depthDistribution, final double[] qualityDistribution,
                                           final int sampleSize, final double logOddsThreshold) {
        return hetSNPSensitivity(depthDistribution, qualityDistribution, sampleSize, logOddsThreshold, true);
    }

    /**
     * @param depthDistribution   the probability of depth n is depthDistribution[n] for n = 0, 1. . . N - 1
     * @param qualityDistribution the probability of quality q is qualityDistribution[q] for q = 0, 1. . . Q
     * @param sampleSize          sample size is the number of random sums of quality scores for each m
     * @param logOddsThreshold    is the log_10 of the likelihood ratio required to call a SNP,
     *                            for example 5 if the variant likelihood must be 10^5 times greater.
     * @param withLogging         true to output log messages, false otherwise.
     */
    public static double hetSNPSensitivity(final double[] depthDistribution, final double[] qualityDistribution,
                                           final int sampleSize, final double logOddsThreshold, final boolean withLogging) {
        final int N = Math.min(depthDistribution.length, MAX_CONSIDERED_DEPTH_HET_SENS + 1);

        if (withLogging) log.info("Creating Roulette Wheel");
        final RouletteWheel qualitySampler = new RouletteWheel(qualityDistribution);

        //qualitySums[m] is a random sample of sums of m quality scores, for m = 0, 1, N - 1
        if (withLogging) log.info("Calculating quality sums from quality sampler");
        final List<ArrayList<Integer>> qualitySums = qualitySampler.sampleCumulativeSums(N, sampleSize, withLogging);

        //if a quality sum of m qualities exceeds the quality sum threshold for n total reads, a SNP is called
        final ArrayList<Double> qualitySumThresholds = new ArrayList<>(N);
        final double LOG_10 = Math.log10(2);

        for (int n = 0; n < N; n++) qualitySumThresholds.add(10 * (n * LOG_10 + logOddsThreshold));

        //probabilityToExceedThreshold[m][n] is the probability that the sum of m quality score
        //exceeds the nth quality sum threshold
        if (withLogging) log.info("Calculating theoretical het sensitivity");
        final List<ArrayList<Double>> probabilityToExceedThreshold = proportionsAboveThresholds(qualitySums, qualitySumThresholds);
        final List<ArrayList<Double>> altDepthDistribution = hetAltDepthDistribution(N);
        double result = 0.0;
        for (int n = 0; n < N; n++) {
            for (int m = 0; m <= n; m++) {
                result += depthDistribution[n] * altDepthDistribution.get(n).get(m) * probabilityToExceedThreshold.get(m).get(n);
            }
        }
        return result;
    }

    //given L lists of lists and N thresholds, count the proportion of each list above each threshold
    public static List<ArrayList<Double>> proportionsAboveThresholds(final List<ArrayList<Integer>> lists, final List<Double> thresholds) {
        final ArrayList<ArrayList<Double>> result = new ArrayList<>();

        for (final ArrayList<Integer> list : lists) {
            final ArrayList<Double> newRow = new ArrayList<>(Collections.nCopies(thresholds.size(), 0.0));
            Collections.sort(list);
            int n = 0;
            int j = 0;  //index within the ordered sample
            while (n < thresholds.size() && j < list.size()) {
                if (thresholds.get(n) > list.get(j)) j++;
                else newRow.set(n++, (double) (list.size() - j) / list.size());
            }
            result.add(newRow);
        }
        return result;
    }

    //Utility function for making table of binomial distribution probabilities nCm * (0.5)^n
    //for n = 0, 1 . . . N - 1 and m = 0, 1. . . n
    public static List<ArrayList<Double>> hetAltDepthDistribution(final int N) {
        final List<ArrayList<Double>> table = new ArrayList<>();
        for (int n = 0; n < N; n++) {
            final ArrayList<Double> nthRow = new ArrayList<>();

            //add the 0th element, then elements 1 through n - 1, then the nth.
            //Note that nCm = (n-1)C(m-1) * (n/m)
            nthRow.add(Math.pow(0.5, n));
            for (int m = 1; m < n; m++) nthRow.add((n * 0.5 / m) * table.get(n - 1).get(m - 1));
            if (n > 0) nthRow.add(nthRow.get(0));

            table.add(nthRow);
        }
        return table;
    }

    /*
    Perform random draws from {0, 1. . . N - 1} according to a list of relative probabilities.

    We use an O(1) stochastic acceptance algorithm -- see Physica A, Volume 391, Page 2193 (2012) --
    which works well when the ratio of maximum weight to average weight is not large.
     */
    public static class RouletteWheel {
        final private List<Double> probabilities;
        final private int N;
        private int count = 0;
        private Random rng;

        RouletteWheel(final double[] weights) {
            rng = new Random(RANDOM_SEED);
            N = weights.length;

            probabilities = new ArrayList<>();
            final double wMax = MathUtil.max(weights);

            if (wMax == 0) {
                throw new PicardException("Quality score distribution is empty.");
            }

            for (final double w : weights) {
                probabilities.add(w / wMax);
            }
        }

        public int draw() {
            while (true) {
                final int n = (int) (N * rng.nextDouble());
                count++;
                if (rng.nextDouble() < probabilities.get(n)) {
                    count = 0;
                    return n;
                } else if (count >= SAMPLING_MAX) {
                    count = 0;
                    return 0;
                }
            }
        }

        //get samples of sums of 0, 1, 2,. . .  N - 1 draws
        public List<ArrayList<Integer>> sampleCumulativeSums(final int maxNumberOfSummands, final int sampleSize, final boolean withLogging) {
            final List<ArrayList<Integer>> result = new ArrayList<>();
            for (int m = 0; m < maxNumberOfSummands; m++) result.add(new ArrayList<>());

            for (int iteration = 0; iteration < sampleSize; iteration++) {
                int cumulativeSum = 0;
                for (int m = 0; m < maxNumberOfSummands; m++) {
                    result.get(m).add(cumulativeSum);
                    cumulativeSum += draw();
                }
                if (withLogging && iteration % 1000 == 0) {
                    log.info(iteration + " sampling iterations completed");
                }
            }
            return result;
        }
    }

    public static double[] normalizeHistogram(final Histogram<Integer> histogram) {
        if (histogram == null) throw new PicardException("Histogram is null and cannot be normalized");

        final double histogramSumOfValues = histogram.getSumOfValues();
        final double[] normalizedHistogram = new double[(int) histogram.getMax() + 1];

        final Iterator<Integer> keySet = histogram.keySet().iterator();

        for (Iterator<Integer> it = keySet; it.hasNext(); ) {
            Integer key = it.next();
            normalizedHistogram[key] = histogram.get(key).getValue() / histogramSumOfValues;
        }

        double s = 0;
        for(int i = 0;i < (int) histogram.getMax() + 1;i++) {
            s += normalizedHistogram[i];
        }

        return normalizedHistogram;
    }

    /**
     * Determines if a variant would be called under the particular conditions of a given total depth, alt depth,
     * average base qualities, allele fraction of variant and log odds threshold necessary to call variant.
     * @param totalDepth Depth at the site to be called, both alt and ref.
     * @param altDepth Number of alt bases at this site.
     * @param sumOfAltQualities Average Phred-scaled quality of bases
     * @param alleleFraction Allele fraction we are attempting to detect
     * @param logOddsThreshold Log odds threshold necessary for variant to be called
     * @return Whether or not the model would call a variant given the parameters
     */
    @VisibleForTesting
     static boolean isCalled(final int totalDepth, final int altDepth, final double sumOfAltQualities, final double alleleFraction, final double logOddsThreshold) {
        final double threshold = 10.0 * (altDepth * -Math.log10(alleleFraction) + (totalDepth - altDepth) * -Math.log10(1.0 - alleleFraction) + logOddsThreshold);

        return sumOfAltQualities > threshold;
    }

    public TheoreticalSensitivity() {
    }

    /**
     * Calculates the theoretical sensitivity with a given Phred-scaled quality score distribution at a constant
     * depth.
     * @param depth Depth to compute sensitivity at
     * @param qualityHistogram Phred-scaled quality score histogram
     * @param logOddsThreshold Log odd threshold necessary for variant to be called
     * @param sampleSize sampleSize is the total number of simulations to run
     * @param alleleFraction the allele fraction to evaluate sensitivity at
     * @param randomSeed random number seed to use for random number generator
     * @param overlapProbability Probability a position on a read overlaps with its mate
     * @param maxOverlapBaseQuality Maximum base quality that can result from summing overlapping base qualities, often capped by PCR error rate
     * @return Theoretical sensitivity for the given arguments at a constant depth.
     */
    public static double sensitivityAtConstantDepth(final int depth, final Histogram<Integer> qualityHistogram,
                                                    final double logOddsThreshold, final int sampleSize,
                                                    final double alleleFraction, final long randomSeed,
                                                    final double overlapProbability, final int maxOverlapBaseQuality) {
        // If the depth is 0 at a particular locus, the sensitivity is trivially 0.0.
        if (depth == 0) {
            return 0.0;
        }

        final RouletteWheel qualityRW = new RouletteWheel(trimDistribution(normalizeHistogram(qualityHistogram)));
        final Random randomNumberGenerator = new Random(randomSeed);
        final RandomGenerator rg = new Well19937c(randomSeed);
        final BinomialDistribution bd = new BinomialDistribution(rg, depth, alleleFraction);

        // Calculate mean and deviation of quality score distribution to enable Gaussian sampling below
        final double averageQuality = qualityHistogram.getMean();
        final double standardDeviationQuality = qualityHistogram.getStandardDeviation();

        int calledVariants = 0;
        // Sample simulated variants, and count the number that would get called.  The ratio
        // of the number called to the total sampleSize is the sensitivity.
        for (int sample = 0; sample < sampleSize; sample++) {
            final int altDepth = bd.sample();

            int sumOfQualities;
            // If the number of alt reads is "small" we draw from the actual base quality distribution.
            sumOfQualities = 0;
            for(int i = 0;i < altDepth;i++) {
                if(rg.nextDouble() > overlapProbability) {
                    sumOfQualities += Math.min(qualityRW.draw(), maxOverlapBaseQuality);
                } else {
                    // Simulate overlapping reads by sampling from the quality distribution twice.
                    // Since PCR errors affect overlapping reads, the total contribution to the sum
                    // of qualities cannot exceed the PCR error rate.
                    sumOfQualities += Math.min(qualityRW.draw() + qualityRW.draw(), maxOverlapBaseQuality);
                }
            }


            if (isCalled(depth, altDepth, (double) sumOfQualities, alleleFraction, logOddsThreshold)) {
                calledVariants++;
            }
            //System.out.println("AF = " + alleleFraction + " altDepth = " + altDepth + " " + sumOfQualities + " " + isCalled(depth, altDepth, (double) sumOfQualities, alleleFraction, logOddsThreshold));
        }
        return (double) calledVariants / sampleSize;
    }

    /**
     * Calculates the theoretical sensitivity with a given Phred-scaled quality score distribution at a constant
     * depth.
     * @param depth Depth to compute sensitivity at
     * @param qualityHistogram Phred-scaled quality score histogram
     * @param logOddsThreshold Log odds threshold necessary for variant to be called
     * @param sampleSize the total number of simulations to run
     * @param alleleFraction the allele fraction to evaluate sensitivity at
     * @return Theoretical sensitivity for the given arguments at a constant depth.
     */
    private static double sensitivityAtConstantDepth(final int depth, final Histogram<Integer> qualityHistogram,
                                                     final double logOddsThreshold, final int sampleSize,
                                                     final double alleleFraction, final double overlapProbability,
                                                     final int pcrErrorRate) {
        return sensitivityAtConstantDepth(depth, qualityHistogram, logOddsThreshold, sampleSize, alleleFraction,
                RANDOM_SEED, overlapProbability, pcrErrorRate);
    }

    /**
     * Calculates the theoretical sensitivity with a given Phred-scaled quality score distribution and depth
     * distribution.
     * @param depthHistogram Depth histogram to compute theoretical sensitivity over
     * @param qualityHistogram Phred-scaled quality score histogram
     * @param sampleSize the total number of simulations to run
     * @param logOddsThreshold Log odds threshold necessary for variant to be called
     * @param alleleFraction the allele fraction to evaluate sensitivity at
     * @param overlapProbability the probability two bases from the same fragment are overlapping
     * @param pcrErrorRate Upper bound for the effective base quality two bases in an overlapping read
     * @return Theoretical sensitivity for the given arguments over a particular depth distribution.
     */
    public static double theoreticalSensitivity(final Histogram<Integer> depthHistogram,
                                                final Histogram<Integer> qualityHistogram, final int sampleSize,
                                                final double logOddsThreshold, final double alleleFraction,
                                                final double overlapProbability, final int pcrErrorRate) {
        return theoreticalSensitivity(depthHistogram, qualityHistogram, sampleSize, logOddsThreshold,
           alleleFraction, overlapProbability, pcrErrorRate, false);
    }

    public static double theoreticalSensitivity(final Histogram<Integer> depthHistogram,
                                                final Histogram<Integer> qualityHistogram, final int sampleSize,
                                                final double logOddsThreshold, final double alleleFraction,
                                                final double overlapProbability, final int pcrErrorRate,
                                                final boolean useTrapezoidRule) {
        if (alleleFraction > 1.0 || alleleFraction < 0.0) {
            throw new IllegalArgumentException("Allele fractions must be between 0 and 1.");
        }

        final double[] depthDistribution = normalizeHistogram(depthHistogram);

        // Integrate sensitivity over depth distribution
        double sensitivity = 0.0;
        int currentDepth = 0;

        if(useTrapezoidRule) {
            double right = 0;
            while (currentDepth < depthDistribution.length) {
                double deltaDepthProbability = 0.0;
                // Accumulate a portion of the depth distribution to compute theoretical sensitivity over.
                // This helps prevent us from spending lots of compute over coverages
                // that occur with low probability and don't contribute much to sensitivity anyway, but
                // it complicates things a bit by having a variable deltaDepthProbability which
                // amount of the depth distribution to use with the trapezoid rule integration step.
                while (deltaDepthProbability == 0 && currentDepth < depthDistribution.length ||
                        deltaDepthProbability < DEPTH_BIN_WIDTH && currentDepth < depthDistribution.length &&
                                depthDistribution[currentDepth] < DEPTH_BIN_WIDTH / 2.0) {
                    deltaDepthProbability += depthDistribution[currentDepth];
                    currentDepth++;
                }
                // I think we overshoot by 1


                // Calculate sensitivity for a particular depth, and use trapezoid rule to integrate sensitivity.
                final double left = right;

                System.out.println("Current Depth = " + currentDepth);
                right = sensitivityAtConstantDepth(currentDepth, qualityHistogram, logOddsThreshold, sampleSize,
                        alleleFraction, overlapProbability, pcrErrorRate);
                sensitivity += deltaDepthProbability * (left + right) / 2.0;
            }
        }
        else {
            for(currentDepth = 0;currentDepth < depthDistribution.length;currentDepth++) {
                if(depthDistribution[currentDepth] != 0) {
                    System.out.println("Depth = " + currentDepth);
                    sensitivity += sensitivityAtConstantDepth(currentDepth, qualityHistogram, logOddsThreshold, sampleSize,
                            alleleFraction, overlapProbability, pcrErrorRate) * depthDistribution[currentDepth];
                }
            }
        }
        return sensitivity;
    }

    /**
     * Removes trailing zeros in a distribution.  The purpose of this function is to prevent other
     * functions from evaluating in regions where the distribution has zero probability.
     * If the array consists of only zeros or is empty, an empty array is returned.
     *
     * @param distribution Distribution of base qualities
     * @return Distribution of base qualities removing any trailing zeros
     */
    static double[] trimDistribution(final double[] distribution) {
        int endOfDistribution = distribution.length;
        if (distribution.length == 0) {
            return distribution;
        }

        // Locate last element in array that is non-zero if it exists.
        // If there are no non-zero elements, endOfDistribution should be 0.
        while (distribution[endOfDistribution - 1] == 0) {
            endOfDistribution--;
            if (endOfDistribution == 0) break;
        }

        // Remove trailing zeros and return.
        return Arrays.copyOfRange(distribution, 0, endOfDistribution);
    }

    /**
     * This is a utility function to calculate the metrics specific to running
     * theoretical sensitivity over several different allele fractions.
     * @param simulationSize Number of simulations to run at each depth.
     * @param depthHistogram Histogram of depth distribution.
     * @param baseQHistogram Histogram of Phred-scaled quality scores.
     * @param alleleFractions List of allele fractions to measure theoretical sensitivity over.
     */
    public static List<TheoreticalSensitivityMetrics> calculateSensitivities(final int simulationSize,
                                                                             final Histogram<Integer> depthHistogram,
                                                                             final Histogram<Integer> baseQHistogram,
                                                                             final List<Double> alleleFractions,
                                                                             final double overlapProbability,
                                                                             final int pcrErrorRate) {

        final List<TheoreticalSensitivityMetrics> metricsOverVariousAlleleFractions = new ArrayList<>();
        final double logOddsThreshold = 6.2; // This threshold is used because it is the value used for MuTect2.

        // For each allele fraction in alleleFractions calculate theoretical sensitivity and add the results
        // to the histogram sensitivityHistogram.
        for (final double alleleFraction : alleleFractions) {
            final TheoreticalSensitivityMetrics theoreticalSensitivityMetrics = new TheoreticalSensitivityMetrics();
            theoreticalSensitivityMetrics.ALLELE_FRACTION = alleleFraction;
            theoreticalSensitivityMetrics.THEORETICAL_SENSITIVITY = TheoreticalSensitivity.theoreticalSensitivity(
                    depthHistogram, baseQHistogram, simulationSize, logOddsThreshold, alleleFraction, overlapProbability, pcrErrorRate);
            theoreticalSensitivityMetrics.THEORETICAL_SENSITIVITY_Q = QualityUtil.getPhredScoreFromErrorProbability((1 - theoreticalSensitivityMetrics.THEORETICAL_SENSITIVITY));
            theoreticalSensitivityMetrics.SAMPLE_SIZE = simulationSize;
            theoreticalSensitivityMetrics.LOG_ODDS_THRESHOLD = logOddsThreshold;
            metricsOverVariousAlleleFractions.add(theoreticalSensitivityMetrics);
        }

        return metricsOverVariousAlleleFractions;
    }
}
