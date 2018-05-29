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

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.util.*;
import java.util.stream.IntStream;

/**
 * Created by davidben on 5/18/15.
 */
public class TheoreticalSensitivityTest {

    private final static File TEST_DIR = new File("testdata/picard/analysis/TheoreticalSensitivity/");
    private final static File DEPTH = new File(TEST_DIR, "Solexa332667_DepthDist.histo");
    private final static File BASEQ = new File(TEST_DIR, "Solexa332667_BaseQ.histo");

    @Test
    public void testRouletteWheel() throws Exception {

        //test that a deterministic roulette wheel only gives one value
        final double[] deterministicWeights = {0.0, 1.0, 0.0};
        final TheoreticalSensitivity.RouletteWheel deterministicWheel = new TheoreticalSensitivity.RouletteWheel(deterministicWeights);
        for (int n = 0; n < 10; n++) Assert.assertEquals(deterministicWheel.draw(), 1);

        //test the sums of this deterministic wheel: a sum of n 1's equals n
        final List<ArrayList<Integer>> deterministicSums = deterministicWheel.sampleCumulativeSums(10, 1, false);
        for (int n = 0; n < 10; n++) Assert.assertEquals(deterministicSums.get(n).get(0), (Integer) n);
    }

    @Test
    public void testProportionsAboveThresholds() throws Exception {
        final List<ArrayList<Integer>> sums = new ArrayList<>();
        sums.add(new ArrayList<>(Arrays.asList(0, 0, 0)));
        sums.add(new ArrayList<>(Arrays.asList(10, 10)));
        sums.add(new ArrayList<>(Arrays.asList(5, 11, -2, 4)));
        final List<Double> thresholds = Arrays.asList(-1.0, 1.0, 6.0);
        Assert.assertEquals(sums.size(), 3);
        Assert.assertEquals(thresholds.size(), 3);

        final List<ArrayList<Double>> proportions = TheoreticalSensitivity.proportionsAboveThresholds(sums, thresholds);
        Assert.assertEquals(proportions.size(), 3);

        Assert.assertEquals(proportions.get(0).get(0), (double) 3 / 3);
        Assert.assertEquals(proportions.get(0).get(1), (double) 0 / 3);
        Assert.assertEquals(proportions.get(0).get(2), (double) 0 / 3);
        Assert.assertEquals(proportions.get(1).get(0), (double) 2 / 2);
        Assert.assertEquals(proportions.get(1).get(1), (double) 2 / 2);
        Assert.assertEquals(proportions.get(1).get(2), (double) 2 / 2);
        Assert.assertEquals(proportions.get(2).get(0), (double) 3 / 4);
        Assert.assertEquals(proportions.get(2).get(1), (double) 3 / 4);
        Assert.assertEquals(proportions.get(2).get(2), (double) 1 / 4);
    }

    @Test
    public void testHetAltDepthDistribution() throws Exception {
        final int N = 6;
        final double p = 0.5;
        final List<ArrayList<Double>> distribution = TheoreticalSensitivity.hetAltDepthDistribution(N);

        for (int n = 0; n < N - 1; n++) {
            for (int m = 0; m <= n; m++) {

                final long binomialCoefficient = CombinatoricsUtils.binomialCoefficient(n, m);

                Assert.assertEquals(distribution.get(n).get(m), binomialCoefficient * Math.pow(p, n));
            }
        }
    }

    //test that large-sample sums from the RouletteWheel converge to a normal distribution
    //using the empirical CDF as measured by proportionsAboveThresholds
    @Test
    public void testCentralLimitTheorem() throws Exception {
        //use a RouletteWheel that gives 0, 1, 2 with equal probability
        final double[] weights = {1.0, 1.0, 1.0};
        final TheoreticalSensitivity.RouletteWheel wheel = new TheoreticalSensitivity.RouletteWheel(weights);

        final int sampleSize = 1000;
        final int numSummands = 100;

        //the mean and standard deviation of a single roulette draw and of many draws
        final double muSingleDraw = 1.0;
        final double sigmaSingleDraw = Math.sqrt(2.0 / 3.0);
        final double mu = numSummands * muSingleDraw;
        final double sigma = Math.sqrt(numSummands) * sigmaSingleDraw;

        //test the sums of this deterministic wheel: a sum of n 1's equals n
        final List<ArrayList<Integer>> sums = wheel.sampleCumulativeSums(numSummands, sampleSize, false);
        //we only want the last set of sums, those with numSummands summands
        sums.subList(0, sums.size() - 1).clear();

        Assert.assertEquals(sums.size(), 1);

        //test whether the number of elements within one standard deviation agrees with the normal distribution
        final List<Double> thresholds = Arrays.asList(mu - sigma, mu + sigma);

        //sums is 1 x sampleSize, thresholds is a 2-vector, so proportions is 1 x 2
        final List<ArrayList<Double>> proportions = TheoreticalSensitivity.proportionsAboveThresholds(sums, thresholds);
        final double empiricalProportionWithinOneSigma = proportions.get(0).get(0) - proportions.get(0).get(1);

        //the proportion within one sigma for the normal distribution
        //hence whether any element falls within one sigma is a Bernoulli variable
        final double theoreticalProportionWithinOneSigma = 0.682689492;
        final double samplingStandardDeviationOfProportion = Math.sqrt(theoreticalProportionWithinOneSigma * (1 - theoreticalProportionWithinOneSigma) / sampleSize);

        Assert.assertEquals(empiricalProportionWithinOneSigma, theoreticalProportionWithinOneSigma, 5 * samplingStandardDeviationOfProportion);
    }

    //Put it all together for deterministic quality and depths
    @Test
    public void testDeterministicQualityAndDepth() throws Exception {
        final double logOddsThreshold = 0.0;
        final double tolerance = 0.001;
        final int sampleSize = 1; //quality is deterministic, hence no sampling error
        for (int q = 5; q < 10; q++) {
            for (int n = 5; n < 10; n++) {
                final double minAltCount = 10 * n * Math.log10(2) / q;  //alts required to call when log odds ratio threshold = 1
                double expectedResult = 0.0;

                final List<ArrayList<Double>> altCountProbabilities = TheoreticalSensitivity.hetAltDepthDistribution(n + 1);
                for (int altCount = n; altCount > minAltCount; altCount--) {
                    expectedResult += altCountProbabilities.get(n).get(altCount);
                }

                //deterministic weights that always yield q are 0.0 for 0 through q - 1 and 1.0 for q
                final double[] qualityDistribution = new double[q + 1];
                Arrays.fill(qualityDistribution, 0L);
                qualityDistribution[qualityDistribution.length - 1] = 1L;
                final double[] depthDistribution = new double[n + 1];
                Arrays.fill(depthDistribution, 0L);
                depthDistribution[depthDistribution.length - 1] = 1L;

                final double result = TheoreticalSensitivity.hetSNPSensitivity(depthDistribution, qualityDistribution, sampleSize, logOddsThreshold);
                Assert.assertEquals(result, expectedResult, tolerance);
            }
        }
    }

    @Test
    public void testHetSensDistributions() throws Exception {
        //Expect theoretical sens to be close to .9617 for Solexa-332667
        final double tolerance = 0.02;
        final double expectedResult = .9617;
        final int maxDepth = 500;
        final double[] depthDistribution = new double[maxDepth + 1];
        final double[] qualityDistribution = new double[50];

        final Scanner scanDepth = new Scanner(DEPTH);
        for (int i = 0; scanDepth.hasNextDouble(); i++) {
            depthDistribution[i] = scanDepth.nextDouble();
        }
        final Scanner scanBaseQ = new Scanner(BASEQ);
        for (int j = 0; scanBaseQ.hasNextDouble(); j++) {
            qualityDistribution[j] = scanBaseQ.nextDouble();
        }

        final int sampleSize = 1_000;
        final double logOddsThreshold = 3.0;
        final double result = TheoreticalSensitivity.hetSNPSensitivity(depthDistribution, qualityDistribution, sampleSize, logOddsThreshold);
        Assert.assertEquals(result, expectedResult, tolerance);
    }

    @DataProvider(name = "hetSensDataProvider")
    public Object[][] hetSensDataProvider() {
        final File wgsMetricsFile = new File(TEST_DIR, "test_Solexa-332667.wgs_metrics");
        final File targetedMetricsFile = new File(TEST_DIR, "test_25103070136.targeted_pcr_metrics");

        //These magic numbers come from a separate implementation of the code in R.
        return new Object[][]{
                {0.897_342_54, wgsMetricsFile},
                {0.956_186_66, targetedMetricsFile}
        };
    }

    @Test(dataProvider = "hetSensDataProvider")
    public void testHetSensTargeted(final double expected, final File metricsFile) throws Exception {
        final double tolerance = 0.000_000_01;

        final MetricsFile<?, Integer> metrics = new MetricsFile<>();
        try (final FileReader metricsFileReader = new FileReader(metricsFile)) {
            metrics.read(metricsFileReader);
        }

        final List<Histogram<Integer>> histograms = metrics.getAllHistograms();
        final Histogram<Integer> depthHistogram = histograms.get(0);
        final Histogram<Integer> qualityHistogram = histograms.get(1);

        final double[] depthDistribution = TheoreticalSensitivity.normalizeHistogram(depthHistogram);
        final double[] qualityDistribution = TheoreticalSensitivity.normalizeHistogram(qualityHistogram);

        final int sampleSize = 1_000;
        final double logOddsThreshold = 3.0;

        final double result = TheoreticalSensitivity.hetSNPSensitivity(depthDistribution, qualityDistribution, sampleSize, logOddsThreshold);
        Assert.assertEquals(result, expected, tolerance);
    }

    @DataProvider(name = "TheoreticalSensitivityConstantDepthDataProvider")
    public Object[][] fractionalAlleleSensDataProvider() {
        final File wgsMetricsFile = new File(TEST_DIR, "test_Solexa-332667.wgs_metrics");
        final File targetedMetricsFile = new File(TEST_DIR, "test_25103070136.targeted_pcr_metrics");

        // These magic numbers come from a separate implementation of the code in R.
        return new Object[][]{
                // Expected sensitivity, metrics file, allele fraction, constant depth, sample size.
                {1.00, wgsMetricsFile, .5, 30, 10000, 0.01},
                {0.78, targetedMetricsFile, .1, 30, 10000, 0.02},
                {0.26, targetedMetricsFile, 0.1, 10, 10000, 0.01}
        };
    }

    @Test(dataProvider = "TheoreticalSensitivityConstantDepthDataProvider")
    public void testSensitivityAtConstantDepth(final double expected, final File metricsFile, final double alleleFraction, final int depth, final int sampleSize, final double tolerance) throws Exception {
        // This tests Theoretical Sensitivity assuming a uniform depth with a distribution of base quality scores.
        // Because this only tests sensitivity at a constant depth, we use this for testing the code at high depths.
        final MetricsFile<?, Integer> metrics = new MetricsFile<>();
        try (final FileReader metricsFileReader = new FileReader(metricsFile)) {
            metrics.read(metricsFileReader);
        }

        final List<Histogram<Integer>> histograms = metrics.getAllHistograms();
        final Histogram<Integer> qualityHistogram = histograms.get(1);

        // We ensure that even using different random seeds we converge to roughly the same value.
        for (int i = 0; i < 3; i++) {
            double result = TheoreticalSensitivity.sensitivityAtConstantDepth(depth, qualityHistogram, 3, sampleSize, alleleFraction, i);
            Assert.assertEquals(result, expected, tolerance);
        }
    }

    @DataProvider(name = "TheoreticalSensitivityDataProvider")
    public Object[][] arbFracSensDataProvider() {
        final File wgsMetricsFile = new File(TEST_DIR, "test_Solexa-332667.wgs_metrics");

        // This test acts primarily as an integration test.  The sample sizes
        // are not quite large enough to converge properly, but is used for the purpose of
        // keeping the compute time of the tests short.
        return new Object[][]{
                {0.90, wgsMetricsFile, 0.5, 400},
                {0.77, wgsMetricsFile, 0.3, 400},
                {0.29, wgsMetricsFile, 0.1, 500},
                {0.08, wgsMetricsFile, 0.05, 500},
        };
    }

    @Test(dataProvider = "TheoreticalSensitivityDataProvider")
    public void testSensitivity(final double expected, final File metricsFile, final double alleleFraction, final int sampleSize) throws Exception {
        // This tests Theoretical Sensitivity using distributions on both base quality scores
        // and the depth histogram.

        // We use a pretty forgiving tolerance here because for these tests
        // we are not using large enough sample sizes to converge.
        final double tolerance = 0.02;

        final MetricsFile<?, Integer> metrics = new MetricsFile<>();
        try (final FileReader metricsFileReader = new FileReader(metricsFile)) {
            metrics.read(metricsFileReader);
        }

        final List<Histogram<Integer>> histograms = metrics.getAllHistograms();
        final Histogram<Integer> depthHistogram = histograms.get(0);
        final Histogram<Integer> qualityHistogram = histograms.get(1);

        final double result = TheoreticalSensitivity.theoreticalSensitivity(depthHistogram, qualityHistogram, sampleSize, 3, alleleFraction);

        Assert.assertEquals(result, expected, tolerance);
    }

    @DataProvider(name = "equivalanceHetVsArbitrary")
    public Object[][] equivalenceHetVsFull() {
        final File wgsMetricsFile = new File(TEST_DIR, "test_Solexa-332667.wgs_metrics");
        final File targetedMetricsFile = new File(TEST_DIR, "test_25103070136.targeted_pcr_metrics");

        return new Object[][]{
                // The sample sizes chosen here for these tests are smaller than what would normally be used
                // in order to keep the test time low.  It should be noted that for larger sample sizes
                // the values converge.
                {wgsMetricsFile, 0.02, 500},
                {targetedMetricsFile, 0.01, 500}
        };
    }

    @Test(dataProvider = "equivalanceHetVsArbitrary")
    public void testHetVsArbitrary(final File metricsFile, final double tolerance, final int sampleSize) throws Exception {
        // This test compares Theoretical Sensitivity for arbitrary allele fractions with the theoretical het sensitivity
        // model.  Since allele fraction of 0.5 is equivalent to a het, these should provide the same answer.
        final MetricsFile<?, Integer> metrics = new MetricsFile<>();
        try (final FileReader metricsFileReader = new FileReader(metricsFile)) {
            metrics.read(metricsFileReader);
        }

        final List<Histogram<Integer>> histograms = metrics.getAllHistograms();
        final Histogram<Integer> depthHistogram = histograms.get(0);
        final Histogram<Integer> qualityHistogram = histograms.get(1);

        final double[] qualityDistribution = TheoreticalSensitivity.normalizeHistogram(qualityHistogram);
        final double[] depthDistribution = TheoreticalSensitivity.normalizeHistogram(depthHistogram);

        final double resultFromTS = TheoreticalSensitivity.theoreticalSensitivity(depthHistogram, qualityHistogram, sampleSize, 3, 0.5);
        final double resultFromTHS = TheoreticalSensitivity.hetSNPSensitivity(depthDistribution, qualityDistribution, sampleSize, 3);

        Assert.assertEquals(resultFromTS, resultFromTHS, tolerance);
    }

    @DataProvider(name = "callingThresholdDataProvider")
    public Object[][] callingThreshold() {
        return new Object[][]{
                // These values were tested with an independent implementation in R.
                // Test a transition due to a change in the logOddsThreshold
                {100, 10, 10 * 20, .1, 5.8, true},
                {100, 10, 10 * 20, .1, 5.9, false},

                // Test a transition due to change in average base quality from 20 to 21
                {100, 10, 10 * 21, .1, 6.2, true},
                {100, 10, 10 * 20, .1, 6.2, false},

                // Test a transition due to change in total depth
                {115, 10, 10 * 21, .1, 6.2, false},
                {114, 10, 10 * 21, .1, 6.2, true}
        };
    }

    @Test(dataProvider = "callingThresholdDataProvider")
    public void testCallingThreshold(final int totalDepth, final int altDepth, final double sumOfAltQualities, final double alleleFraction, final double logOddsThreshold, final boolean expectedCall) {
        Assert.assertEquals(TheoreticalSensitivity.isCalled(totalDepth, altDepth, sumOfAltQualities, alleleFraction, logOddsThreshold), expectedCall);
    }

    @DataProvider(name = "sumOfGaussiansDataProvider")
    public Object[][] sumOfGaussians() {
        final File wgsMetricsFile = new File(TEST_DIR, "test_Solexa-332667.wgs_metrics");
        final File targetedMetricsFile = new File(TEST_DIR, "test_25103070136.targeted_pcr_metrics");

        // When we sum more base qualities from a particular distribution, it should look increasingly Gaussian.
        return new Object[][]{
                {wgsMetricsFile, 500, 0.03},
                {wgsMetricsFile, 20, 0.05},
                {wgsMetricsFile, 10, 0.10},
                {targetedMetricsFile, 500, 0.03},
                {targetedMetricsFile, 20, 0.05},
                {targetedMetricsFile, 10, 0.10}
        };
    }

    @Test(dataProvider = "sumOfGaussiansDataProvider")
    public void testDrawSumOfQScores(final File metricsFile, final int altDepth, final double tolerance) throws Exception {
        final MetricsFile<TheoreticalSensitivityMetrics, Integer> metrics = new MetricsFile<>();
        try (final FileReader metricsFileReader = new FileReader(metricsFile)) {
            metrics.read(metricsFileReader);
        }

        final List<Histogram<Integer>> histograms = metrics.getAllHistograms();

        final Histogram<Integer> qualityHistogram = histograms.get(1);
        final TheoreticalSensitivity.RouletteWheel qualityRW = new TheoreticalSensitivity.RouletteWheel(TheoreticalSensitivity.trimDistribution(TheoreticalSensitivity.normalizeHistogram(qualityHistogram)));

        final Random randomNumberGenerator = new Random(51);

        // Calculate mean and deviation of quality score distribution to enable Gaussian sampling below
        final double averageQuality = qualityHistogram.getMean();
        final double standardDeviationQuality = qualityHistogram.getStandardDeviation();

        for (int k = 0; k < 1; k++) {
            int sumOfQualitiesFull = IntStream.range(0, altDepth).map(n -> qualityRW.draw()).sum();
            int sumOfQualities = TheoreticalSensitivity.drawSumOfQScores(altDepth, averageQuality, standardDeviationQuality, randomNumberGenerator.nextGaussian());

            Assert.assertEquals(sumOfQualitiesFull, sumOfQualities, sumOfQualitiesFull * tolerance);
        }
    }
}
