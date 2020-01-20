/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
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

import htsjdk.samtools.util.Log;
import picard.pedigree.Sex;
import picard.util.MathUtil;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/**
 * Abstract class that encapsulates the logic behind determining sex from sequence data.
 * Extending classes will have to implement two functions, one that returns a list of SampleXy instances, and another that returns
 * the starting centroids for the k-means algorithm
 * The algorithm cluasters the data, looking for 2 clusters, corresponding to the two sexes.
 *
 * @author Yossi Farjoun
 */
abstract public class SexInferenceEngine {

    final protected Set<String> MALE_CHROMS = new TreeSet<>();
    final protected Set<String> FEMALE_CHROMS = new TreeSet<>();
    final protected Set<String> SEX_CHROMS = new TreeSet<>();

    final Log log = Log.getInstance(SexInferenceEngine.class);

    protected SexInferenceEngine(final Set<String> MALE_CHROMS, final Set<String> FEMALE_CHROMS) {

        this.MALE_CHROMS.addAll(MALE_CHROMS);
        this.FEMALE_CHROMS.addAll(FEMALE_CHROMS);

        this.SEX_CHROMS.addAll(FEMALE_CHROMS);
        this.SEX_CHROMS.addAll(MALE_CHROMS);

    }

    /**
     * Determine the sex of each of the samples.
     * performs a 2-cluster k-means clustering of those values to assign samples to either sex
     * based on the fact that females will have density ~ 0 on chrY and that the density on chrX that is approximately double
     * that of males.
     *
     * @return a Map of sample name to sex.
     */
    public Map<String, Sex> determineSexes() {
        final Map<String, Sex> sampleSexes = new HashMap<>();
        log.info("getting Chromosomal Coverage");
        final List<SampleXy> samples = getSexChromCoverageDensity();

        // Now cluster the coverage values.
        final double[][] data = new double[samples.size()][];
        for (int i = 0; i < data.length; ++i) {
            data[i] = new double[]{samples.get(i).xDensity, samples.get(i).yDensity};
        }
        log.info("Clustering sex data");
        final int[] assignments = MathUtil.kMeansCluster(data, getCentroids());

        // Lastly assign gender based on cluster membership and emit debug messaging about called genders
        final NumberFormat fmt = new DecimalFormat("00.000");
        log.info("Assigning sex to samples");
        for (int i = 0; i < assignments.length; ++i) {
            final SampleXy sample = samples.get(i);
            final Sex sex = assignments[i] == 0 ? Sex.Female : Sex.Male;
            sampleSexes.put(sample.sample, sex);
            log.debug(sample.sample + "\t" + data[i][0] + "\t" + data[i][1] + "\t" + sex);
        }
        return sampleSexes;
    }

    /**
     * Used by determineSexes
     *
     * @return A list of SampleXy's - objects containing sample names and a measure of coverage on the sex chromosome.
     */
    abstract protected List<SampleXy> getSexChromCoverageDensity();

    /**
     * Get coordinates of the points with which to begin k means clustering in determineSexes
     *
     * @return Double array
     */
    abstract protected double[][] getCentroids();

    /*A little class to encapsulate sample and the density of reads on each gender chromosome. */
    static protected class SampleXy {
        final public String sample;
        private double xDensity;
        private double yDensity;

        SampleXy(final String sample, final double x, final double y) {
            this.sample = sample;
            this.xDensity = x;
            this.yDensity = y;
        }

        void setxDensity(final double den) {
            xDensity = den;
        }

        void setyDensity(final double den) {
            yDensity = den;
        }

        public double getxDensity() {
            return xDensity;
        }

        public double getyDensity() {
            return yDensity;
        }
    }
}
