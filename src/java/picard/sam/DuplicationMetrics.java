/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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

package picard.sam;

import htsjdk.samtools.metrics.MetricsFile;
import org.apache.commons.jexl2.internal.ArrayListWrapper;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.util.FastMath;
import static picard.util.MathUtil.pNormalizeVector;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.util.Histogram;

import java.io.*;
import java.util.Collection;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

/**
 * Metrics that are calculated during the process of marking duplicates
 * within a stream of SAMRecords.
 */
public class DuplicationMetrics extends MetricBase {
    /** The library on which the duplicate marking was performed. */
    public String LIBRARY;

    /**
     * The number of mapped reads examined which did not have a mapped mate pair,
     * either because the read is unpaired, or the read is paired to an unmapped mate.
     */
    public long UNPAIRED_READS_EXAMINED;

    /** The number of mapped read pairs examined. (Primary, non-supplemental) */
    public long READ_PAIRS_EXAMINED;

    /** The number of reads that were either secondary or supplementary */
    public long SECONDARY_OR_SUPPLEMENTARY_RDS;

    /** The total number of unmapped reads examined. (Primary, non-supplemental) */
    public long UNMAPPED_READS;

    /** The number of fragments that were marked as duplicates. */
    public long UNPAIRED_READ_DUPLICATES;

    /** The number of read pairs that were marked as duplicates. */
    public long READ_PAIR_DUPLICATES;

    /**
     * The number of read pairs duplicates that were caused by optical duplication.
     * Value is always < READ_PAIR_DUPLICATES, which counts all duplicates regardless of source.
     */
    public long READ_PAIR_OPTICAL_DUPLICATES;

    /** The percentage of mapped sequence that is marked as duplicate. */
    public Double PERCENT_DUPLICATION;

    /** The estimated number of unique molecules in the library based on PE duplication. */
    public Long ESTIMATED_LIBRARY_SIZE;

    /** Grid search matrix for parameter optimization */
    public double[][] GRID_SEARCH_MTRX;

    /** Polya draw search interval */
    public double[] POLYA_DRAW_SEARCH_INT;

    /** Library size multiple search interval */
    public double[] LIBRARY_SIZE_MULTIPLE_SEARCH_INT;

    /** Interval of additional read pairs multiples for ROI projections */
    public double[] SEQUENCING_MULTIPLE_INTERVAL;

    /** Matrix of projected occupancies based on sequencing multiple */
    public double[][] PROJECTED_OCCUPANCY_MATRIX;

    /** projection of unique molecules based on sequencing multiple */
    public double[] UNIQUE_MOLECULE_COUNT_PROJECTION;

    /** optimization results */
    public double[] OPTIMIZATION_RESULTS;

    /**
     * Fills in the ESTIMATED_LIBRARY_SIZE based on the paired read data examined where
     * possible and the PERCENT_DUPLICATION.
     */
    public void calculateDerivedMetrics() {
        this.ESTIMATED_LIBRARY_SIZE = estimateLibrarySize(this.READ_PAIRS_EXAMINED - this.READ_PAIR_OPTICAL_DUPLICATES,
                this.READ_PAIRS_EXAMINED - this.READ_PAIR_DUPLICATES);

        PERCENT_DUPLICATION = (UNPAIRED_READ_DUPLICATES + READ_PAIR_DUPLICATES *2) /(double) (UNPAIRED_READS_EXAMINED + READ_PAIRS_EXAMINED *2);
    }

    /**
     * Estimates the size of a library based on the number of paired end molecules observed
     * and the number of unique pairs observed.
     *
     * Based on the Lander-Waterman equation that states:
     *     C/X = 1 - exp( -N/X )
     * where
     *     X = number of distinct molecules in library
     *     N = number of read pairs
     *     C = number of distinct fragments observed in read pairs
     */
    public static Long estimateLibrarySize(final long readPairs, final long uniqueReadPairs) {
        final long readPairDuplicates = readPairs - uniqueReadPairs;

        if (readPairs > 0 && readPairDuplicates > 0) {
            long n = readPairs;
            long c = uniqueReadPairs;

            double m = 1.0, M = 100.0;

            if (c >= n || f(m*c, c, n) < 0) {
                throw new IllegalStateException("Invalid values for pairs and unique pairs: "
                        + n + ", " + c);

            }

            while( f(M*c, c, n) >= 0 ) M *= 10.0;

            for (int i=0; i<40; i++ ) {
                double r = (m+M)/2.0;
                double u = f( r * c, c, n );
                if ( u == 0 ) break;
                else if ( u > 0 ) m = r;
                else if ( u < 0 ) M = r;
            }

            return (long) (c * (m+M)/2.0);
        }
        else {
            return null;
        }
    }

    /** Method that is used in the computation of estimated library size. */
    private static double f(double x, double c, double n) {
        return c/x - 1 + Math.exp(-n/x);
    }

    /**
     * Estimates the ROI (return on investment) that one would see if a library was sequenced to
     * x higher coverage than the observed coverage.
     *
     * @param estimatedLibrarySize the estimated number of molecules in the library
     * @param x the multiple of sequencing to be simulated (i.e. how many X sequencing)
     * @param pairs the number of pairs observed in the actual sequencing
     * @param uniquePairs the number of unique pairs observed in the actual sequencing
     * @return a number z <= x that estimates if you had pairs*x as your sequencing then you
     *         would observe uniquePairs*z unique pairs.
     */
    public static double estimateRoi(long estimatedLibrarySize, double x, long pairs, long uniquePairs) {
        return estimatedLibrarySize * ( 1 - Math.exp(-(x*pairs)/estimatedLibrarySize) ) / uniquePairs;
    }

    /**
     * Calculates a histogram using the estimateRoi method to estimate the effective yield
     * doing x sequencing for x=1..10.
     */
    public Histogram<Double> calculateRoiHistogram() {
        if (ESTIMATED_LIBRARY_SIZE == null) {
            try {
                calculateDerivedMetrics();
                if (ESTIMATED_LIBRARY_SIZE == null) return null;
            }
            catch (IllegalStateException ise) { return null; }
        }

        long uniquePairs = READ_PAIRS_EXAMINED - READ_PAIR_DUPLICATES;
        Histogram<Double> histo = new Histogram<Double>();

        for (double x=1; x<=100; x+=1) {
            histo.increment(x, estimateRoi(ESTIMATED_LIBRARY_SIZE, x, READ_PAIRS_EXAMINED, uniquePairs));
        }

        return histo;
    }

    public void setSequencingMultipleInterval(List<String> passedInterval) {
        this.SEQUENCING_MULTIPLE_INTERVAL = new double[passedInterval.size()];
        for (int i=0; i<passedInterval.size(); i++){
            this.SEQUENCING_MULTIPLE_INTERVAL[i] = Double.parseDouble(passedInterval.get(i));
        }
    }

    public void getInputHistFromMetricsFile (java.io.File INPUT_METRICS_FILE) {

        final MetricsFile<DuplicationMetrics, Double> metricsOutput = new MetricsFile<DuplicationMetrics, Double>();
        try{
            int g = 1;
            metricsOutput.read(new FileReader(INPUT_METRICS_FILE));
        }
        catch (final FileNotFoundException ex) {
            System.err.println("Metrics file not found: " + ex);
        }

        
        for (final Histogram<Double> histo : metricsOutput.getAllHistograms()) {
            String label = histo.getValueLabel();
            if (label.equals("non_optical_pairs")) {
                //histToArray(histo);
                getNonOptDupHistogram(histo);
            }
        }

    }

    public void histToArray(Histogram<Double> histo) {
        final double histoCount = histo.getCount();
        Collection<Histogram.Bin> vals = histo.values();
        int[] countArray = new int[vals.size()];
        int[] sizeArray = new int[vals.size()];
        int count = 0;
        HashMap<Integer, Integer> hm = new HashMap();
        for (final Histogram<Double>.Bin bin : histo.values()) {
            countArray[count] = (int) bin.getValue();
            sizeArray[count] = (int) bin.getIdValue();
            hm.put(sizeArray[count], countArray[count]);
            count += 1;
        }
        int g=1;
    }

    public Histogram<Double> getNonOptDupHistogram(Histogram<Double> histo) {
        final double histoCount = histo.getCount();
        Collection<Histogram.Bin> vals = histo.values();
        int[] countArray = new int[vals.size()];
        int[] sizeArray = new int[vals.size()];
        int count = 0;
        HashMap<Integer, Integer> hm = new HashMap();
        for (final Histogram<Double>.Bin bin : histo.values()) {
            countArray[count] = (int) bin.getValue();
            sizeArray[count] = (int) bin.getIdValue();
            hm.put(sizeArray[count], countArray[count]);
            count += 1;
        }
        /*int[] countArray = {1990,75,3,1,1,0,0};
        int[] sizeArray = {1,2,3,4,5,6,7};
        HashMap<Integer, Integer> hm = new HashMap();
        for (int i=0; i<countArray.length; i++){
            hm.put(sizeArray[i], countArray[i]);
        }*/
        // values for test array
        int maxSize = Arrays.stream(sizeArray).max().getAsInt();
        // extend the array to contain 0 counts for entries with zero occupancy
        double[] observedOccupancy = new double[maxSize];
        for (int i=0; i<maxSize; i++) {
            if (hm.containsKey(i+1)) {
                observedOccupancy[i] = (double) hm.get(i + 1);
            } else {
                observedOccupancy[i] = 0;
            }
        }
        int m = 0; // non-optical pairs observed
        int u = 0;
        for (int k=0; k<maxSize; k++){
            m += (k+1)*observedOccupancy[k];
            u += observedOccupancy[k];
        }
        double[] normCntVec = pNormalizeVector(observedOccupancy);
        // coarse grid search
        double[] seqMultipleRange1 = makeInterval(1,6,.2,true);
        double[] drawRange1 = makeInterval(3,20,1,true);
        System.out.println("completing coarse grid search:"+ seqMultipleRange1.length + "x" +drawRange1.length);
        double[][] paramMtrx = performGridSearchforKLMin(normCntVec,m,seqMultipleRange1,drawRange1);
        int[] coordMinCoarse = getMatrixMin(paramMtrx);
        //refined grid search
        double cgLibMin = seqMultipleRange1[coordMinCoarse[0]];
        double cgDrawMin = drawRange1[coordMinCoarse[1]];
        this.LIBRARY_SIZE_MULTIPLE_SEARCH_INT = makeInterval(.5*cgLibMin,2*cgLibMin,cgLibMin/20,false);
        this.POLYA_DRAW_SEARCH_INT = makeInterval(FastMath.log(cgDrawMin)-2,FastMath.log(cgDrawMin)+1.2,.1,true);
        System.out.println("completing refined grid search:"+ this.LIBRARY_SIZE_MULTIPLE_SEARCH_INT.length + "x" + this.POLYA_DRAW_SEARCH_INT.length);
        this.GRID_SEARCH_MTRX = performGridSearchforKLMin(normCntVec,
                m,
                this.LIBRARY_SIZE_MULTIPLE_SEARCH_INT,
                this.POLYA_DRAW_SEARCH_INT);
        int[] coordMin = getMatrixMin(this.GRID_SEARCH_MTRX);
        int librarySizeEstimate = (int) (m*this.LIBRARY_SIZE_MULTIPLE_SEARCH_INT[coordMin[0]]);
        int urnDrawEstimate = (int) this.POLYA_DRAW_SEARCH_INT[coordMin[1]];
        System.out.println("indices:"+ coordMin[0] + ":" +coordMin[1] +
                ", librarySizeEstimate = " + librarySizeEstimate +
                ", urnDrawEstimate = " + urnDrawEstimate);
        double[] sampledOccupancyEstimate = estimateSampledLibraryOccupancy(librarySizeEstimate,
                urnDrawEstimate,m);
        this.OPTIMIZATION_RESULTS = new double[5];
        this.OPTIMIZATION_RESULTS[0]=coordMin[0];
        this.OPTIMIZATION_RESULTS[1]=coordMin[1];
        this.OPTIMIZATION_RESULTS[2]=librarySizeEstimate;
        this.OPTIMIZATION_RESULTS[3]=urnDrawEstimate;
        this.OPTIMIZATION_RESULTS[4]=m;
        System.out.println(Arrays.toString(sampledOccupancyEstimate));
        System.out.println(Arrays.toString(normCntVec));
        // make projections
        //this.SEQUENCING_MULTIPLE_INTERVAL = makeInterval(0,10,1,false);
        double[] estimatedFractUnique = makeOccupancyProjections(m,librarySizeEstimate,urnDrawEstimate,this.SEQUENCING_MULTIPLE_INTERVAL);
        this.UNIQUE_MOLECULE_COUNT_PROJECTION = new double[estimatedFractUnique.length];
        for (int i=0; i<estimatedFractUnique.length; i++){
            this.UNIQUE_MOLECULE_COUNT_PROJECTION[i] = m*this.SEQUENCING_MULTIPLE_INTERVAL[i]*estimatedFractUnique[i];
        }
        return histo;
    }

    public double[] makeOccupancyProjections(int m, int librarySizeEstimate,
                                             int urnDrawEstimate,
                                             double[] seqMultipleInterval) {
        double[] sampledOccupancyEstimate = estimateSampledLibraryOccupancy(librarySizeEstimate,
                urnDrawEstimate,m);
        double[][] probProjMtrx = new double[seqMultipleInterval.length][sampledOccupancyEstimate.length]; // matrix representing occupancy probabilities
        this.PROJECTED_OCCUPANCY_MATRIX = new double[seqMultipleInterval.length][sampledOccupancyEstimate.length];
        double[] estimatedFractUnique = new double[seqMultipleInterval.length];
        for (int i=0; i<seqMultipleInterval.length; i++){
            probProjMtrx[i] = estimateSampledLibraryOccupancy(librarySizeEstimate,
                    urnDrawEstimate, (int) (m*seqMultipleInterval[i]));
            double tDup = 0;
            double tTot = 0;
            for (int j=0; j<sampledOccupancyEstimate.length; j++){
                tDup += probProjMtrx[i][j]*(j);
                tTot += probProjMtrx[i][j]*(j+1);
            }
            estimatedFractUnique[i] = 1-(tDup/tTot); // fraction unique
            for (int j=0; j<probProjMtrx[i].length; j++) {
                this.PROJECTED_OCCUPANCY_MATRIX[i][j] = seqMultipleInterval[i]*m*(1-(tDup/tTot))*probProjMtrx[i][j];
            }
        }
        return estimatedFractUnique;
    }


    public double[] estimateLibraryOccupancy(int N, int n, int m){
        /**
         * This function estimates the theoretical library occupancy before sampling under a set of Polya paramaters
         * returns a vector of probabilities for k=1...40
         * @param N Library size estimate
         * @param n number of polya draws
         * @param m number of observed samples
         */
        double[] libEstimate = new double[40];
        double betaDen = Beta.logBeta(1, (double) N-1);
        for (int k=0; k<40; k++ ) {
            double nCk = CombinatoricsUtils.binomialCoefficientLog(n,k);
            double betaNum = Beta.logBeta(k+1, (n-k)+N-1);
            double lgProb = nCk+betaNum-betaDen; // log sum of probabilities
            libEstimate[k] = Math.exp(lgProb);
            int g=1;
        }
        double libSum = 0;
        for (double x : libEstimate)libSum += x;
        return libEstimate;
    }

    public double[] estimateSampledLibraryOccupancy(int N, int n, int m){
        /**
         * This function estimates the library occupancy after sampling under a set of Polya paramaters. The
         * model assumes sampling with replacement under binomial sampling.
         * returns a vector of probabilities for k=1...40
         */
        double[] libOccupancyEstimate = estimateLibraryOccupancy(N,n,m);
        double[] sampledOccupancyEstimate = new double[libOccupancyEstimate.length];
        double M = 0;
        for (int k=0; k<(libOccupancyEstimate.length); k++){
            M += N*(k+1)*libOccupancyEstimate[k];
        }
        for (int k=0; k<(libOccupancyEstimate.length); k++ ) {
            double[] t = new double[libOccupancyEstimate.length];
            double tSum = 0;
            for (int j=0; j<(libOccupancyEstimate.length); j++ ) {
                double p = ((double) j / M);
                if (p>1){
                    p = 1;
                }
                BinomialDistribution binom = new BinomialDistribution(m,p);
                t[j] = libOccupancyEstimate[j]*binom.probability(k);
                tSum += t[j];
            }
            sampledOccupancyEstimate[k] = tSum;
        }
        double sampSum = 0;
        for (double x : sampledOccupancyEstimate)sampSum += x;
        return sampledOccupancyEstimate;
    }

    public double[] makeInterval(double start, double stop, double step, boolean makeExponential) {
        double k = start;
        int nIter = (int) ((int) (stop-start)/step);
        double[] intval = new double[nIter];
        for(int i = 0; i < nIter; i++){
            k += step;
            if (makeExponential) {
                intval[i] = Math.exp(k);
            }
            else {
                intval[i] = k;
            }

        }
        return intval;
    }

    public double[][] performGridSearchforKLMin(double [] observedOccupancy, int m,
                                                double[] seqMultipleRange, double[] drawRange){
        /**
         * Perform grid search to find minimum of KL divergence
         */
        double[][] paramMtrx = new double[seqMultipleRange.length][drawRange.length];
        for (int i = 0; i<seqMultipleRange.length; i++){
            for (int j = 0; j<drawRange.length; j++){
                paramMtrx[i][j] = occupancyKLDivergence(observedOccupancy, (int) drawRange[j],
                        (int) (m*seqMultipleRange[i]), m, observedOccupancy.length);
            }
        }
        return paramMtrx;
    }

    public double occupancyKLDivergence(double [] observedOccupancy,int n, int N, int m, int maxK){
        /**
         * calculate KL divergence between observed occupancy and expected occupancy
         * @param observedOccupancy vector of observed occupancy
         * @param N Library size estimate
         * @param n number of polya draws
         * @param m number of observed samples
         * @param maxK maximum value of k to which KL will be calculated
         */
        double[] sampledOccupancyEstimate = estimateSampledLibraryOccupancy(N,n,m);
        int ix; // max index for evaluation
        if (sampledOccupancyEstimate.length < maxK) {
            ix = sampledOccupancyEstimate.length;
        }
        else {
            ix = maxK;
        }
        double klDiverge = 0;
        for (int i=0; i<ix; i++){
            klDiverge += observedOccupancy[i] * FastMath.log(observedOccupancy[i]/sampledOccupancyEstimate[i]);
        }
        return klDiverge;
    }

    public int[] getMatrixMin (double[][] A) {
        double minVal = A[0][0];
        int[] minCoord = {0,0};
        for ( int i = 0; i < A.length; i++ ){
            for ( int j = 0; j < A[i].length; j++ ) {
                if (A[i][j] < minVal) {
                    minVal = A[i][j];
                    minCoord[0] = i;
                    minCoord[1] = j;
                }
            }
        }
        return minCoord;
    }

    public void writeMatricesToFile (String filename) {
        try {
            BufferedWriter outputWriter = null;
            outputWriter = new BufferedWriter(new FileWriter(filename,true));
            // write unique reads result matrix
            /*int nLines = 30;  // arbitrary lines to write
            for (int i = 0; i < nLines; ++i) {
                outputWriter.newLine();
            }*/
            outputWriter.write("Polya Urn Model - KL divergence grid search (Library size multiple x number of urn draws)");
            outputWriter.newLine();
            outputWriter.write("-, ");
            outputWriter.write(Arrays.toString(this.POLYA_DRAW_SEARCH_INT)
                    .replace("[","")
                    .replace("]",""));
            outputWriter.newLine();
            for (int i = 0; i < this.GRID_SEARCH_MTRX.length; i++) {
                outputWriter.write(this.LIBRARY_SIZE_MULTIPLE_SEARCH_INT[i] + ", ");
                outputWriter.write(Arrays.toString(this.GRID_SEARCH_MTRX[i])
                        .replace("[","")
                        .replace("]",""));
                outputWriter.newLine();
            }
            // write projection of unique molecules
            outputWriter.newLine();
            outputWriter.write("Optimization results");
            outputWriter.newLine();
            outputWriter.write("-, ");
            outputWriter.write(Arrays.toString(this.OPTIMIZATION_RESULTS)
                    .replace("[","")
                    .replace("]",""));
            // write projection matrix
            outputWriter.newLine();
            outputWriter.write("Count projections (sequencing multiple x occupancy)");
            outputWriter.newLine();
            outputWriter.write("-, ");
            outputWriter.write(Arrays.toString(makeInterval(0,this.PROJECTED_OCCUPANCY_MATRIX[0].length,1,false))
                    .replace("[","")
                    .replace("]",""));
            outputWriter.newLine();
            for (int i = 0; i < this.PROJECTED_OCCUPANCY_MATRIX.length; i++) {
                outputWriter.write(this.SEQUENCING_MULTIPLE_INTERVAL[i] + ", ");
                outputWriter.write(Arrays.toString(this.PROJECTED_OCCUPANCY_MATRIX[i])
                        .replace("[","")
                        .replace("]",""));
                outputWriter.newLine();
            }
            // write projection of unique molecules
            outputWriter.newLine();
            outputWriter.write("unique molecule projections");
            outputWriter.newLine();
            outputWriter.write("sequencing-multiples, ");
            outputWriter.write(Arrays.toString(this.SEQUENCING_MULTIPLE_INTERVAL)
                    .replace("[","")
                    .replace("]",""));
            outputWriter.newLine();
            outputWriter.write("unique-molecule-projections, ");
            outputWriter.write(Arrays.toString(this.UNIQUE_MOLECULE_COUNT_PROJECTION)
                    .replace("[","")
                    .replace("]",""));
            outputWriter.flush();
            outputWriter.close();
        } catch (final IOException e) {
            System.out.println("problem writing result matrix");
        }
    }

    // Main method used for debugging the derived metrics
    // Usage = DuplicationMetrics READ_PAIRS READ_PAIR_DUPLICATES
    public static void main(String[] args) {
        DuplicationMetrics m = new DuplicationMetrics();
        m.READ_PAIRS_EXAMINED  = Integer.parseInt(args[0]);
        m.READ_PAIR_DUPLICATES = Integer.parseInt(args[1]);
        m.calculateDerivedMetrics();
        System.out.println("Percent Duplication: " + m.PERCENT_DUPLICATION);
        System.out.println("Est. Library Size  : " + m.ESTIMATED_LIBRARY_SIZE);
        System.out.println();

        System.out.println("X Seq\tX Unique");
        for (Histogram<Double>.Bin bin : m.calculateRoiHistogram().values()) {
            System.out.println(bin.getId() + "\t" + bin.getValue());
        }

//        DuplicationMetrics m = new DuplicationMetrics();
//        m.READ_PAIRS_EXAMINED  = Long.parseLong(args[0]);
//        m.READ_PAIR_DUPLICATES = Long.parseLong(args[1]);
//        final long UNIQUE_READ_PAIRS = m.READ_PAIRS_EXAMINED - m.READ_PAIR_DUPLICATES;
//        final double xCoverage = Double.parseDouble(args[2]);
//        final double uniqueXCoverage = xCoverage * ((double) UNIQUE_READ_PAIRS / (double) m.READ_PAIRS_EXAMINED);
//        final double oneOverCoverage = 1 / xCoverage;
//
//        m.calculateDerivedMetrics();
//        System.out.println("Percent Duplication: " + m.PERCENT_DUPLICATION);
//        System.out.println("Est. Library Size  : " + m.ESTIMATED_LIBRARY_SIZE);
//        System.out.println();
//
//
//        System.out.println("Coverage\tUnique Coverage\tDuplication");
//        for (double d = oneOverCoverage; (int) (d*xCoverage)<=50; d+=oneOverCoverage) {
//            double coverage = d * xCoverage;
//            double uniqueCoverage = uniqueXCoverage * m.estimateRoi(m.ESTIMATED_LIBRARY_SIZE, d, m.READ_PAIRS_EXAMINED, UNIQUE_READ_PAIRS);
//            double duplication = (coverage - uniqueCoverage) / coverage;
//            System.out.println(coverage + "\t" + uniqueCoverage + "\t" + duplication);
//        }
    }
}
