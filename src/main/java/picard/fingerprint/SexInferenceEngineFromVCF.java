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

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Class that implements  SexInferencer for an input type of VCF. The main issue here is that we want to extract average depth information from a VCF.
 * This is done by looking at a more-or-less equally spaced set of variants, and collecting statistics about the coverage that they have.
 */
public class SexInferenceEngineFromVCF extends SexInferenceEngine {

    private final OverlapDetector<Interval> parDetector = new OverlapDetector<>(0, 0);
    private final File vcf;
    private final int autosomalVariants;
    private final int allosomalVariants;
    private final ProgressLogger progressLogger;

    /**
     * Some VCFs do not have variants (called) on Y. Consequently, when determining sex, we wish to give less weight to the Y Coverage than
     * the X coverage. We do this by the scaling the Y coverage values by a factor of yCoverageShrinkingFactor.
     */
    final double yCoverageShrinkFactor;

    public SexInferenceEngineFromVCF(final Set<String> MALE_CHROMS, final Set<String> FEMALE_CHROMS, final File vcf, final int autosomalVariants, final int allosomalVariants, final File pseudoautosomalIntervalList, final double yCoverageShrinkFactor) {
        super(MALE_CHROMS, FEMALE_CHROMS);

        this.yCoverageShrinkFactor = yCoverageShrinkFactor;
        this.allosomalVariants = allosomalVariants;
        this.autosomalVariants = autosomalVariants;
        this.vcf = vcf;

        progressLogger = new ProgressLogger(log, 10000, "examined", "variants");

        final IntervalList intervalList;
        if (pseudoautosomalIntervalList != null) {
            IOUtil.assertFileIsReadable(pseudoautosomalIntervalList);
            intervalList = IntervalList.fromFile(pseudoautosomalIntervalList);
            parDetector.addAll(intervalList.getIntervals(), intervalList.getIntervals());
        }
    }

    // Check that the variant does lie on the pseudoautosomal region.
    // Used to filter out variants for coverage sampling.
    private boolean isPseudoautosomal(final String chrom, final int pos) {
        // skip the overlap detection if not within a sex chromosome
        return SEX_CHROMS.contains(chrom) &&
                !parDetector.getOverlaps(new Interval(chrom, pos, pos)).isEmpty();
    }

    /**
     * A class used to enumerate chromosome "types" - male, female, and non-sex - and collect coverage sampling stats on each type.
     */
    private enum chromStats {
        MALE_CHROM,
        FEMALE_CHROM,
        NON_SEX_CHROME;

        private int numVariantsUsed;
        private Map<String, Long> coverageStats;

        public static void initializeChromStats(final List<String> sampleNames) {
            for (final chromStats type : chromStats.values()) {
                type.coverageStats = new HashMap<>();
                type.numVariantsUsed = 0;
                for (final String sampleName : sampleNames) {
                    type.coverageStats.put(sampleName, 0L);
                }
            }
        }

        public void incrementVariantsUsed() {
            numVariantsUsed++;
        }

        public void addCoverage(final String sampleName, final int n) {
            coverageStats.put(sampleName, coverageStats.get(sampleName) + n);
        }

        /**
         * Calculate Average Coverage Per Variant (on the calling chromosome type)
         *
         * @return double
         */
        public double getAcpv(final String sampleName) {
            return numVariantsUsed == 0 ? 0 : coverageStats.get(sampleName) / (double) numVariantsUsed;
        }

    }

    private chromStats getChromType(final String chromName) {
        if (MALE_CHROMS.contains(chromName)) {
            return chromStats.MALE_CHROM;
        }
        if (FEMALE_CHROMS.contains(chromName)) {
            return chromStats.FEMALE_CHROM;
        }

        return chromStats.NON_SEX_CHROME;
    }

    private boolean isSexChrom(final String chromName) {
        return MALE_CHROMS.contains(chromName) || FEMALE_CHROMS.contains(chromName);
    }

    /**
     * @return A list of XySample object, one for each sample, containing the ratio
     * (average coverage on the female chromosomes)/(average coverage off the sex chromosomes),
     * and the same for male chromosomes.
     */
    @Override
    protected List<SampleXy> getSexChromCoverageDensity() {

        final List<SampleXy> samples = new ArrayList<>();

        final VCFFileReader reader = new VCFFileReader(vcf);
        final VCFHeader header = reader.getFileHeader();
        final SAMSequenceDictionary dict = reader.getFileHeader().getSequenceDictionary();
        final List<String> sampleNames = header.getGenotypeSamples();

        // initialize sampling coverage statistics.

        chromStats.initializeChromStats(sampleNames);

        // a map from true/false to size of ALLO/AUTO somal regions in the genome
        final Map<Boolean, Long> RegionSize = new HashMap<>(2);
        RegionSize.put(Boolean.TRUE, 0L);
        RegionSize.put(Boolean.FALSE, 0L);

        // collect total length of autosomal and allosomal regions
        for (final SAMSequenceRecord chrom : dict.getSequences()) {
            final Boolean isSex = isSexChrom(chrom.getSequenceName());
            RegionSize.put(isSex, chrom.getSequenceLength() + RegionSize.get(isSex));
        }

        for (final SAMSequenceRecord chrom : dict.getSequences()) {
            if (isSexChrom(chrom.getSequenceName())) {
                sampleDPs(chrom, (int) (RegionSize.get(Boolean.TRUE) / allosomalVariants), reader);
            } else {
                sampleDPs(chrom, (int) (RegionSize.get(Boolean.FALSE) / autosomalVariants), reader);
            }
        }

        reader.close();

        /*
         * Populate the list of sampleXy's.
         */
        for (final String sampleName : sampleNames) {
            final SampleXy xy = new SampleXy(sampleName,
                    chromStats.FEMALE_CHROM.getAcpv(sampleName) / chromStats.NON_SEX_CHROME.getAcpv(sampleName),
                    chromStats.MALE_CHROM.getAcpv(sampleName) / (yCoverageShrinkFactor * chromStats.NON_SEX_CHROME.getAcpv(sampleName)));
            samples.add(xy);
        }

        return samples;
    }

    /**
     * Sample coverage depth on a given chromosome at approximately evenly spaced intervals. This method updates covStats.
     *
     * @param chrom     The name of the chromosome where the estimate is sought.
     * @param blockSize The size of the block to use in the traversal
     * @param reader    VCFFilereader of the VCF.
     */
    private void sampleDPs(final SAMSequenceRecord chrom, final int blockSize, final VCFFileReader reader) {
        final int chromLength = chrom.getSequenceLength();

        int pos = 1;
        while (pos < chromLength) {
            progressLogger.record(chrom.getSequenceName(), pos);

            if (isPseudoautosomal(chrom.getSequenceName(), pos)) {
                //skip pseudoautosomal region
                final Interval parRegion = parDetector.getOverlaps(new Interval(chrom.getSequenceName(), pos, pos)).iterator().next();
                pos += parRegion.length();
            } else {
                //not pseudoautosomal
                final CloseableIterator<VariantContext> it = reader.query(chrom.getSequenceName(), pos, Math.min(pos + blockSize, chromLength));
                if (it.hasNext()) {
                    final VariantContext gc = it.next();
                    processVariantDepth(gc);
                }
                pos += blockSize;
            }
        }
    }

    /**
     * Extracts the coverage depth over a given context for every sample in the vcf. Updates covStats
     *
     * @param context a VariantContext
     */
    private void processVariantDepth(final VariantContext context) {
        final chromStats type = getChromType(context.getChr());
        type.incrementVariantsUsed();
        for (final Genotype g : context.getGenotypes()) {
            final String name = g.getSampleName();
            type.addCoverage(name, Math.max(0, g.getDP()));
        }
    }

    /**
     * Centroids from which to begin sex chromosome coverage density k means clustering for sex determination.
     * The starting centroid for the male cluster reflects the fact that we expect to not always have variants on the Y chromosome,
     * thus we give the coverage on Y less weight than the coverage on X (see yCoverageShrinkingFactor) when determining sex.
     * Consequently, we expect yDensities to be between 0 and 0.25, and we choose the centroid to be (0.5,0.125) (instead of (0.5,0.5) as in
     * ReconstructTrios).
     *
     * @return Double array
     */
    @Override
    protected double[][] getCentroids() {
        return new double[][]{{1.0, 0.0}, {0.5, 0.125}};
    }
}
