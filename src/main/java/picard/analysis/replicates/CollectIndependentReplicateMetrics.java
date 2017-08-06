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

package picard.analysis.replicates;


import htsjdk.samtools.DuplicateSet;
import htsjdk.samtools.DuplicateSetIterator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.FilteringSamIterator;
import htsjdk.samtools.filter.MappingQualityFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryOrSupplementaryFilter;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.ComparableTuple;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.FilteringVariantContextIterator;
import htsjdk.variant.variantcontext.filter.GenotypeQualityFilter;
import htsjdk.variant.variantcontext.filter.HeterozygosityFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.variantcontext.filter.SnpFilter;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Alpha;
import picard.filter.CountingPairedFilter;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static picard.cmdline.StandardOptionDefinitions.MINIMUM_MAPPING_QUALITY_SHORT_NAME;

/**
 * A CLP that, given a BAM and a VCF with genotypes of the same sample, estimates the rate of independent replication of reads within the bam.
 * That is, it estimates the fraction of the reads which look like duplicates (in the MarkDuplicates sense of the word) but are actually
 * independent observations of the data. In the presence of Unique Molecular Identifiers (UMIs), various metrics are collected regarding the
 * utility of the UMI's for the purpose of increasing coverage.
 * <p>
 * The estimation is based on duplicate-sets of size 2 and 3 and gives separate estimates from each. The assumption is that the duplication
 * rate (biological or otherwise) is independent of the duplicate-set size. A significant difference between the two rates may be an indication that
 * this assumption is incorrect.
 * <p>
 * The duplicate sets are found using the mate-cigar tag (MC) which is added by {@link picard.sam.MergeBamAlignment} , or {@link picard.sam.FixMateInformation}.
 * This program will not work without the MC tag.
 * <p>
 * Explanation of the calculation behind the estimation can be found in the {@link IndependentReplicateMetric} class.
 * <p>
 * The calculation Assumes a diploid organism (more accurately, assumes that only two alleles can appear at a HET site and that
 * these two alleles will appear at equal probabilities. It requires as input a VCF with genotypes for the sample in question.
 *
 * NOTE: This class is very much in alpha stage, and still under heavy development (feel free to join!)
 *
 *
 * @author Yossi Farjoun
 *
 */

@CommandLineProgramProperties(
        summary = "Estimates the rate of independent replication rate of reads within a bam. \n" +
                "That is, it estimates the fraction of the reads which would be marked as duplicates but " +
                "are actually biological replicates, independent observations of the data. ",
        oneLineSummary = "Estimates the rate of independent replication of reads within a bam.",
        programGroup = Alpha.class
)
public class CollectIndependentReplicateMetrics extends CommandLineProgram {
    private static final int DOUBLETON_SIZE = 2, TRIPLETON_SIZE = 3;

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input (indexed) BAM file.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Write metrics to this file")
    public File OUTPUT;

    @Argument(shortName = "MO", doc = "Write the confusion matrix (of UMIs) to this file", optional = true)
    public File MATRIX_OUTPUT;

    @Argument(shortName = "V", doc = "Input VCF file")
    public File VCF;

    @Argument(shortName = "GQ", doc = "minimal value for the GQ field in the VCF to use variant site.", optional = true)
    public Integer MINIMUM_GQ = 90;

    @Argument(shortName = MINIMUM_MAPPING_QUALITY_SHORT_NAME, doc = "minimal value for the mapping quality of the reads to be used in the estimation.", optional = true)
    public Integer MINIMUM_MQ = 40;

    @Argument(shortName = "BQ", doc = "minimal value for the base quality of a base to be used in the estimation.", optional = true)
    public Integer MINIMUM_BQ = 17;

    @Argument(shortName = StandardOptionDefinitions.SAMPLE_ALIAS_SHORT_NAME,
            doc = "Name of sample to look at in VCF. Can be omitted if VCF contains only one sample.", optional = true)
    public String SAMPLE = null;

    @Argument(doc = "Number of sets to examine before stopping.", optional = true)
    public Integer STOP_AFTER = 0;

    @Argument(doc = "Barcode SAM tag.", optional = true)
    public String BARCODE_TAG = "RX";

    @Argument(doc = "Barcode Quality SAM tag.", optional = true)
    public String BARCODE_BQ = "QX";

    @Argument(shortName = "MBQ", doc = "minimal value for the base quality of all the bases in a molecular barcode, for it to be used.", optional = true)
    public Integer MINIMUM_BARCODE_BQ = 30;

    private static final Log log = Log.getInstance(CollectIndependentReplicateMetrics.class);

    @Override
    protected int doWork() {

        IOUtil.assertFileIsReadable(VCF);
        IOUtil.assertFileIsReadable(INPUT);

        IOUtil.assertFileIsWritable(OUTPUT);
        if (MATRIX_OUTPUT != null) IOUtil.assertFileIsWritable(MATRIX_OUTPUT);

        final VCFFileReader vcf = new VCFFileReader(VCF, false);

        final VCFHeader vcfFileHeader = vcf.getFileHeader();
        final List<String> samples = vcfFileHeader.getSampleNamesInOrder();

        if (SAMPLE == null) {
            if (samples.size() != 1) {
                throw new IllegalArgumentException("When sample is null, VCF must have exactly 1 sample. found " + samples.size());
            } else {
                SAMPLE = samples.get(0);
                log.info("No SAMPLE given, using sample from VCF: ", SAMPLE);
            }
        } else if (!samples.contains(SAMPLE)) {
            throw new IllegalArgumentException("When sample is not null, VCF must contain supplied sample. Cannot find sample " + SAMPLE + " in vcf.");
        }

        final Histogram<ComparableTuple<String, String>> umiConfusionMatrix = new Histogram<>("ConfusionUMI", "Count");
        final Histogram<ComparableTuple<String, String>> umiConfusionMatrixEditDistance = new Histogram<>("ConfusionUMI", "EditDistance");

        final IndependentReplicateMetric metric = new IndependentReplicateMetric();

        final Histogram<Byte> umiEditDistanceInDiffBiDups = new Histogram<>("editDistance", "diffAllelesCount");
        final Histogram<Byte> umiEditDistanceInSameBiDups = new Histogram<>("editDistance", "sameAllelesCount");
        final Histogram<Byte> alleleBalanceCount = new Histogram<>("alleleBalance", "alleleBalanceCount");

        // get the intervals that correspond to het sites in the VCF
        final SortedMap<QueryInterval, List<Allele>> intervalAlleleMap = getQueryIntervalsMap(VCF);
        final Iterator<QueryInterval> queryIntervalIterator = intervalAlleleMap.keySet().iterator();

        log.info("Found " + intervalAlleleMap.size() + " heterozygous sites in VCF.");

        // get an iterator to reads that overlap the heterozygous sites
        final SamReader in = SamReaderFactory.makeDefault().open(INPUT);

        log.info("Querying BAM for sites.");

        final SAMRecordIterator samRecordIterator = in.query(intervalAlleleMap.keySet().toArray(new QueryInterval[intervalAlleleMap.size()]), false);
        final List<SamRecordFilter> samFilters = CollectionUtil.makeList(
                new AlignedFilter(true),
                new CountingPairedFilter(),
                new SecondaryOrSupplementaryFilter(),
                new MappingQualityFilter(MINIMUM_MQ)
        );

        final FilteringSamIterator filteredSamRecordIterator = new FilteringSamIterator(samRecordIterator, new AggregateFilter(samFilters));
        log.info("Queried BAM, getting duplicate sets.");

        // get duplicate iterator from iterator above
        final DuplicateSetIterator duplicateSets = new DuplicateSetIterator(filteredSamRecordIterator, in.getFileHeader());

        QueryInterval queryInterval = null;

        log.info("Starting iteration on reads");
        final ProgressLogger progress = new ProgressLogger(log, 10000000, "examined", "duplicate sets");

        IndependentReplicateMetric locusData = new IndependentReplicateMetric();
        boolean useLocus = true;
        boolean newLocus = false;
        int thirdAlleleInfos = 0;
        Allele badAllele = null;
        String offendingReadName = null;

        set:
        while (duplicateSets.hasNext()) {

            final DuplicateSet set = duplicateSets.next();
            final SAMRecord setRep = set.getRepresentative();
            final QueryInterval setRepsInterval = queryIntervalFromSamRecord(setRep);

            progress.record(setRep);
            // if the current duplicate set no longer overlaps the query interval then null it (and handle it below)
            // also move to the next variant if the previous variant is bad.
            if (!useLocus || queryInterval != null && isCleanlyBefore(queryInterval, setRepsInterval)) {
                if (!useLocus) {
                    metric.nThreeAllelesSites++;
                    if (++thirdAlleleInfos < 100) {
                        log.debug("Skipping a locus due to third allele: " + badAllele + " but expected " +
                                intervalAlleleMap.get(queryInterval) + " queryInterval " + queryInterval +
                                " offending read name is : " + offendingReadName);
                    }
                }
                queryInterval = null;
            }

            // Iterate until we find the query interval that contains the current duplicate set.

            // Simply polling for the "next" query will not do since the next one might not be covered by any reads, or it may have been
            // covered by past reads (if close enough to previous query interval)
            while (queryIntervalIterator.hasNext() &&
                    (queryInterval == null || isCleanlyBefore(queryInterval, setRepsInterval))) {
                // if we haven't seen either the reference or the alternate in the locus (subject to our stringent filters) do not use locus.
                if (locusData.nReferenceReads == 0 || locusData.nAlternateReads == 0) {
                    useLocus = false;
                    log.debug("will not use this locus due to lack of evidence of het site.");
                }
                // Query interval didn't get killed by 3rd alleles and so we combine the results with the tally
                if (useLocus && newLocus) {
                    metric.merge(locusData);
                    log.debug("merging metric. total nSites so far: " + metric.nSites);
                    //calculate allele balance with faux counts
                    final byte alleleBalance = (byte) Math.round(100D * (locusData.nAlternateReads + 0.5) / (locusData.nAlternateReads + locusData.nReferenceReads + 1));
                    alleleBalanceCount.increment(alleleBalance);
                    // we have merged now, no need to merge the old locus data or update the nSites until out of this while.
                    newLocus = false;
                }
                queryInterval = queryIntervalIterator.next();
                locusData = new IndependentReplicateMetric();
                locusData.nSites = 1;
                useLocus = true;
            }
            // we have a new locus, next time we should perhaps merge
            newLocus = true;

            // shouldn't happen, but being safe.
            if (queryInterval == null) break;

            final int setSize = set.size();

            locusData.nTotalReads += setSize;

            if (setSize > 1) locusData.nDuplicateSets++;
            if (setSize == DOUBLETON_SIZE) {
                locusData.nExactlyDouble++;
            } else if (setSize == TRIPLETON_SIZE) {
                locusData.nExactlyTriple++;
            } else if (setSize > TRIPLETON_SIZE) {  // singletons are only counted in nTotalReads
                locusData.nReadsInBigSets += setSize;
            }

            log.debug("set size is: " + setSize);
            final List<Allele> allelesInVc = intervalAlleleMap.get(queryInterval);

            log.debug("alleles in VC: " + allelesInVc);

            int nRef = 0, nAlt = 0, nOther = 0;
            for (final SAMRecord read : set.getRecords()) {

                // getReadPositionAtReferencePosition gives 1-based offset
                final int offset = read.getReadPositionAtReferencePosition(queryInterval.start) - 1;

                if (offset == -1) {
                    // a labeled continue watch-out!
                    // This could be a deletion OR a clipped end. Get a new set.
                    log.debug("got offset -1, getting new set");
                    continue set;
                }
                // a labeled continue watch-out!
                // need to move to the next set since this set has a low quality base-quality.

                if (read.getBaseQualities()[offset] <= MINIMUM_BQ) {
                    log.debug("got low read quality, getting new set");
                    continue set;
                }

                final Allele allele = Allele.create(read.getReadBases()[offset]);

                if (allelesInVc.get(0).basesMatch(allele)) {
                    nRef++;
                } else if (allelesInVc.get(1).basesMatch(allele)) {
                    nAlt++;
                } else {
                    nOther++;
                    // if other alleles were found, toss out the whole locus! (but read the reads first)
                    useLocus = false;
                    badAllele = allele;
                    offendingReadName = read.getReadName();
                }
            }
            locusData.nAlternateReads += nAlt;
            locusData.nReferenceReads += nRef;

            if ( setSize == 1 || setSize > TRIPLETON_SIZE) continue;
            // From here on there should only be 2 or 3 reads in the set

            final SetClassification classification = classifySet(nRef, nAlt, nOther);

            log.debug("Classification of set is: " + classification);
            if (setSize == DOUBLETON_SIZE) {

                final boolean useBarcodes = !set.getRecords().stream()
                        .map(read -> read.getStringAttribute(BARCODE_BQ))
                        .map(string -> string == null ? "" : string).map(string ->
                        {
                            final byte[] bytes = SAMUtils.fastqToPhred(string);
                            return IntStream.range(0, bytes.length).map(i -> bytes[i]).anyMatch(q -> q < MINIMUM_BARCODE_BQ);
                        }).anyMatch(a -> a);

                log.debug("using barcodes?" + useBarcodes);

                if(useBarcodes) locusData.nGoodBarcodes++; else locusData.nBadBarcodes++;

                final List<String> barcodes = set.getRecords().stream()
                        .map(read -> read.getStringAttribute(BARCODE_TAG))
                        .map(string -> string == null ? "" : string).collect(Collectors.toList());

                log.debug("found UMIs:" + barcodes);
                final boolean hasMultipleOrientations = set.getRecords().stream()
                        .map(SAMRecord::getFirstOfPairFlag) //must be paired, due to filter on sam Iterator
                        .distinct().count() != 1;
                log.debug("reads have multiple orientation?" + hasMultipleOrientations);

                final byte editDistance = calculateEditDistance(barcodes.get(0), barcodes.get(1));

                log.debug("Edit distance between umi: " + editDistance);

                if (useBarcodes && editDistance != 0) {
                    if (hasMultipleOrientations) locusData.nMismatchingUMIsInContraOrientedBiDups++;
                    else locusData.nMismatchingUMIsInCoOrientedBiDups++;
                }

                // sanity check the number of distinct tags
                if (classification == SetClassification.DIFFERENT_ALLELES) {
                    locusData.nDifferentAllelesBiDups++;
                    if (useBarcodes) {
                        umiEditDistanceInDiffBiDups.increment(editDistance);

                        if(editDistance == 0) locusData.nMatchingUMIsInDiffBiDups++; else locusData.nMismatchingUMIsInDiffBiDups++;
                    }

                    // we're going to toss out this locus.
                } else if (classification == SetClassification.MISMATCHING_ALLELE) {
                    locusData.nMismatchingAllelesBiDups++;
                } else { // the classification is either ALTERNATE_ALLELE or REFERENCE_ALLELE if we've reached here
                    if (classification == SetClassification.ALTERNATE_ALLELE) locusData.nAlternateAllelesBiDups++;
                    else locusData.nReferenceAllelesBiDups++;

                    if (useBarcodes) {

                        umiEditDistanceInSameBiDups.increment(editDistance);
                        final ComparableTuple<String, String> key = new ComparableTuple<>(barcodes.get(0), barcodes.get(1));
                        umiConfusionMatrix.increment(key);
                        if (!umiConfusionMatrixEditDistance.containsKey(key)) umiConfusionMatrixEditDistance.increment(key, editDistance);

                        if (editDistance == 0) locusData.nMatchingUMIsInSameBiDups++; else  locusData.nMismatchingUMIsInSameBiDups++;
                    }
                }
            }
            if (setSize == TRIPLETON_SIZE) {
                switch (classification) {
                    case MISMATCHING_ALLELE:
                        locusData.nMismatchingAllelesTriDups++;
                        break;
                    case DIFFERENT_ALLELES:
                        locusData.nDifferentAllelesTriDups++;
                        break;
                    case ALTERNATE_ALLELE:
                        locusData.nAlternateAllelesTriDups++;
                        break;
                    case REFERENCE_ALLELE:
                        locusData.nReferenceAllelesTriDups++;
                        break;
                    default:
                        throw new IllegalStateException("Un possible!");
                }
            }
            if (STOP_AFTER > 0 && progress.getCount() > STOP_AFTER) break;
        }
        if (useLocus && newLocus) {
            metric.merge(locusData);
            log.debug("Merged final metric. nSites:" + metric.nSites);
        } else {
            metric.nThreeAllelesSites++;
            log.debug("didn't merge last metric, due to 3rd allele: nThreeAllelesSites =" + metric.nThreeAllelesSites);
        }

        log.info("Iteration done. Emitting metrics.");

        // Emit metrics
        final MetricsFile<IndependentReplicateMetric, Byte> metricsFile = getMetricsFile();

        metric.calculateDerivedFields();
        metricsFile.addMetric(metric);
        metricsFile.addHistogram(alleleBalanceCount);
        metricsFile.addHistogram(umiEditDistanceInDiffBiDups);
        metricsFile.addHistogram(umiEditDistanceInSameBiDups);

        metricsFile.write(OUTPUT);

        final MetricsFile<?, ComparableTuple<String, String>> confusionMetrics = getMetricsFile();

        if (MATRIX_OUTPUT != null) {
            confusionMetrics.addHistogram(umiConfusionMatrix);
            confusionMetrics.addHistogram(umiConfusionMatrixEditDistance);
            confusionMetrics.write(MATRIX_OUTPUT);
        }

        return 0;
    }

    private enum SetClassification {
        MISMATCHING_ALLELE,
        DIFFERENT_ALLELES,
        REFERENCE_ALLELE,
        ALTERNATE_ALLELE
    }

    /**
     * a small utility to inform if one interval is cleanly before another, meaning that they do not overlap and
     * the first is prior (in genomic order) to the second
     *
     * @param lhs the "first" {@link QueryInterval}
     * @param rhs the "second" {@link QueryInterval}
     * @return true if the to intervals do not intersect _and_ the first is prior to the second in genomic order
     */
    private static boolean isCleanlyBefore(final QueryInterval lhs, final QueryInterval rhs) {
        return !lhs.overlaps(rhs) && lhs.compareTo(rhs) < 0;
    }

    private static SetClassification classifySet(final int nRef, final int nAlt, final int nOther) {
        // if we found any "other" alleles, this is a mismatching set
        if (nOther != 0) return SetClassification.MISMATCHING_ALLELE;

        // if we found both ref and alt alleles, this is a heterogeneous set
        if (nAlt > 0 && nRef > 0) return SetClassification.DIFFERENT_ALLELES;

        // if we found no reference alleles, this is an "alternate" set
        if (nRef == 0) return SetClassification.ALTERNATE_ALLELE;

        // if we found no alternate alleles, this is a "reference" set.
        if (nAlt == 0) return SetClassification.REFERENCE_ALLELE;

        throw new IllegalAccessError("shouldn't be here!");
    }

    private static QueryInterval queryIntervalFromSamRecord(final SAMRecord samRecord) {
        return new QueryInterval(samRecord.getReferenceIndex(), samRecord.getStart(), samRecord.getEnd());
    }

    /** Gives the edit distance between this barcode and another of the same length. */
    private static byte calculateEditDistance(final String lhs, final String rhs) {
        assert(lhs.length()==rhs.length());
        byte tmp = 0;
        for (int i = 0; i < rhs.length(); ++i) {
            if (rhs.charAt(i) != lhs.charAt(i)) ++tmp;
        }
        return tmp;
    }

    private SortedMap<QueryInterval, List<Allele>> getQueryIntervalsMap(final File vcf) {

        final Map<String, Integer> contigIndexMap = new HashMap<>();
        final VCFFileReader vcfReader = new VCFFileReader(vcf, false);

        // We want to look at unfiltered SNP sites for which the sample is genotyped as a het
        // with high quality.
        final CompoundFilter compoundFilter = new CompoundFilter(true);
        compoundFilter.add(new SnpFilter());
        compoundFilter.add(new PassingVariantFilter());
        compoundFilter.add(new GenotypeQualityFilter(MINIMUM_GQ, SAMPLE));
        compoundFilter.add(new HeterozygosityFilter(true, SAMPLE));

        final Iterator<VariantContext> hetIterator = new FilteringVariantContextIterator(vcfReader.iterator(), compoundFilter);

        for (final VCFContigHeaderLine vcfContig : vcfReader.getFileHeader().getContigLines()) {
            contigIndexMap.put(vcfContig.getID(), vcfContig.getContigIndex());
        }

        // return a TreeMap since the keys are comparable, and this will use their order in the iteration
        final SortedMap<QueryInterval, List<Allele>> map = new TreeMap<>();

        while (hetIterator.hasNext()) {
            final VariantContext vc = hetIterator.next();
            map.put(new QueryInterval(contigIndexMap.get(vc.getContig()), vc.getStart(), vc.getEnd()), vc.getGenotype(SAMPLE).getAlleles());
        }

        vcfReader.close();

        return map;
    }
}
