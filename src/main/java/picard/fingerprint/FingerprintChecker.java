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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.NotPrimaryAlignmentFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.*;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import picard.PicardException;

import java.io.File;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import static htsjdk.samtools.SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES;

/**
 * Major class that coordinates the activities involved in comparing genetic fingerprint
 * data whether the source is from a genotyping platform or derived from sequence data.
 *
 * @author Tim Fennell
 */
public class FingerprintChecker {
    public static final double DEFAULT_GENOTYPING_ERROR_RATE = 0.01;
    public static final int DEFAULT_MINIMUM_MAPPING_QUALITY = 10;
    public static final int DEFAULT_MINIMUM_BASE_QUALITY = 20;
    public static final int DEFAULT_MAXIMAL_PL_DIFFERENCE = 30;

    private final HaplotypeMap haplotypes;
    private int minimumBaseQuality      = DEFAULT_MINIMUM_BASE_QUALITY;
    private int minimumMappingQuality   = DEFAULT_MINIMUM_MAPPING_QUALITY;
    private double genotypingErrorRate  = DEFAULT_GENOTYPING_ERROR_RATE;
    private int maximalPLDifference     = DEFAULT_MAXIMAL_PL_DIFFERENCE;

    private boolean allowDuplicateReads = false;
    private double pLossofHet = 0;

    private final Log log = Log.getInstance(FingerprintChecker.class);

    /**
     * Creates a fingerprint checker that will work with the set of haplotypes stored in
     * the supplied file.
     */
    public FingerprintChecker(final File haplotypeData) {
        this.haplotypes = new HaplotypeMap(haplotypeData);
    }

    /** Creates a fingerprint checker that will work with the set of haplotyped provided. */
    public FingerprintChecker(final HaplotypeMap haplotypes) {
        this.haplotypes = haplotypes;
    }

    /** Sets the minimum base quality for bases used when computing a fingerprint from sequence data. */
    public void setMinimumBaseQuality(final int minimumBaseQuality) {
        this.minimumBaseQuality = minimumBaseQuality;
    }

    /** Sets the minimum mapping quality for reads used when computing fingerprints from sequence data. */
    public void setMinimumMappingQuality(final int minimumMappingQuality) {
        this.minimumMappingQuality = minimumMappingQuality;
    }

    /** Sets the assumed genotyping error rate used when accurate error rates are not available. */
    public void setGenotypingErrorRate(final double genotypingErrorRate) {
        this.genotypingErrorRate = genotypingErrorRate;
    }
    /** Sets the maximal difference in PL scores considered when reading PLs from a VCF. */
    public void setmaximalPLDifference(final int maximalPLDifference) {
        this.maximalPLDifference = maximalPLDifference;
    }

    public SAMFileHeader getHeader(){
        return haplotypes.getHeader();
    }
    /**
     * Sets whether duplicate reads should be allowed when calling genotypes from SAM files. This is
     * useful when comparing read groups within a SAM file and individual read groups show artifactually
     * high duplication (e.g. a single-ended read group mixed in with paired-end read groups).
     * @param allowDuplicateReads
     */
    public void setAllowDuplicateReads(final boolean allowDuplicateReads) {
        this.allowDuplicateReads = allowDuplicateReads;
    }

    //sets the value of the probability that a genotype underwent a Loss of Hetrozygosity (for Tumors)
    public void setpLossofHet(final double pLossofHet) {
        this.pLossofHet = pLossofHet;
    }

    /**
     * Loads genotypes from the supplied file into one or more Fingerprint objects and returns them in a
     * Map of Sample->Fingerprint.
     *
     * @param fingerprintFile - VCF file containing genotypes for one or more samples
     * @param specificSample - null to load genotypes for all samples contained in the file or the name
     *                         of an individual sample to load (and exclude all others).
     * @return a Map of Sample name to Fingerprint
     */
    public Map<String, Fingerprint> loadFingerprints(final File fingerprintFile, final String specificSample) {
        final VCFFileReader reader = new VCFFileReader(fingerprintFile, false);
        final CloseableIterator<VariantContext> iterator = reader.iterator();  

        SequenceUtil.assertSequenceDictionariesEqual(this.haplotypes.getHeader().getSequenceDictionary(),
                                                        VCFFileReader.getSequenceDictionary(fingerprintFile));

        final Map<String, Fingerprint> fingerprints = new HashMap<>();
        Set<String> samples = null;

        for( final VariantContext ctx : reader) {
            // Setup the sample names set if needed

            if (samples == null) {
                if (specificSample != null) {
                    samples = new HashSet<>();
                    samples.add(specificSample);
                } else {
                    samples = ctx.getSampleNames();
                    if (samples == null) {
                        log.warn("No samples found in file " + fingerprintFile + ". Skipping file.");
                        return Collections.emptyMap();
                    }
                }

                samples.forEach(s -> fingerprints.put(s, new Fingerprint(s, fingerprintFile, null)));
            }

            if (isUsableSnp(ctx)) {
                final HaplotypeBlock h = this.haplotypes.getHaplotype(ctx.getContig(), ctx.getStart());
                final Snp snp = this.haplotypes.getSnp(ctx.getContig(), ctx.getStart());
                if (h == null) continue;

                // Check the alleles from the file against the expected set of genotypes
                {
                    boolean allelesOk = true;
                    for (final Allele allele : ctx.getAlleles()) {
                        final byte[] bases = allele.getBases();
                        if (bases.length > 1 || (bases[0] != snp.getAllele1() && bases[0] != snp.getAllele2())) {
                            allelesOk = false;
                        }
                    }
                    if (!allelesOk) {
                        log.warn("Problem with genotype file '" + fingerprintFile.getName() + "': Alleles "
                                + ctx.getAlleles() + " do not match to alleles for SNP " + snp
                                + " with alleles " + snp.getAlleleString());
                        continue ;
                    }
                }

                for (final String sample : samples) {
                    final Fingerprint fp = fingerprints.get(sample);

                    //PLs are preferred over GTs
                    //TODO: this code is replicated in various places (ReconstructTriosFromVCF for example). Needs refactoring.
                    //TODO: add a way to force using GTs when both are available (why?)

                    // Get the genotype for the sample and check that it is useful
                    final Genotype genotype = ctx.getGenotype(sample);
                    if (genotype == null) {
                        throw new IllegalArgumentException("Cannot find sample " + sample + " in provided file: " + fingerprintFile);
                    }
                    if (genotype.hasPL()) {

                        final HaplotypeProbabilitiesFromGenotypeLikelihoods hFp = new HaplotypeProbabilitiesFromGenotypeLikelihoods(h);
                        //do not modify the PL array directly fragile!!!!!
                        final int[] pls = genotype.getPL();
                        final int[] newPLs = new int[pls.length];
                        for (int i = 0; i < pls.length; i++) {
                            newPLs[i] = Math.min(maximalPLDifference, pls[i]);
                        }
                        hFp.addToLogLikelihoods(snp, ctx.getAlleles(), GenotypeLikelihoods.fromPLs(newPLs).getAsVector());
                        fp.add(hFp);
                    } else {

                        if (genotype.isNoCall()) continue;

                        // TODO: when multiple genotypes are available for a Haplotype check that they
                        // TODO: agree. Not urgent since DownloadGenotypes already does this.
                        if (fp.containsKey(h)) continue;

                        final boolean hom = genotype.isHom();
                        final byte allele = StringUtil.toUpperCase(genotype.getAllele(0).getBases()[0]);

                        final double halfError = this.genotypingErrorRate / 2;
                        final double accuracy = 1 - this.genotypingErrorRate;
                        final double[] probs = new double[]{
                                (hom && allele == snp.getAllele1()) ? accuracy : halfError,
                                (!hom) ? accuracy : halfError,
                                (hom && allele == snp.getAllele2()) ? accuracy : halfError
                        };

                        fp.add(new HaplotypeProbabilitiesFromGenotype(snp, h, probs[0], probs[1], probs[2]));
                    }
                }
            }
        }

        return fingerprints;
    }

    /**
     * Quick method to check and see if the variant context represents a usable SNP variant. Unfortunately
     * ctx.isSnp doesn't always work if the genotype(s) are all monomorphic and the alternate allele isn't
     * listed.
     */
    public static boolean isUsableSnp(final VariantContext ctx) {
        if (ctx.isFiltered()) return false;
        if (ctx.isIndel()) return false;
        if (ctx.isMixed()) return false;

        // Also check that all alleles are length 1
        for (final Allele a : ctx.getAlleles()) {
            if (a.length() != 1) return false;
        }

        return true;
    }

    /**
     * Takes a set of fingerprints and returns an IntervalList containing all the loci that
     * can be productively examined in sequencing data to compare to one or more of the
     * fingerprints.
     */
    public IntervalList getLociToGenotype(final Collection<Fingerprint> fingerprints) {
        final IntervalList intervals = new IntervalList(this.haplotypes.getHeader());

        for (final Fingerprint fp : fingerprints) {
            for (final HaplotypeProbabilities genotype : fp.values()) {
                final HaplotypeBlock h = genotype.getHaplotype();
                for (final Snp snp : h.getSnps()) {
                    intervals.add(new Interval(snp.getChrom(), snp.getPos(), snp.getPos(), false, snp.getName()));
                }
            }
        }

        return intervals.uniqued();
    }

    public Map<FingerprintIdDetails, Fingerprint> fingerprintVcf(final File vcfFile) {
        Map<FingerprintIdDetails, Fingerprint> fpIdMap = new HashMap<>();

        final Map< String, Fingerprint> sampleFpMap = loadFingerprints(vcfFile, null);

        sampleFpMap.entrySet().stream().forEach(entry -> {
            FingerprintIdDetails fpId = new FingerprintIdDetails();
            fpId.sample = entry.getKey();
            fpId.file = vcfFile.getAbsolutePath();


            fpIdMap.put(fpId, entry.getValue());
        });
        return fpIdMap;
    }

    /**
     * Generates a Fingerprint per read group in the supplied SAM file using the loci provided in
     * the interval list.
     */
    public Map<FingerprintIdDetails, Fingerprint> fingerprintSamFile(final File samFile, final IntervalList loci) {
        final SamReader in = SamReaderFactory.makeDefault()
                                             .enable(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES)
                                             .open(samFile);

        SequenceUtil.assertSequenceDictionariesEqual(this.haplotypes.getHeader().getSequenceDictionary(),
                                                     in.getFileHeader().getSequenceDictionary());

        final SamLocusIterator iterator = new SamLocusIterator(in, loci, in.hasIndex());
        iterator.setEmitUncoveredLoci(true);
        iterator.setMappingQualityScoreCutoff(this.minimumMappingQuality);
        iterator.setQualityScoreCutoff(this.minimumBaseQuality);

        // In some cases it is useful to allow duplicate reads to be used - the most common is in single-end
        // sequence data where the duplicate marking may have been overly aggressive, and there is useful
        // non-redundant data in the reads marked as "duplicates'.
        if (this.allowDuplicateReads) {
            final List<SamRecordFilter> filters = new ArrayList<>(1);
            filters.add(new NotPrimaryAlignmentFilter());
            iterator.setSamFilters(filters);
        }

        final Map<FingerprintIdDetails, Fingerprint> fingerprintsByReadGroup = new HashMap<>();
        final Collection<FingerprintIdDetails> rgs = in.getFileHeader().getReadGroups().stream().map(rg->
        {FingerprintIdDetails id = new FingerprintIdDetails(rg.getPlatformUnit(), samFile.getAbsolutePath());
         id.library = rg.getLibrary();
         id.sample  = rg.getSample();
         return id;}).collect(Collectors.toSet());

        for (final FingerprintIdDetails rg : rgs) {
            final Fingerprint fingerprint = new Fingerprint(rg.sample,
                                                            samFile,
                                                            rg.platformUnit);
            fingerprintsByReadGroup.put(rg, fingerprint);

            for (final HaplotypeBlock h : this.haplotypes.getHaplotypes()) {
                fingerprint.add(new HaplotypeProbabilitiesFromSequence(h));
            }
        }

        // Set of read/template names from which we have already sampled a base and a qual. Since we assume
        // that all evidence for a haplotype is independent we can't sample two or more bases from a single
        // read or read-pair because they would not be independent!
        final Set<String> usedReadNames = new HashSet<>(10000);

        // Now go through the data at each locus and figure stuff out!
        for (final SamLocusIterator.LocusInfo info : iterator) {
            // TODO: Filter out the locus if the allele balance doesn't make sense for either a
            // TODO: 50/50 het or a hom with some errors; in HS data with deep coverage any base
            // TODO: with major strand bias could cause errors

            // Find the matching Snp and HaplotypeProbs
            final HaplotypeBlock haplotypeBlock = this.haplotypes.getHaplotype(info.getSequenceName(), info.getPosition());
            final Snp snp = this.haplotypes.getSnp(info.getSequenceName(), info.getPosition());

            for (final SamLocusIterator.RecordAndOffset rec : info.getRecordAndPositions()) {
                final SAMReadGroupRecord rg = rec.getRecord().getReadGroup();
                if (rg == null) {
                    final PicardException e = new PicardException("Found read with no readgroup: " + rec.getRecord().getReadName() + " in file: " + samFile);
                    log.error(e);
                    throw e;
                }
                final FingerprintIdDetails details = new FingerprintIdDetails(rg, samFile.getAbsolutePath());
                if (!fingerprintsByReadGroup.containsKey(details)) {
                    final PicardException e = new PicardException("Unknown read group: " + rg + " in file: " + samFile);
                    log.error(e);
                    throw e;
                } else {
                    final String readName = rec.getRecord().getReadName();
                    if (!usedReadNames.contains(readName)) {
                        final HaplotypeProbabilitiesFromSequence probs = (HaplotypeProbabilitiesFromSequence) fingerprintsByReadGroup.get(details).get(haplotypeBlock);
                        final byte base = StringUtil.toUpperCase(rec.getReadBase());
                        final byte qual = rec.getBaseQuality();

                        probs.addToProbs(snp, base, qual);
                        usedReadNames.add(readName);
                    }
                }
            }
        }

        return fingerprintsByReadGroup;
    }

    /**
     * Generates a per-sample Fingerprint for the contaminant in the supplied SAM file.
     * Data is aggregated by sample, not read-group.
     */
    public Map<String, Fingerprint> identifyContaminant(final File samFile, final double contamination, final int locusMaxReads) {
        final SamReader in = SamReaderFactory.makeDefault().enable(CACHE_FILE_BASED_INDEXES).open(samFile);
        SequenceUtil.assertSequenceDictionariesEqual(this.haplotypes.getHeader().getSequenceDictionary(),
                in.getFileHeader().getSequenceDictionary());

        final SamLocusIterator iterator = new SamLocusIterator(in, haplotypes.getIntervalList(), in.hasIndex());
        iterator.setEmitUncoveredLoci(true);
        iterator.setMappingQualityScoreCutoff(this.minimumMappingQuality);
        iterator.setQualityScoreCutoff(this.minimumBaseQuality);

        // In some cases it is useful to allow duplicate reads to be used - the most common is in single-end
        // sequence data where the duplicate marking may have been overly aggressive, and there is useful
        // non-redundant data in the reads marked as "duplicates'.
        if (this.allowDuplicateReads) {
            final List<SamRecordFilter> filters = new ArrayList<>(1);
            filters.add(new NotPrimaryAlignmentFilter());
            iterator.setSamFilters(filters);
        }

        final Map<String, Fingerprint> fingerprintsBySample = new HashMap<>();
        for (final SAMReadGroupRecord rg : in.getFileHeader().getReadGroups()) {
            if (!fingerprintsBySample.containsKey(rg.getSample())) {
                final Fingerprint fingerprint = new Fingerprint(rg.getSample(),
                        samFile,
                        rg.getSample());

                for (final HaplotypeBlock h : this.haplotypes.getHaplotypes()) {
                    fingerprint.add(new HaplotypeProbabilitiesFromContaminatorSequence(h, contamination));
                }
                fingerprintsBySample.put(rg.getSample(), fingerprint);
            }
        }

        // Set of read/template names from which we have already sampled a base and a qual. Since we assume
        // that all evidence for a haplotype is independent we can't sample two or more bases from a single
        // read or read-pair because they would not be independent!
        final Set<String> usedReadNames = new HashSet<>(10000);

        // Now go through the data at each locus and figure stuff out!
        for (final SamLocusIterator.LocusInfo info : iterator) {

            // Find the matching Snp and HaplotypeProbs
            final HaplotypeBlock haplotypeBlock = this.haplotypes.getHaplotype(info.getSequenceName(), info.getPosition());
            final Snp snp = this.haplotypes.getSnp(info.getSequenceName(), info.getPosition());

            // randomly select locusMaxReads elements from the list
            final List<SamLocusIterator.RecordAndOffset> recordAndOffsetList = randomSublist(info.getRecordAndPositions(), locusMaxReads);

            for (final SamLocusIterator.RecordAndOffset rec : recordAndOffsetList) {
                final SAMReadGroupRecord rg = rec.getRecord().getReadGroup();
                if (rg == null || !fingerprintsBySample.containsKey(rg.getSample())) {
                    final PicardException e = new PicardException("Unknown sample: " + (rg != null ? rg.getSample() : "(null readgroup)"));
                    log.error(e);
                    throw e;
                } else {
                    final String readName = rec.getRecord().getReadName();
                    if (!usedReadNames.contains(readName)) {
                        final HaplotypeProbabilitiesFromContaminatorSequence probs =
                                (HaplotypeProbabilitiesFromContaminatorSequence) fingerprintsBySample.get(rg.getSample()).get(haplotypeBlock);
                        final byte base = StringUtil.toUpperCase(rec.getReadBase());
                        final byte qual = rec.getBaseQuality();

                        probs.addToProbs(snp, base, qual);
                        usedReadNames.add(readName);
                    }
                }
            }
        }

        return fingerprintsBySample;
    }

    /**
     * A small utility function to choose n random elements (un-shuffled) from a list
     *
     * @param list A list of elements
     * @param n a number of elements requested from list
     * @return a list of n randomly chosen (but in the original order) elements from list.
     * If the list has less than n elements it is returned in its entirety.
     */
    protected static <T> List<T> randomSublist(final List<T> list, final int n) {
        int availableElements = list.size();
        if (availableElements <= n) return list;

        int stillNeeded = n;
        final Random rg = new Random();
        final List<T> shortList = new ArrayList<>(n);
        for (final T aList : list) {
            if (rg.nextDouble() < stillNeeded / (double) availableElements) {
                shortList.add(aList);
                stillNeeded--;
            }
            if (stillNeeded == 0 ) break; // fast out if do not need more elements
            availableElements--;
        }

        return shortList;
    }


    @Deprecated
    public Map<FingerprintIdDetails, Fingerprint> fingerprintSamFiles(final Collection<File> files, final int threads,
                                                                   final int waitTime, final TimeUnit waitTimeUnit) {
        return fingerprintFiles(files,  threads, waitTime, waitTimeUnit);
    }

        /**
         * Fingerprints one or more SAM/BAM/VCF files at all available loci within the haplotype map, using multiple threads
         * to speed up the processing.
         */
    public Map<FingerprintIdDetails, Fingerprint> fingerprintFiles(final Collection<File> files, final int threads,
                                                                    final int waitTime, final TimeUnit waitTimeUnit) {

        // Generate fingerprints from each file
        final AtomicInteger filesRead = new AtomicInteger(0);
        final ExecutorService executor = Executors.newFixedThreadPool(threads);
        final IntervalList intervals = this.haplotypes.getIntervalList();
        final Map<FingerprintIdDetails, Fingerprint> retval = new ConcurrentHashMap<>();

        for (final File f : files) {
            executor.submit(() -> {
                try {
                    if (CheckFingerprint.isBamOrSamFile(f)) {
                        retval.putAll(fingerprintSamFile(f, intervals));
                    } else {
                        retval.putAll(fingerprintVcf(f));
                    }

                    if (filesRead.incrementAndGet() % 100 == 0) {
                        log.info("Processed " + filesRead.get() + " out of " + files.size());
                    }
                } catch(Exception e){
                    log.warn("Exception thrown in thread:" + e.getMessage());
                    throw e;
                }
            });
        }
        executor.shutdown();
        try { executor.awaitTermination(waitTime, waitTimeUnit); }
        catch (final InterruptedException ie) { log.warn(ie, "Interrupted while waiting for executor to terminate."); }

        return retval;
    }

    /**
     * Takes a collection of fingerprints and, assuming that they are independent, merged the fingerprints
     * by samples and totals up the probabilities.
     */
    static public SortedMap<String, Fingerprint> mergeFingerprintsBySample(final Collection<Fingerprint> inputs) {
        final SortedMap<String, Fingerprint> sampleFps = new TreeMap<>();
        for (final Fingerprint fp : inputs) {
            Fingerprint sampleFp = sampleFps.get(fp.getSample());
            if (sampleFp == null) {
                sampleFp = new Fingerprint(fp.getSample(), null, fp.getSample());
                sampleFps.put(fp.getSample(), sampleFp);
            }

            sampleFp.merge(fp);
        }

        return sampleFps;
    }


    /**
     * Top level method to take a set of one or more SAM files and one or more Genotype files and compare
     * each read group in each SAM file to each set of fingerprint genotypes.
     *
     * @param samFiles the list of SAM files to fingerprint
     * @param genotypeFiles the list of genotype files from which to pull fingerprint genotypes
     * @param specificSample an optional single sample who's genotypes to load from the supplied files
     * @param ignoreReadGroups aggregate data into one fingerprint per file, instead of splitting by RG
     */
    public List<FingerprintResults> checkFingerprints(final List<File> samFiles,
                                                      final List<File> genotypeFiles,
                                                      final String specificSample,
                                                      final boolean ignoreReadGroups) {
        // Load the fingerprint genotypes
        final List<Fingerprint> expectedFingerprints = new LinkedList<>();
        for (final File f : genotypeFiles) {
            expectedFingerprints.addAll(loadFingerprints(f, specificSample).values());
        }

        if (expectedFingerprints.isEmpty()) {
            throw new IllegalStateException("Could not find any fingerprints in: " + genotypeFiles);
        }

        final List<FingerprintResults> resultsList = new ArrayList<>();
        final IntervalList intervals = getLociToGenotype(expectedFingerprints);

        // Fingerprint the SAM files and calculate the results
        for (final File f : samFiles) {
            final Map<FingerprintIdDetails, Fingerprint> fingerprintsByReadGroup = fingerprintSamFile(f, intervals);

            if (ignoreReadGroups) {
                final Fingerprint combinedFp = new Fingerprint(specificSample, f, null);
                fingerprintsByReadGroup.values().forEach(combinedFp::merge);

                final FingerprintResults results = new FingerprintResults(f, null, specificSample);
                for (final Fingerprint expectedFp : expectedFingerprints) {
                    final MatchResults result = calculateMatchResults(combinedFp, expectedFp, 0, pLossofHet);
                    results.addResults(result);
                }

                resultsList.add(results);

            } else {
                for (final FingerprintIdDetails rg : fingerprintsByReadGroup.keySet()) {
                    final FingerprintResults results = new FingerprintResults(f, rg.platformUnit, rg.sample);
                    for (final Fingerprint expectedFp : expectedFingerprints) {
                        final MatchResults result = calculateMatchResults(fingerprintsByReadGroup.get(rg), expectedFp, 0, pLossofHet);
                        results.addResults(result);
                    }

                    resultsList.add(results);
                }
            }
        }

        return resultsList;
    }

    /**
     * Top level method to take a set of one or more observed genotype (VCF) files and one or more expected genotype (VCF) files and compare
     * one or more sample in the observed genotype file with one or more in the expected file and generate results for each set.
     *
     * @param observedGenotypeFiles The list of genotype files containing observed calls, from which to pull fingerprint genotypes
     * @param expectedGenotypeFiles The list of genotype files containing expected calls, from which to pull fingerprint genotypes
     * @param observedSample an optional single sample whose genotypes to load from the observed genotype file (if null, use all)
     * @param expectedSample an optional single sample whose genotypes to load from the expected genotype file (if null, use all)
     */
    public List<FingerprintResults> checkFingerprints(final List<File> observedGenotypeFiles,
                                                      final List<File> expectedGenotypeFiles,
                                                      final String observedSample,
                                                      final String expectedSample) {

        // Load the expected fingerprint genotypes
        final List<Fingerprint> expectedFingerprints = new ArrayList<>();
        for (final File f : expectedGenotypeFiles) {
            expectedFingerprints.addAll(loadFingerprints(f, expectedSample).values());
        }

        if (expectedFingerprints.isEmpty()) {
            throw new IllegalStateException("Could not find any fingerprints in: " + expectedGenotypeFiles);
        }

        final List<FingerprintResults> resultsList = new ArrayList<>();

        for (final File f : observedGenotypeFiles) {
            final Map<String, Fingerprint> observedFingerprintsBySample = loadFingerprints(f, observedSample);
            if (observedFingerprintsBySample.isEmpty()) {
                throw new IllegalStateException("Found no fingerprints in observed genotypes file: " + observedGenotypeFiles);
            }

            for (final String sample : observedFingerprintsBySample.keySet()) {
                final FingerprintResults results = new FingerprintResults(f, null, sample);
                for (Fingerprint expectedFp : expectedFingerprints) {
                    final MatchResults result = calculateMatchResults(observedFingerprintsBySample.get(sample), expectedFp, 0, pLossofHet);
                    results.addResults(result);
                }
                resultsList.add(results);
            }
        }
        return resultsList;
    }


    /**
     * Compares two fingerprints and calculates a MatchResults object which contains detailed
     * information about the match (or mismatch) between fingerprints including the LOD score
     * for whether or not the two are likely from the same sample.
     *
     * If comparing sequencing data to genotype data then the sequencing data should be passed
     * as the observedFp and the genotype data as the expectedFp in order to get the best output.
     *
     * In the cases where the most likely genotypes from the two fingerprints do not match the
     * lExpectedSample is Max(actualpExpectedSample, minPExpected).
     */
    public static MatchResults calculateMatchResults(final Fingerprint observedFp, final Fingerprint expectedFp, final double minPExpected, final double pLoH) {
        final List<LocusResult> locusResults = new ArrayList<>();

        double llThisSample  = 0;
        double llOtherSample = 0;

        double lodExpectedSampleTumorNormal = 0;
        double lodExpectedSampleNormalTumor = 0;

        final double lminPExpected = Math.log10(minPExpected);

        for (final HaplotypeProbabilities probs2 : expectedFp.values()) {
            final HaplotypeBlock haplotypeBlock = probs2.getHaplotype();
            final HaplotypeProbabilities probs1 = observedFp.get(haplotypeBlock);
            if (probs1 == null) continue;

            final HaplotypeProbabilityOfNormalGivenTumor prob1AssumingDataFromTumor = new HaplotypeProbabilityOfNormalGivenTumor(probs1, pLoH);
            final HaplotypeProbabilityOfNormalGivenTumor prob2AssumingDataFromTumor = new HaplotypeProbabilityOfNormalGivenTumor(probs2, pLoH);

            // If one is from genotype data we'd like to report the output relative
            // to the genotyped SNP instead of against a random SNP from the haplotype
            final Snp snp = probs2.getRepresentativeSnp();
            final DiploidGenotype externalGenotype = probs2.getMostLikelyGenotype(snp);
            final LocusResult lr = new LocusResult(snp,
                    externalGenotype,
                    probs1.getMostLikelyGenotype(snp),
                    probs1.getObsAllele1(),
                    probs1.getObsAllele2(),
                    probs1.getLodMostProbableGenotype(),
                    // expected sample log-likelihood
                    probs1.shiftedLogEvidenceProbabilityGivenOtherEvidence(probs2),
                    // random sample log-likelihood
                    probs1.shiftedLogEvidenceProbability(),

            // probs1 is tumor probs2 is normal, correct sample lod
            prob1AssumingDataFromTumor.shiftedLogEvidenceProbabilityGivenOtherEvidence(probs2) -
                    prob1AssumingDataFromTumor.shiftedLogEvidenceProbability(),
                    // probs1 is normal probs2 is tumor, correct sample lod
                    probs1.shiftedLogEvidenceProbabilityGivenOtherEvidence(prob2AssumingDataFromTumor) -
                            probs1.shiftedLogEvidenceProbability());


            locusResults.add(lr);

            if (probs1.hasEvidence() && probs2.hasEvidence()) {
                final double lRandom = lr.lRandomSample();
                //TODO: what's the mathematics behind the lminPexpected?
                final double lExpected = Math.max(lminPExpected, lr.lExpectedSample());

                llThisSample  += lExpected;
                llOtherSample += lRandom;
                lodExpectedSampleTumorNormal += lr.getLodExpectedSampleTumorNormal();
                lodExpectedSampleNormalTumor += lr.getLodExpectedSampleNormalTumor();
            }
        }

        // TODO: prune the set of LocusResults for things that are too close together?
        return new MatchResults(expectedFp.getSource(), expectedFp.getSample(), llThisSample, llOtherSample, lodExpectedSampleTumorNormal, lodExpectedSampleNormalTumor, locusResults);
    }

    /**
     * Compares two fingerprints and calculates a MatchResults object which contains detailed
     * information about the match (or mismatch) between fingerprints including the LOD score
     * for whether or not the two are likely from the same sample.
     *
     * If comparing sequencing data to genotype data then the sequencing data should be passed
     * as the observedFp and the genotype data as the expectedFp in order to get the best output.
     */
    public static MatchResults calculateMatchResults(final Fingerprint observedFp, final Fingerprint expectedFp) {
        return calculateMatchResults(observedFp, expectedFp, 0, 0);
    }
}
