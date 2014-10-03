/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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
package picard.vcf;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.vcf.VCFFileReader;
import picard.PicardException;
import picard.cmdline.CommandLineParser;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.Usage;

import java.io.File;
import java.util.Iterator;
import java.util.List;

/**
 * Calculates the concordance between genotype data for two samples in two different VCFs - one representing the truth (or reference)
 * the other being the call.  The concordance is broken into separate results sections for SNPs and indels.  Summary and detailed statistics
 * are reported
 *
 * @author Tim Fennell
 * @author George Grant
 */
public class GenotypeConcordance extends CommandLineProgram {
    @Usage
    public final String USAGE =
            CommandLineParser.getStandardUsagePreamble(getClass()) +
                    "Calculates the concordance between genotype data for two samples in two different VCFs - one representing the truth (or reference) " +
                    "the other being the call.  The concordance is broken into separate results sections for SNPs and indels.  Summary and detailed statistics are reported\n\n" +
                    "Note that for any pair of variants to compare, a GATK.subcontext is generated for just the specified sample's genotypes and only variants that are " +
                    "NOT MNP, Symbolic, nor Mixed in either variant are used.  Variants (loci) that are invariant in both samples are not used.";

    @Option(shortName = "TV", doc="The VCF containing the truth sample")
    public File TRUTH_VCF;

    @Option(shortName = "CV", doc="The VCF containing the call sample")
    public File CALL_VCF;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Basename for the two metrics files that are to be written." +
            " Resulting files will be <OUTPUT>" + SUMMARY_METRICS_FILE_EXTENSION + "  and <OUTPUT>" + DETAILED_METRICS_FILE_EXTENSION + ".")
    public File OUTPUT;

    @Option(shortName="DM", doc="Name for the debug (variant-categories) metrics file.", optional = true)
    public File DEBUG_VARIANT_METRICS_FILE;

    @Option(shortName = "TS", doc="The name of the truth sample within the truth VCF")
    public String TRUTH_SAMPLE;

    @Option(shortName = "CS", doc="The name of the call sample within the call VCF")
    public String CALL_SAMPLE;

    @Option(doc="One or more interval list files that will be used to limit the genotype concordance.")
    public List<File> INTERVALS;

    @Option(doc="If true, multiple interval lists will be intersected. If false multiple lists will be unioned.")
    public boolean INTERSECT_INTERVALS = true;

    @Option(doc="Genotypes below this genotype quality will have genotypes classified as LowGq.")
    public int MIN_GQ = 0;

    @Option(doc="Genotypes below this depth will have genotypes classified as LowDp.")
    public int MIN_DP = 0;

    @Option(doc="If true, output all rows in detailed statistics even when count == 0.  When false only output rows with non-zero counts.")
    public boolean OUTPUT_ALL_ROWS = false;

    @Option(doc="If accessing less than this percentage of the file, use the VCF index.", optional=true)
    public int USE_VCF_INDEX_PCT = 25;

    private final Log log = Log.getInstance(GenotypeConcordance.class);
    private final ProgressLogger progress = new ProgressLogger(log, 10000, "checked", "variants");

    public static final String SUMMARY_METRICS_FILE_EXTENSION = ".summary_metrics.txt";
    public static final String DETAILED_METRICS_FILE_EXTENSION = ".detailed_metrics.txt";

    protected ConcordanceResults snpCounter;
    public ConcordanceResults getSnpCounter() { return snpCounter; }

    protected ConcordanceResults indelCounter;
    public ConcordanceResults getIndelCounter() { return indelCounter; }

    // TODO: add optimization if the samples are in the same file
    // TODO: add option for auto-detect pairs based on same sample name
    // TODO: allow multiple sample-pairs in one pass

    public static void main(final String[] args) {
        new GenotypeConcordance().instanceMainWithExit(args);
    }

    protected String[] customCommandLineValidation() {
        if (USE_VCF_INDEX_PCT < 0 || USE_VCF_INDEX_PCT > 100) {
            return new String[]{"USE_VCF_INDEX_PCT " + USE_VCF_INDEX_PCT + " out of range (0-100)"};
        }
        return super.customCommandLineValidation();
    }

    @Override protected int doWork() {
        IOUtil.assertFileIsReadable(TRUTH_VCF);
        IOUtil.assertFileIsReadable(CALL_VCF);
        final File summaryMetricsFile = new File(OUTPUT + SUMMARY_METRICS_FILE_EXTENSION);
        final File detailedMetricsFile = new File(OUTPUT + DETAILED_METRICS_FILE_EXTENSION);
        IOUtil.assertFileIsWritable(summaryMetricsFile);
        IOUtil.assertFileIsWritable(detailedMetricsFile);

        if (DEBUG_VARIANT_METRICS_FILE != null) IOUtil.assertFileIsWritable(DEBUG_VARIANT_METRICS_FILE);

        final boolean usingIntervals = this.INTERVALS != null && this.INTERVALS.size() > 0;
        IntervalList intervals = null;
        SAMSequenceDictionary intervalsSamSequenceDictionary = null;
        long genomeBaseCount = 0;
        boolean useIndex = false;
        if (usingIntervals) {
            log.info("Loading up region lists.");
            for (final File f : INTERVALS) {
                IOUtil.assertFileIsReadable(f);
                final IntervalList tmpIntervalList = IntervalList.fromFile(f);
                if (genomeBaseCount == 0) {         // Don't count the reference length more than once.
                    intervalsSamSequenceDictionary = tmpIntervalList.getHeader().getSequenceDictionary();
                    genomeBaseCount = intervalsSamSequenceDictionary.getReferenceLength();
                }

                if (intervals == null)        intervals = tmpIntervalList;
                else if (INTERSECT_INTERVALS) intervals = IntervalList.intersection(intervals, tmpIntervalList);
                else intervals =              IntervalList.union(intervals, tmpIntervalList);
            }
            if (intervals != null) {
                intervals = intervals.uniqued();
                // Optimization.  Only use the VCF index if we are using < USE_VCF_INDEX_PCT% of the file
                useIndex = USE_VCF_INDEX_PCT > (100 * intervals.getBaseCount() / genomeBaseCount);
            }
        }

        final VCFFileReader truthReader = new VCFFileReader(TRUTH_VCF, useIndex);
        final VCFFileReader callReader = new VCFFileReader(CALL_VCF, useIndex);

        // Check that the samples actually exist in the files!
        if (!truthReader.getFileHeader().getGenotypeSamples().contains(TRUTH_SAMPLE)) {
            throw new PicardException("File " + TRUTH_VCF.getAbsolutePath() + " does not contain genotypes for sample " + TRUTH_SAMPLE);
        }
        if (!callReader.getFileHeader().getGenotypeSamples().contains(CALL_SAMPLE)) {
            throw new PicardException("File " + CALL_VCF.getAbsolutePath() + " does not contain genotypes for sample " + CALL_SAMPLE);
        }

        // Verify that both VCFs have the same Sequence Dictionary
        SequenceUtil.assertSequenceDictionariesEqual(truthReader.getFileHeader().getSequenceDictionary(), callReader.getFileHeader().getSequenceDictionary());

        if (usingIntervals) {
            // If using intervals, verify that the sequence dictionaries agree with those of the VCFs
            SequenceUtil.assertSequenceDictionariesEqual(intervalsSamSequenceDictionary, truthReader.getFileHeader().getSequenceDictionary());
        } else {
            // If not using intervals, we need to determine the reference sequence length from the VCF's sequence dictionary
            genomeBaseCount = truthReader.getFileHeader().getSequenceDictionary().getReferenceLength();
        }

        // Build the pair of iterators over the regions of interest
        final Iterator<VariantContext> truthIterator, callIterator;
        if (usingIntervals) {
            truthIterator = new ByIntervalListVariantContextIterator(truthReader, intervals);
            callIterator = new ByIntervalListVariantContextIterator(callReader, intervals);
        }
        else {
            truthIterator = truthReader.iterator();
            callIterator = callReader.iterator();
        }

        // Now do the iteration and count things up
        final PairedVariantSubContextIterator iterator = new PairedVariantSubContextIterator(truthIterator, TRUTH_SAMPLE, callIterator, CALL_SAMPLE, truthReader.getFileHeader().getSequenceDictionary());
        snpCounter   = new ConcordanceResults(TRUTH_SAMPLE, CALL_SAMPLE);
        indelCounter = new ConcordanceResults(TRUTH_SAMPLE, CALL_SAMPLE);

        final VariantContextTypeResults variantContextTypeResults = new VariantContextTypeResults(TRUTH_SAMPLE, CALL_SAMPLE);

        log.info("Starting iteration over variants.");
        while (iterator.hasNext()) {
            final VcTuple tuple = iterator.next();
            final VariantCallState truthState = determineState(tuple.truthVariantContext, TRUTH_SAMPLE, MIN_GQ, MIN_DP);
            final VariantCallState callState = determineState(tuple.callVariantContext, CALL_SAMPLE, MIN_GQ, MIN_DP);

            final VariantContext.Type truthVariantContextType = (tuple.truthVariantContext == null) ? VariantContext.Type.NO_VARIATION : tuple.truthVariantContext.getType();
            final VariantContext.Type callVariantContextType  = (tuple.callVariantContext  == null) ? VariantContext.Type.NO_VARIATION  : tuple.callVariantContext.getType();

            variantContextTypeResults.add(truthVariantContextType, callVariantContextType);

            final boolean allelesAgree = doAllelesAgree(tuple.truthVariantContext, truthState, tuple.callVariantContext, callState);
            if ((truthVariantContextType == VariantContext.Type.SNP) || ((truthVariantContextType == VariantContext.Type.NO_VARIATION) && (callVariantContextType == VariantContext.Type.SNP))) {
                snpCounter.add(truthState, callState, allelesAgree);
            }
            else if ((truthVariantContextType == VariantContext.Type.INDEL) || ((truthVariantContextType == VariantContext.Type.NO_VARIATION) && (callVariantContextType == VariantContext.Type.INDEL))) {
                indelCounter.add(truthState, callState, allelesAgree);
            }

            final VariantContext variantContextForLogging = tuple.truthVariantContext != null ? tuple.truthVariantContext : tuple.callVariantContext;
            progress.record(variantContextForLogging.getChr(), variantContextForLogging.getStart());
        }

        // Calculate and store the summary-level metrics
        final MetricsFile<GenotypeConcordanceSummaryMetrics,?> genotypeConcordanceSummaryMetricsFile = getMetricsFile();
        GenotypeConcordanceSummaryMetrics m = new GenotypeConcordanceSummaryMetrics("SNP", snpCounter);
        m.FALSE_POS_PER_MB = 1e6 * snpCounter.numFalsePositives() / genomeBaseCount;
        genotypeConcordanceSummaryMetricsFile.addMetric(m);
        m = new GenotypeConcordanceSummaryMetrics("Indel", indelCounter);
        m.FALSE_POS_PER_MB = 1e6 * indelCounter.numFalsePositives() / genomeBaseCount;
        genotypeConcordanceSummaryMetricsFile.addMetric(m);
        genotypeConcordanceSummaryMetricsFile.write(summaryMetricsFile);

        // Calculate and store the detailed metrics for both SNP and indels
        final MetricsFile<GenotypeConcordanceDetailMetrics,?> genotypeConcordanceDetailMetrics = getMetricsFile();
        outputDetailMetricsFile(genotypeConcordanceDetailMetrics, snpCounter, "SNP");
        outputDetailMetricsFile(genotypeConcordanceDetailMetrics, indelCounter, "Indel");
        genotypeConcordanceDetailMetrics.write(detailedMetricsFile);

        if (DEBUG_VARIANT_METRICS_FILE != null) {
            final MetricsFile<GenotypeConcordanceVariantTypeDetailsMetrics,?> genotypeConcordanceVariantTypeDetailsMetricsFile = getMetricsFile();
            outputVariantContextTypeDetailsTable(genotypeConcordanceVariantTypeDetailsMetricsFile, variantContextTypeResults);
            genotypeConcordanceVariantTypeDetailsMetricsFile.write(DEBUG_VARIANT_METRICS_FILE);
        }

        return 0;
    }

    /**
     * Outputs the detailed statistics tables for SNP and Indel match categories.
     **/
    private void outputDetailMetricsFile(final MetricsFile<GenotypeConcordanceDetailMetrics,?> genotypeConcordanceDetailMetricsFile, final ConcordanceResults counter, final String eventType) {
        for (final VariantCallState truthSate : VariantCallState.values()) {
            for (final VariantCallState callState : VariantCallState.values()) {
                for (final boolean altAllelesAgree : new boolean[] { true, false}) {
                    final long count = counter.getCount(truthSate, callState, altAllelesAgree);
                    if (count > 0 || OUTPUT_ALL_ROWS) {
                        final GenotypeConcordanceDetailMetrics detailMetrics = new GenotypeConcordanceDetailMetrics();
                        detailMetrics.EVENT_TYPE = eventType;
                        detailMetrics.TRUTH_SAMPLE = counter.getTruthSample();
                        detailMetrics.CALL_SAMPLE = counter.getCallSample();
                        detailMetrics.TRUTH_STATE = truthSate.name();
                        detailMetrics.CALL_STATE = callState.name();
                        detailMetrics.ALT_ALLELES_AGREE = altAllelesAgree;
                        detailMetrics.COUNT = count;
                        genotypeConcordanceDetailMetricsFile.addMetric(detailMetrics);
                    }
                }
            }
        }
    }

    /**
     * Outputs the detailed statistics tables for SNP and Indel match categories.
     **/
    private void outputVariantContextTypeDetailsTable(final MetricsFile<GenotypeConcordanceVariantTypeDetailsMetrics,?> genotypeConcordanceVariantTypeDetailsMetricsFile, final VariantContextTypeResults counter) {
        for (final VariantContext.Type truthType : VariantContext.Type.values()) {
            for (final VariantContext.Type callType : VariantContext.Type.values()) {
                final long count = counter.getCount(truthType, callType);
                if (count > 0 || OUTPUT_ALL_ROWS) {
                    final GenotypeConcordanceVariantTypeDetailsMetrics typeMetrics = new GenotypeConcordanceVariantTypeDetailsMetrics();
                    typeMetrics.TRUTH_SAMPLE = snpCounter.getTruthSample();
                    typeMetrics.CALL_SAMPLE = snpCounter.getCallSample();
                    typeMetrics.TRUTH_TYPE = truthType.name();
                    typeMetrics.CALL_TYPE = callType.name();
                    typeMetrics.COUNT = count;
                    genotypeConcordanceVariantTypeDetailsMetricsFile.addMetric(typeMetrics);
                }
            }
        }
    }


    /** Determines the classification for a single sample at a single locus. */
    final VariantCallState determineState(final VariantContext ctx, final String sample, final int minGq, final int minDp) {
        // Site level checks
        if (ctx == null) return VariantCallState.NoVariant;
        else if (ctx.isFiltered()) return VariantCallState.FilteredVariant;

        // Genotype level checks
        final Genotype gt = ctx.getGenotype(sample);
        if (gt.isNoCall())           return VariantCallState.NoCall;
        else if (gt.isFiltered())    return VariantCallState.FilteredGenotype;
        else if ((gt.getGQ() != -1) && (gt.getGQ() < minGq))    return VariantCallState.LowGq;
        else if ((gt.getDP() != -1) && (gt.getDP() < minDp))    return VariantCallState.LowDp;
        else if (gt.isHet())         return VariantCallState.Het;
        else if (gt.isHomRef())      return VariantCallState.HomRef;
        else if (gt.isHomVar())      return VariantCallState.HomVar;

        throw new IllegalStateException("Could not classify variant: " + gt);
    }

    protected boolean doAllelesAgree(final VariantContext truthVariantContext, final VariantCallState truthCallState,
                                      final VariantContext callVariantContext, final VariantCallState callCallState) {
        boolean altAllelesAgree = true;
        if ((truthVariantContext != null) && (callVariantContext != null)) {
            if (!truthVariantContext.getReference().equals(callVariantContext.getReference())) {
                return false;       // This can happen for indels.
//                throw new PicardException("Ref alleles differ between VCF " + TRUTH_VCF.getAbsolutePath() + " Variant Context: " + truthVariantContext +
//                        " and " + CALL_VCF.getAbsolutePath() + " Variant Context: " + callVariantContext);
            }
            // Only checking for altAllele disagreement for HomVar and Het
            if (((truthCallState == VariantCallState.HomVar) && (callCallState == VariantCallState.HomVar)) ||
                    ((truthCallState == VariantCallState.Het) && (callCallState == VariantCallState.Het))) {
                // If here, either truth or call variant context have more than one alternate allele.  Check the genotype
                final List<Allele> truthSampleAlleles = truthVariantContext.getGenotype(TRUTH_SAMPLE).getAlleles();
                assert truthSampleAlleles.size() == 2;
                final List<Allele> callSampleAlleles = callVariantContext.getGenotype(CALL_SAMPLE).getAlleles();
                assert callSampleAlleles.size() == 2;

                final Allele truthAllele1 = truthSampleAlleles.get(0);
                final Allele callAllele1 = callSampleAlleles.get(0);

                if (truthCallState == VariantCallState.HomVar) {
                    altAllelesAgree = truthAllele1.equals(callAllele1);
                } else {
                    final Allele callAllele2 = callSampleAlleles.get(1);
                    final Allele truthAllele2 = truthSampleAlleles.get(1);
                    altAllelesAgree = ((truthAllele1.equals(callAllele1) && truthAllele2.equals(callAllele2)) ||
                            (truthAllele1.equals(callAllele2) && truthAllele2.equals(callAllele1)));
                }
            }
        }
        return altAllelesAgree;
    }


    public class GenotypeConcordanceVariantTypeDetailsMetrics extends MetricBase {
        /** The name of the 'truth' sample */
        public String TRUTH_SAMPLE;

        /** The name of the 'call' sample */
        public String CALL_SAMPLE;

        /** The type of the 'truth' sample (i.e. SNP, INDEL, NO_VARIATION) */
        public String TRUTH_TYPE;

        /** The type of the 'call' sample (i.e. SNP, INDEL, NO_VARIATION) */
        public String CALL_TYPE;

        /** The number of types of type TRUTH_STATE and CALL_STATE for the EVENT_TYPE and SAMPLEs */
        public long COUNT;
    }
}

/** Little class to hold a pair of VariantContexts that are in sync with one another. */
class VcTuple {
    public final VariantContext truthVariantContext;
    public final VariantContext callVariantContext;

    VcTuple(final VariantContext truthVariantContext, final VariantContext callVariantContext) {
        this.truthVariantContext = truthVariantContext;
        this.callVariantContext = callVariantContext;
    }
}

/** Iterator that takes a pair of iterators over VariantContexts and iterates over them in tandem. */
class PairedVariantSubContextIterator implements Iterator<VcTuple> {
    private final PeekableIterator<VariantContext> truthIterator;
    private final String truthSample;
    private final PeekableIterator<VariantContext> callIterator;
    private final String callSample;
    private final VariantContextComparator comparator;

    PairedVariantSubContextIterator(final Iterator<VariantContext> truthIterator, final String truthSample,
                                 final Iterator<VariantContext> callIterator, final String callSample,
                                 final SAMSequenceDictionary dict) {
        this.truthIterator = new PeekableIterator<VariantContext>(truthIterator);
        this.truthSample = truthSample;
        this.callIterator = new PeekableIterator<VariantContext>(callIterator);
        this.callSample = callSample;
        this.comparator = new VariantContextComparator(dict);
    }

    @Override
    public boolean hasNext() {
        return this.truthIterator.hasNext() || this.callIterator.hasNext();
    }

    @Override
    public VcTuple next() {
        if (!hasNext()) throw new IllegalStateException("next() called while hasNext() is false.");

        final VariantContext truthVariantContext = this.truthIterator.hasNext() ? this.truthIterator.peek() : null;
        final VariantContext callVariantContext = this.callIterator.hasNext() ? this.callIterator.peek() : null;

        // If one or the other is null because there is no next, just return a one-sided tuple
        if (truthVariantContext == null) {
            return new VcTuple(null, this.callIterator.next().subContextFromSample(callSample));
        }
        else if (callVariantContext == null) {
            return new VcTuple(this.truthIterator.next().subContextFromSample(truthSample), null);
        }

        // Otherwise check the ordering and do the right thing
        final int ordering = this.comparator.compare(truthVariantContext, callVariantContext);
        if (ordering == 0) {
            return new VcTuple(this.truthIterator.next().subContextFromSample(truthSample), this.callIterator.next().subContextFromSample(callSample));
        }
        else if (ordering < 0) {
            return new VcTuple(this.truthIterator.next().subContextFromSample(truthSample), null);
        }
        else {
            return new VcTuple(null, this.callIterator.next().subContextFromSample(callSample));
        }
    }

    @Override public void remove() {
        throw new UnsupportedOperationException();
    }
}




