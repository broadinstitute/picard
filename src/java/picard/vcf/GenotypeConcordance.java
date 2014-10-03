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
import picard.vcf.GenotypeConcordanceStates.*;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import static htsjdk.variant.variantcontext.VariantContext.Type.*;

/**
 * Calculates the concordance between genotype data for two samples in two different VCFs - one being considered the truth (or reference)
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
                    "Calculates the concordance between genotype data for two samples in two different VCFs - one being considered the truth (or reference) " +
                    "the other being considered the call.  The concordance is broken into separate results sections for SNPs and indels.  Summary and detailed statistics are reported\n\n" +
                    "Note that for any pair of variants to compare, only the alleles for the samples under interrogation are considered " +
                    "and MNP, Symbolic, and Mixed classes of variants are not included.";

    @Option(shortName = "TV", doc="The VCF containing the truth sample")
    public File TRUTH_VCF;

    @Option(shortName = "CV", doc="The VCF containing the call sample")
    public File CALL_VCF;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Basename for the two metrics files that are to be written." +
            " Resulting files will be <OUTPUT>" + SUMMARY_METRICS_FILE_EXTENSION + "  and <OUTPUT>" + DETAILED_METRICS_FILE_EXTENSION + ".")
    public File OUTPUT;

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

    @Option(doc="If true, use the VCF index, else iterate over the entire VCF.", optional = true)
    public boolean USE_VCF_INDEX = false;

    private final Log log = Log.getInstance(GenotypeConcordance.class);
    private final ProgressLogger progress = new ProgressLogger(log, 10000, "checked", "variants");

    public static final String SUMMARY_METRICS_FILE_EXTENSION = ".summary_metrics.txt";
    public static final String DETAILED_METRICS_FILE_EXTENSION = ".detailed_metrics.txt";

    protected GenotypeConcordanceCounts snpCounter;
    public GenotypeConcordanceCounts getSnpCounter() { return snpCounter; }

    protected GenotypeConcordanceCounts indelCounter;
    public GenotypeConcordanceCounts getIndelCounter() { return indelCounter; }

    // TODO: add optimization if the samples are in the same file
    // TODO: add option for auto-detect pairs based on same sample name
    // TODO: allow multiple sample-pairs in one pass

    public static void main(final String[] args) {
        new GenotypeConcordance().instanceMainWithExit(args);
    }

    @Override protected int doWork() {
        IOUtil.assertFileIsReadable(TRUTH_VCF);
        IOUtil.assertFileIsReadable(CALL_VCF);
        final File summaryMetricsFile = new File(OUTPUT + SUMMARY_METRICS_FILE_EXTENSION);
        final File detailedMetricsFile = new File(OUTPUT + DETAILED_METRICS_FILE_EXTENSION);
        IOUtil.assertFileIsWritable(summaryMetricsFile);
        IOUtil.assertFileIsWritable(detailedMetricsFile);

        final boolean usingIntervals = this.INTERVALS != null && this.INTERVALS.size() > 0;
        IntervalList intervals = null;
        SAMSequenceDictionary intervalsSamSequenceDictionary = null;
        if (usingIntervals) {
            log.info("Loading up region lists.");
            long genomeBaseCount = 0;
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
            }
            log.info("Finished loading up region lists.");
        }

        final VCFFileReader truthReader = new VCFFileReader(TRUTH_VCF, USE_VCF_INDEX);
        final VCFFileReader callReader = new VCFFileReader(CALL_VCF, USE_VCF_INDEX);

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
        final PairedVariantSubContextIterator pairedIterator = new PairedVariantSubContextIterator(truthIterator, TRUTH_SAMPLE, callIterator, CALL_SAMPLE, truthReader.getFileHeader().getSequenceDictionary());
        snpCounter   = new GenotypeConcordanceCounts();
        indelCounter = new GenotypeConcordanceCounts();

        // A map to keep track of the count of Truth/Call States which we could not successfully classify
        final Map<String, Integer> unClassifiedStatesMap = new HashMap<String, Integer>();
        int numMismatchedRefAlleles = 0;

        log.info("Starting iteration over variants.");
        while (pairedIterator.hasNext()) {
            final VcTuple tuple = pairedIterator.next();

            final VariantContext.Type truthVariantContextType;
            final Allele truthRef;
            if (tuple.truthVariantContext != null) {
                truthVariantContextType = tuple.truthVariantContext.getType();
                truthRef = tuple.truthVariantContext.getReference();
            } else {
                truthVariantContextType = NO_VARIATION;
                truthRef = null;
            }
            final VariantContext.Type callVariantContextType;
            final Allele callRef;
            if (tuple.callVariantContext != null) {
                callVariantContextType = tuple.callVariantContext.getType();
                callRef = tuple.callVariantContext.getReference();
            } else {
                callVariantContextType = NO_VARIATION;
                callRef = null;
            }

            // A flag to keep track of whether we have been able to successfully classify the Truth/Call States.
            // Unclassified include MIXED/MNP/Symbolic...
            boolean stateClassified = false;
            boolean misMatchedRefAlleles = false;
            if (truthRef != null && callRef != null) {
                if (!truthRef.equals(callRef)) {
                    log.warn("Ref alleles disagree between: " + tuple.truthVariantContext + " and " + tuple.callVariantContext + " - skipping");
                    misMatchedRefAlleles = true;
                    numMismatchedRefAlleles++;
                }
            }
            if (!misMatchedRefAlleles) {
                final TruthAndCallStates truthAndCallStates = determineState(tuple.truthVariantContext, TRUTH_SAMPLE, tuple.callVariantContext, CALL_SAMPLE, MIN_GQ, MIN_DP);
                if (truthVariantContextType == SNP) {
                    if ((callVariantContextType == SNP) || (callVariantContextType == MIXED) || (callVariantContextType == NO_VARIATION)) {
                        // Note.  If truth is SNP and call is MIXED, the event will be logged in the indelCounter, with row = MIXED
                        snpCounter.increment(truthAndCallStates);
                        stateClassified = true;
                    }
                }
                else if (truthVariantContextType == INDEL) {
                    // Note.  If truth is Indel and call is MIXED, the event will be logged in the indelCounter, with row = MIXED
                    if ((callVariantContextType == INDEL) || (callVariantContextType == MIXED) || (callVariantContextType == NO_VARIATION)) {
                        indelCounter.increment(truthAndCallStates);
                        stateClassified = true;
                    }
                }
                else if (truthVariantContextType == MIXED) {
                    // Note.  If truth is MIXED and call is SNP, the event will be logged in the snpCounter, with column = MIXED
                    if (callVariantContextType == SNP) {
                        snpCounter.increment(truthAndCallStates);
                        stateClassified = true;
                    }
                    // Note.  If truth is MIXED and call is INDEL, the event will be logged in the snpCounter, with column = MIXED
                    else if (callVariantContextType == INDEL) {
                        indelCounter.increment(truthAndCallStates);
                        stateClassified = true;
                    }
                }
                else if (truthVariantContextType == NO_VARIATION) {
                    if (callVariantContextType == SNP) {
                        snpCounter.increment(truthAndCallStates);
                        stateClassified = true;
                    }
                    else if (callVariantContextType == INDEL) {
                        indelCounter.increment(truthAndCallStates);
                        stateClassified = true;
                    }
                }
                if (!stateClassified) {
                    final String condition = truthVariantContextType + " " + callVariantContextType;
                    Integer count = unClassifiedStatesMap.get(condition);
                    if (count == null) count = 0;
                    unClassifiedStatesMap.put(condition, ++count);
                }
            }

            final VariantContext variantContextForLogging = tuple.truthVariantContext != null ? tuple.truthVariantContext : tuple.callVariantContext;
            progress.record(variantContextForLogging.getChr(), variantContextForLogging.getStart());
        }

        // Calculate and store the summary-level metrics
        final MetricsFile<GenotypeConcordanceSummaryMetrics,?> genotypeConcordanceSummaryMetricsFile = getMetricsFile();
        GenotypeConcordanceSummaryMetrics summaryMetrics = new GenotypeConcordanceSummaryMetrics(SNP, snpCounter, TRUTH_SAMPLE, CALL_SAMPLE);
        genotypeConcordanceSummaryMetricsFile.addMetric(summaryMetrics);
        summaryMetrics = new GenotypeConcordanceSummaryMetrics(INDEL, indelCounter, TRUTH_SAMPLE, CALL_SAMPLE);
        genotypeConcordanceSummaryMetricsFile.addMetric(summaryMetrics);
        genotypeConcordanceSummaryMetricsFile.write(summaryMetricsFile);

        // Calculate and store the detailed metrics for both SNP and indels
        final MetricsFile<GenotypeConcordanceDetailMetrics,?> genotypeConcordanceDetailMetrics = getMetricsFile();
        outputDetailMetricsFile(SNP, genotypeConcordanceDetailMetrics, snpCounter, TRUTH_SAMPLE, CALL_SAMPLE);
        outputDetailMetricsFile(INDEL, genotypeConcordanceDetailMetrics, indelCounter, TRUTH_SAMPLE, CALL_SAMPLE);
        genotypeConcordanceDetailMetrics.write(detailedMetricsFile);

        for (final String condition : unClassifiedStatesMap.keySet()) {
            log.info("Uncovered truth/call Variant Context Type Counts: " + condition + " " + unClassifiedStatesMap.get(condition));
        }
        if (numMismatchedRefAlleles > 0) {
            log.warn("Number of variant comparisons with mismatched REF alleles: " + numMismatchedRefAlleles);
        }

        return 0;
    }

    /**
    * Outputs the detailed statistics tables for SNP and Indel match categories.
    **/
    private void outputDetailMetricsFile(final VariantContext.Type variantType, final MetricsFile<GenotypeConcordanceDetailMetrics,?> genotypeConcordanceDetailMetricsFile,
                                         final GenotypeConcordanceCounts counter, final String truthSampleName, final String callSampleName) {
        for (final TruthState truthState : TruthState.values()) {
            for (final CallState callState : CallState.values()) {
                final int count = counter.getCount(truthState, callState);
                if (count > 0 || OUTPUT_ALL_ROWS) {
                    final GenotypeConcordanceDetailMetrics detailMetrics = new GenotypeConcordanceDetailMetrics();
                    detailMetrics.VARIANT_TYPE = variantType;
                    detailMetrics.TRUTH_SAMPLE = truthSampleName;
                    detailMetrics.CALL_SAMPLE = callSampleName;
                    detailMetrics.TRUTH_STATE = truthState;
                    detailMetrics.CALL_STATE = callState;
                    detailMetrics.COUNT = count;
                    genotypeConcordanceDetailMetricsFile.addMetric(detailMetrics);
                }
            }
        }
    }

    /**
     * A method to determine the truth and call states for a pair of variant contexts representing truth and call.
     * A variety of variant and genotype-level checks are first used to determine if either of the the variant contexts
     * are filtered and after that a comparison of the called genotype alleles to determine appropriate truth and call state
     *
     * Note that this method does NOT check for SNP versus Indel.  It is assumed that that check is done by the caller and the results
     * of this method are stored by SNP/Indel.
     * Note that if a variant context has BOTH GQ and DP less than the specified threshold, then it will be of Truth/Call State LOW_GQ
     *
     * @param truthContext A variant context representing truth
     * @param truthSample The name of the truth sample
     * @param callContext A variant context representing the call
     * @param callSample The name of the call sample
     * @param minGq Threshold for filtering by genotype attribute GQ
     * @param minDp Threshold for filtering by genotype attribute DP
     * @return TruthAndCallStates object containing the TruthState and CallState determined here.
     */
    final TruthAndCallStates determineState(final VariantContext truthContext, final String truthSample, final VariantContext callContext, final String callSample, final int minGq, final int minDp) {
        TruthState truthState = null;
        CallState callState = null;

        // TODO: what about getPloidy()

        Genotype truthGenotype = null, callGenotype = null;

        // Site level checks
        if (truthContext == null) truthState = TruthState.MISSING;
        else if (truthContext.isMixed()) truthState = TruthState.IS_MIXED;
        else if (truthContext.isFiltered()) truthState = TruthState.VC_FILTERED;
        else {
            // Genotype level checks
            truthGenotype = truthContext.getGenotype(truthSample);
            if (truthGenotype.isNoCall())           truthState = TruthState.NO_CALL;
            else if (truthGenotype.isFiltered())    truthState = TruthState.GT_FILTERED;
            else if ((truthGenotype.getGQ() != -1) && (truthGenotype.getGQ() < minGq)) truthState = TruthState.LOW_GQ;
            else if ((truthGenotype.getDP() != -1) && (truthGenotype.getDP() < minDp)) truthState = TruthState.LOW_DP;
            // Note.  Genotype.isMixed means that it is called on one chromosome and NOT on the other
            else if ((truthGenotype.isMixed())) truthState = TruthState.NO_CALL;
        }

        // Site level checks
        if (callContext == null) callState = CallState.MISSING;
        else if (callContext.isMixed()) callState = CallState.IS_MIXED;
        else if (callContext.isFiltered()) callState = CallState.VC_FILTERED;
        else {
            // Genotype level checks
            callGenotype = callContext.getGenotype(callSample);
            if (callGenotype.isNoCall())           callState = CallState.NO_CALL;
            else if (callGenotype.isFiltered())    callState = CallState.GT_FILTERED;
            else if ((callGenotype.getGQ() != -1) && (callGenotype.getGQ() < minGq)) callState = CallState.LOW_GQ;
            else if ((callGenotype.getDP() != -1) && (callGenotype.getDP() < minDp)) callState = CallState.LOW_DP;
                // Note.  Genotype.isMixed means that it is called on one chromosome and NOT on the other
            else if ((callGenotype.isMixed())) callState = CallState.NO_CALL;
        }

        // initialize the reference
        final Allele truthRef = (truthContext != null) ? truthContext.getReference() : null;
        final Allele callRef  = (callContext != null) ?  callContext.getReference() : null;
        if (truthRef != null && callRef != null) {
            if (!truthRef.equals(callRef)) {
                throw new IllegalStateException("Ref alleles mismatch between: " + truthContext + " and " + callContext);
            }
        }

        final OrderedSet<Allele> allAlleles = new OrderedSet<Allele>();

        if (truthContext != null || callContext != null) {
            allAlleles.smartAdd(truthContext == null ? callContext.getReference() : truthContext.getReference()); // zeroth allele;
        }

        if (null == truthState) {
            if (truthGenotype.getAlleles().size() != 2) {
                throw new IllegalStateException("Genotype for Variant Context: " + truthContext + " does not have exactly 2 alleles");
            }
            allAlleles.smartAdd(truthGenotype.getAllele(0));
            allAlleles.smartAdd(truthGenotype.getAllele(1));
        }

        /**
         *  if either of the call alleles is in allAlleles, with index > 1, we need to make sure that allele has index 1.
         *  this is because of the following situations:
         *
         *      REF TRUTH   CALL-GT TRUTH-STATE     CALL-STATE
         *      A   C/G     C/A     HET_VAR1_VAR2   HET_REF_VAR1
         *      A   G/C     C/A     HET_VAR1_VAR2   HET_REF_VAR1
         *      A   G/C     G/A     HET_VAR1_VAR2   HET_REF_VAR1
         *      A   G/C     G/A     HET_VAR1_VAR2   HET_REF_VAR1
         *
         *  so, in effect, the order of the alleles in the TRUTH doesn't determine the assignment of allele to Var1 and Var2,
         *  only once the call is known can this assignment be made.
         */

        if (null == callState) {
            if (callGenotype.getAlleles().size() != 2) {
                throw new IllegalStateException("Genotype for Variant Context: " + callContext + " does not have exactly 2 alleles");
            }
            if (allAlleles.indexOf(callGenotype.getAllele(0)) > 1 || allAlleles.indexOf(callGenotype.getAllele(1)) > 1) {
                allAlleles.remove(2);
                allAlleles.remove(1);
                allAlleles.smartAdd(truthGenotype.getAllele(1));
                allAlleles.smartAdd(truthGenotype.getAllele(0));
            }

            allAlleles.smartAdd(callGenotype.getAllele(0));
            allAlleles.smartAdd(callGenotype.getAllele(1));
        }

        // Truth
        if (null == truthState) {
            final int allele0idx = allAlleles.indexOf(truthGenotype.getAllele(0));
            final int allele1idx = allAlleles.indexOf(truthGenotype.getAllele(1));

            if (allele0idx == allele1idx) { //HOM
                truthState = TruthState.getHom(allele0idx);
            } else { //HET
                truthState = TruthState.getVar(allele0idx, allele1idx);
            }
        }

        // Call
        if (null == callState) {
            final int allele0idx = allAlleles.indexOf(callGenotype.getAllele(0));
            final int allele1idx = allAlleles.indexOf(callGenotype.getAllele(1));

            if (allele0idx == allele1idx) { //HOM
                callState = CallState.getHom(allele0idx);
            } else { //HET
                callState = CallState.getHet(allele0idx, allele1idx);
            }

            if (null == callState) {
                throw new IllegalStateException("This should never happen...  Could not classify the call variant: " + callGenotype);
            }
        }

        return new TruthAndCallStates(truthState, callState);
    }
}

/** like a list, but if you ask for an index of an item, it will first add that item.
also, same item cannot be added more than once (like a set)
 */
class OrderedSet<T> extends ArrayList<T> {

    public int smartIndexOf(final T o) {
        smartAdd(o);
        return super.indexOf(o);
    }

    public boolean smartAdd(final T o) {
        if (!this.contains(o)) {
            return add(o);
        }
        return false;
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




