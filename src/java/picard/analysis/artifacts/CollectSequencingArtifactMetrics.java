package picard.analysis.artifacts;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.FailsVendorReadQualityFilter;
import htsjdk.samtools.filter.InsertSizeFilter;
import htsjdk.samtools.filter.MappingQualityFilter;
import htsjdk.samtools.filter.NotPrimaryAlignmentFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.IntervalListReferenceSequenceMask;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.analysis.SinglePassSamProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;
import picard.util.DbSnpBitSetUtil;
import picard.analysis.artifacts.SequencingArtifactMetrics.*;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static htsjdk.samtools.util.CodeUtil.getOrElse;

/**
 * Quantify substitution errors caused by mismatched base pairings during various
 * stages of sample / library prep.
 *
 * We measure two distinct error types - artifacts that are introduced before
 * the addition of the read1/read2 adapters ("pre adapter") and those that are
 * introduced after target selection ("bait bias"). For each of these, we provide
 * summary metrics as well as detail metrics broken down by reference context
 * (the ref bases surrounding the substitution event).
 *
 * For a deeper explanation, see Costello et al. 2013:
 * http://www.ncbi.nlm.nih.gov/pubmed/23303777
 *
 * @author mattsooknah
 *
 */
@CommandLineProgramProperties(
        usage = CollectSequencingArtifactMetrics.USAGE_SUMMARY + CollectSequencingArtifactMetrics.USAGE_DETAILS,
        usageShort = CollectSequencingArtifactMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectSequencingArtifactMetrics extends SinglePassSamProgram {
    static final String USAGE_SUMMARY = "Collect metrics to quantify single-base sequencing artifacts.";
    static final String USAGE_DETAILS = "This tool examines two sources of sequencing errors resulting from hybrid selection protocols:" +
            " <a href='https://www.broadinstitute.org/gatk/guide/article?id=6333'>bait-bias</a> and " +
            "<a href='https://www.broadinstitute.org/gatk/guide/article?id=6332'>" +
            "pre-adapter artifacts</a>. For a brief primer on these types of artifacts, see the corresponding GATK Dictionary entries." +
            "<br /><br />" +
            "This tool produces four files; summary and detail metrics files for both pre-adapter and bait-bias artifacts. The detailed " +
            "metrics show the error rates for each type of base substitution within every possible triplet base configuration.  Error " +
            "rates associated with these substitutions are Phred-scaled and provided as quality scores, the lower the value, the more " +
            "likely it is that an alternate base call is due to an artifact. The summary metrics provide likelihood information on the " +
            "\"worst-case\" errors. <br />" +
            "" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CollectSequencingArtifactMetrics \\<br />" +
            "     I=input.bam \\<br />" +
            "     O=artifact_metrics.txt \\<br />" +
            "     R=reference_sequence.fasta" +
            "</pre>" +
            "" +
            "For additional information, please see " +
            "<a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#SequencingArtifactMetrics.PreAdapterDetailMetrics'>the PreAdapterDetailMetrics documentation</a>, the " +
            "<a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#SequencingArtifactMetrics.PreAdapterSummaryMetrics'>the PreAdapterSummaryMetrics documentation</a>, the " +
            "<a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#SequencingArtifactMetrics.BaitBiasDetailMetrics'>the BaitBiasDetailMetrics documentation</a>, and the " +
            "<a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#SequencingArtifactMetrics.BaitBiasSummaryMetrics'>the BaitBiasSummaryMetrics documentation</a>. " +
            "<hr />" ;
    @Option(doc = "An optional list of intervals to restrict analysis to.", optional = true)
    public File INTERVALS;

    @Option(doc = "VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis.", optional = true)
    public File DB_SNP;

    @Option(shortName = "Q", doc = "The minimum base quality score for a base to be included in analysis.")
    public int MINIMUM_QUALITY_SCORE = 20;

    @Option(shortName = "MQ", doc = "The minimum mapping quality score for a base to be included in analysis.")
    public int MINIMUM_MAPPING_QUALITY = 30;

    @Option(shortName = "MIN_INS", doc = "The minimum insert size for a read to be included in analysis.")
    public int MINIMUM_INSERT_SIZE = 60;

    @Option(shortName = "MAX_INS", doc = "The maximum insert size for a read to be included in analysis. Set to 0 to have no maximum.")
    public int MAXIMUM_INSERT_SIZE = 600;

    @Option(shortName = "UNPAIRED", doc = "Include unpaired reads. If set to true then all paired reads will be included as well - " +
            "MINIMUM_INSERT_SIZE and MAXIMUM_INSERT_SIZE will be ignored.")
    public boolean INCLUDE_UNPAIRED = false;

    @Option(shortName = "TANDEM", doc = "Set to true if mate pairs are being sequenced from the same strand, " +
            "i.e. they're expected to face the same direction.")
    public boolean TANDEM_READS = false;

    @Option(doc = "When available, use original quality scores for filtering.")
    public boolean USE_OQ = true;

    @Option(doc = "The number of context bases to include on each side of the assayed base.")
    public int CONTEXT_SIZE = 1;

    @Option(doc = "If specified, only print results for these contexts in the detail metrics output. " +
                  "However, the summary metrics output will still take all contexts into consideration.")
    public Set<String> CONTEXTS_TO_PRINT = new HashSet<String>();

    private static final String UNKNOWN_LIBRARY = "UnknownLibrary";
    private static final String UNKNOWN_SAMPLE = "UnknownSample";

    private File preAdapterSummaryOut;
    private File preAdapterDetailsOut;
    private File baitBiasSummaryOut;
    private File baitBiasDetailsOut;

    private IntervalListReferenceSequenceMask intervalMask;
    private DbSnpBitSetUtil dbSnpMask;
    private SamRecordFilter recordFilter;

    private String currentRefString = null;
    private int currentRefIndex = -1;

    private final Set<String> samples = new HashSet<String>();
    private final Set<String> libraries = new HashSet<String>();
    private final Map<String, ArtifactCounter> artifactCounters = new HashMap<String, ArtifactCounter>();

    public static void main(final String[] args) {
        new CollectSequencingArtifactMetrics().instanceMainWithExit(args);
    }

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> messages = new ArrayList<String>();

        final int contextFullLength = 2 * CONTEXT_SIZE + 1;
        if (CONTEXT_SIZE < 0) messages.add("CONTEXT_SIZE cannot be negative");
        for (final String context : CONTEXTS_TO_PRINT) {
            if (context.length() != contextFullLength) {
                messages.add("Context " + context + " is not the length implied by CONTEXT_SIZE: " + contextFullLength);
            }
        }

        if (MINIMUM_INSERT_SIZE < 0) messages.add("MINIMUM_INSERT_SIZE cannot be negative");
        if (MAXIMUM_INSERT_SIZE < 0) messages.add("MAXIMUM_INSERT_SIZE cannot be negative");
        if (MAXIMUM_INSERT_SIZE > 0 && MAXIMUM_INSERT_SIZE < MINIMUM_INSERT_SIZE) {
            messages.add("MAXIMUM_INSERT_SIZE cannot be less than MINIMUM_INSERT_SIZE unless set to 0");
        }

        return messages.isEmpty() ? null : messages.toArray(new String[messages.size()]);
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        preAdapterSummaryOut = new File(OUTPUT + SequencingArtifactMetrics.PRE_ADAPTER_SUMMARY_EXT);
        preAdapterDetailsOut = new File(OUTPUT + SequencingArtifactMetrics.PRE_ADAPTER_DETAILS_EXT);
        baitBiasSummaryOut = new File(OUTPUT + SequencingArtifactMetrics.BAIT_BIAS_SUMMARY_EXT);
        baitBiasDetailsOut = new File(OUTPUT + SequencingArtifactMetrics.BAIT_BIAS_DETAILS_EXT);

        IOUtil.assertFileIsWritable(preAdapterSummaryOut);
        IOUtil.assertFileIsWritable(preAdapterDetailsOut);
        IOUtil.assertFileIsWritable(baitBiasSummaryOut);
        IOUtil.assertFileIsWritable(baitBiasDetailsOut);

        for (final SAMReadGroupRecord rec : header.getReadGroups()) {
            samples.add(getOrElse(rec.getSample(), UNKNOWN_SAMPLE));
            libraries.add(getOrElse(rec.getLibrary(), UNKNOWN_LIBRARY));
        }

        if (INTERVALS != null) {
            IOUtil.assertFileIsReadable(INTERVALS);
            intervalMask = new IntervalListReferenceSequenceMask(IntervalList.fromFile(INTERVALS).uniqued());
        }

        if (DB_SNP != null) {
            IOUtil.assertFileIsReadable(DB_SNP);
            dbSnpMask = new DbSnpBitSetUtil(DB_SNP, header.getSequenceDictionary());
        }

        // set record-level filters
        final List<SamRecordFilter> filters = new ArrayList<SamRecordFilter>();
        filters.add(new FailsVendorReadQualityFilter());
        filters.add(new NotPrimaryAlignmentFilter());
        filters.add(new DuplicateReadFilter());
        filters.add(new AlignedFilter(true)); // discard unmapped reads
        filters.add(new MappingQualityFilter(MINIMUM_MAPPING_QUALITY));
        if (!INCLUDE_UNPAIRED) {
            final int effectiveMaxInsertSize = (MAXIMUM_INSERT_SIZE == 0) ? Integer.MAX_VALUE : MAXIMUM_INSERT_SIZE;
            filters.add(new InsertSizeFilter(MINIMUM_INSERT_SIZE, effectiveMaxInsertSize));
        }
        recordFilter = new AggregateFilter(filters);

        // set up the artifact counters
        final String sampleAlias = StringUtil.join(",", new ArrayList<String>(samples));
        for (final String library : libraries) {
            artifactCounters.put(library, new ArtifactCounter(sampleAlias, library, CONTEXT_SIZE, TANDEM_READS));
        }
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        // see if the whole read should be skipped
        if (recordFilter.filterOut(rec)) return;

        // check read group + library
        final String library = (rec.getReadGroup() == null) ? UNKNOWN_LIBRARY : getOrElse(rec.getReadGroup().getLibrary(), UNKNOWN_LIBRARY);
        if (!libraries.contains(library)) {
            // should never happen if SAM is valid
            throw new PicardException("Record contains library that is missing from header: " + library);
        }

        // set up some constants that don't change in the loop below
        final int contextFullLength = 2 * CONTEXT_SIZE + 1;
        final ArtifactCounter counter = artifactCounters.get(library);
        final byte[] readBases = rec.getReadBases();
        final byte[] readQuals;
        if (USE_OQ) {
            final byte[] tmp = rec.getOriginalBaseQualities();
            readQuals = tmp == null ? rec.getBaseQualities() : tmp;
        } else {
            readQuals = rec.getBaseQualities();
        }

        // iterate over aligned positions
        for (final AlignmentBlock block : rec.getAlignmentBlocks()) {
            for (int offset = 0; offset < block.getLength(); offset++) {
                // remember, these are 1-based!
                final int readPos = block.getReadStart() + offset;
                final int refPos = block.getReferenceStart() + offset;

                // skip low BQ sites
                final byte qual = readQuals[readPos - 1];
                if (qual < MINIMUM_QUALITY_SCORE) continue;

                // skip N bases in read
                final char readBase = Character.toUpperCase((char)readBases[readPos - 1]);
                if (readBase == 'N') continue;

                /**
                 * Skip regions outside of intervals.
                 *
                 * NB: IntervalListReferenceSequenceMask.get() has side-effects which assume
                 * that successive ReferenceSequence's passed to this method will be in-order
                 * (e.g. it will break if you call acceptRead() with chr1, then chr2, then chr1
                 * again). So this only works if the underlying iteration is coordinate-sorted.
                 */
                if (intervalMask != null && !intervalMask.get(ref.getContigIndex(), refPos)) continue;

                // skip dbSNP sites
                if (dbSnpMask != null && dbSnpMask.isDbSnpSite(ref.getName(), refPos)) continue;

                // skip the ends of the reference
                final int contextStartIndex = refPos - CONTEXT_SIZE - 1;
                if (contextStartIndex < 0 || contextStartIndex + contextFullLength > ref.length()) continue;

                // skip contexts with N bases
                final String context = getRefContext(ref, contextStartIndex, contextFullLength);
                if (context.contains("N")) continue;

                // count the base!
                counter.countRecord(context, readBase, rec);
            }
        }
    }

    private String getRefContext(final ReferenceSequence ref, final int contextStartIndex, final int contextFullLength) {
        // cache the upper-cased string version of this reference so we don't need to create a string for every base in every read
        if (currentRefIndex != ref.getContigIndex()) {
            currentRefString = new String(ref.getBases()).toUpperCase();
            currentRefIndex = ref.getContigIndex();
        }
        return currentRefString.substring(contextStartIndex, contextStartIndex + contextFullLength);
    }

    @Override
    protected void finish() {
        final MetricsFile<PreAdapterSummaryMetrics, Integer> preAdapterSummaryMetricsFile = getMetricsFile();
        final MetricsFile<PreAdapterDetailMetrics, Integer> preAdapterDetailMetricsFile = getMetricsFile();
        final MetricsFile<BaitBiasSummaryMetrics, Integer> baitBiasSummaryMetricsFile = getMetricsFile();
        final MetricsFile<BaitBiasDetailMetrics, Integer> baitBiasDetailMetricsFile = getMetricsFile();

        for (final ArtifactCounter counter : artifactCounters.values()) {
            // build metrics
            counter.finish();

            // write metrics
            preAdapterSummaryMetricsFile.addAllMetrics(counter.getPreAdapterSummaryMetrics());
            baitBiasSummaryMetricsFile.addAllMetrics(counter.getBaitBiasSummaryMetrics());

            for (final PreAdapterDetailMetrics preAdapterDetailMetrics : counter.getPreAdapterDetailMetrics()) {
                if (CONTEXTS_TO_PRINT.isEmpty() || CONTEXTS_TO_PRINT.contains(preAdapterDetailMetrics.CONTEXT)) {
                    preAdapterDetailMetricsFile.addMetric(preAdapterDetailMetrics);
                }
            }
            for (final BaitBiasDetailMetrics baitBiasDetailMetrics : counter.getBaitBiasDetailMetrics()) {
                if (CONTEXTS_TO_PRINT.isEmpty() || CONTEXTS_TO_PRINT.contains(baitBiasDetailMetrics.CONTEXT)) {
                    baitBiasDetailMetricsFile.addMetric(baitBiasDetailMetrics);
                }
            }

        }

        preAdapterDetailMetricsFile.write(preAdapterDetailsOut);
        preAdapterSummaryMetricsFile.write(preAdapterSummaryOut);
        baitBiasDetailMetricsFile.write(baitBiasDetailsOut);
        baitBiasSummaryMetricsFile.write(baitBiasSummaryOut);
    }

    @Override
    protected boolean usesNoRefReads() { return false; }
}
