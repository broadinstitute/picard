package picard.analysis.artifacts;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.*;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.analysis.SinglePassSamProgram;
import picard.analysis.artifacts.SequencingArtifactMetrics.BaitBiasDetailMetrics;
import picard.analysis.artifacts.SequencingArtifactMetrics.BaitBiasSummaryMetrics;
import picard.analysis.artifacts.SequencingArtifactMetrics.PreAdapterDetailMetrics;
import picard.analysis.artifacts.SequencingArtifactMetrics.PreAdapterSummaryMetrics;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import picard.util.DbSnpBitSetUtil;
import picard.util.VariantType;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

import static htsjdk.samtools.util.CodeUtil.getOrElse;
import static picard.cmdline.StandardOptionDefinitions.MINIMUM_MAPPING_QUALITY_SHORT_NAME;

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
        summary = CollectSequencingArtifactMetrics.USAGE_SUMMARY + CollectSequencingArtifactMetrics.USAGE_DETAILS,
        oneLineSummary = CollectSequencingArtifactMetrics.USAGE_SUMMARY,
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@DocumentedFeature
public class CollectSequencingArtifactMetrics extends SinglePassSamProgram {
public static final String USAGE_SUMMARY = "Collect metrics to quantify single-base sequencing artifacts.  ";
static final String USAGE_DETAILS = "<p>This tool examines two sources of sequencing errors associated with hybrid selection "+
"protocols.  These errors are divided into two broad categories, pre-adapter and bait-bias.  Pre-adapter errors can arise from "+
"laboratory manipulations of a nucleic acid sample e.g. shearing and occur prior to the ligation of adapters for PCR "+
"amplification (hence the name pre-adapter).  </p>" +
    
"<p>Bait-bias artifacts occur during or after the target selection step, and correlate with substitution rates that are "+
"'biased', or higher for sites having one base on the reference/positive strand relative to sites having the complementary "+
"base on that strand.  For example, during the target selection step, a (G>T) artifact might result in a higher substitution "+
"rate at sites with a G on the positive strand (and C on the negative), relative to sites with the flip (C positive)/(G negative)." +
"  This is known as the 'G-Ref' artifact. </p>" +
"" +
"<p>For additional information on these types of artifacts, please see the corresponding GATK dictionary entries on "+
"<a href='https://www.broadinstitute.org/gatk/guide/article?id=6333'>bait-bias</a> and "+
"<a href='https://www.broadinstitute.org/gatk/guide/article?id=6332'>pre-adapter artifacts</a>.</p>"+
""+
"<p>This tool produces four files; summary and detail metrics files for both pre-adapter and bait-bias artifacts. The detailed "+
"metrics show the error rates for each type of base substitution within every possible triplet base configuration.  Error " +
"rates associated with these substitutions are Phred-scaled and provided as quality scores, the lower the value, the more " +
"likely it is that an alternate base call is due to an artifact. The summary metrics provide likelihood information on the " +
"'worst-case' errors. </p>" +
""+
"<h4>Usage example:</h4>" +
"<pre>" +
"java -jar picard.jar CollectSequencingArtifactMetrics \\<br />" +
"     I=input.bam \\<br />" +
"     O=artifact_metrics.txt \\<br />" +
"     R=reference_sequence.fasta" +
"</pre>" +

"Please see the metrics at the following links " +
"<a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#SequencingArtifactMetrics.PreAdapterDetailMetrics'>PreAdapterDetailMetrics</a>, "+
"<a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#SequencingArtifactMetrics.PreAdapterSummaryMetrics'>PreAdapterSummaryMetrics</a>, "+
"<a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#SequencingArtifactMetrics.BaitBiasDetailMetrics'>BaitBiasDetailMetrics</a>, and "+
"<a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#SequencingArtifactMetrics.BaitBiasSummaryMetrics'>BaitBiasSummaryMetrics</a> "+
"for complete descriptions of the output metrics produced by this tool. "+
"<hr />"
;
    @Argument(doc = "An optional list of intervals to restrict analysis to.", optional = true)
    public File INTERVALS;

    @Argument(doc = "VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis.", optional = true)
    public File DB_SNP;

    @Argument(shortName = "Q", doc = "The minimum base quality score for a base to be included in analysis.")
    public int MINIMUM_QUALITY_SCORE = 20;

    @Argument(shortName = MINIMUM_MAPPING_QUALITY_SHORT_NAME, doc = "The minimum mapping quality score for a base to be included in analysis.")
    public int MINIMUM_MAPPING_QUALITY = 30;

    @Argument(shortName = "MIN_INS", doc = "The minimum insert size for a read to be included in analysis.")
    public int MINIMUM_INSERT_SIZE = 60;

    @Argument(shortName = "MAX_INS", doc = "The maximum insert size for a read to be included in analysis. Set to 0 to have no maximum.")
    public int MAXIMUM_INSERT_SIZE = 600;

    @Argument(shortName = "UNPAIRED", doc = "Include unpaired reads. If set to true then all paired reads will be included as well - " +
            "MINIMUM_INSERT_SIZE and MAXIMUM_INSERT_SIZE will be ignored.")
    public boolean INCLUDE_UNPAIRED = false;

    @Argument(shortName = "DUPES", doc = "Include duplicate reads. If set to true then all reads flagged as duplicates will be included as well.")
    public boolean INCLUDE_DUPLICATES = false;

    @Argument(shortName = "NON_PF", doc = "Whether or not to include non-PF reads.")
    public boolean INCLUDE_NON_PF_READS = false;

    @Argument(shortName = "TANDEM", doc = "Set to true if mate pairs are being sequenced from the same strand, " +
            "i.e. they're expected to face the same direction.")
    public boolean TANDEM_READS = false;

    @Argument(doc = "When available, use original quality scores for filtering.")
    public boolean USE_OQ = true;

    @Argument(doc = "The number of context bases to include on each side of the assayed base.")
    public int CONTEXT_SIZE = 1;

    @Argument(doc = "If specified, only print results for these contexts in the detail metrics output. " +
                  "However, the summary metrics output will still take all contexts into consideration.", optional = true)
    public Set<String> CONTEXTS_TO_PRINT = new HashSet<>();

    @Argument(shortName = "EXT", doc="Append the given file extension to all metric file names (ex. OUTPUT.pre_adapter_summary_metrics.EXT). None if null", optional=true)
    public String FILE_EXTENSION = null;

    private static final String UNKNOWN_LIBRARY = "UnknownLibrary";
    private static final String UNKNOWN_SAMPLE = "UnknownSample";

    private File preAdapterSummaryOut;
    private File preAdapterDetailsOut;
    private File baitBiasSummaryOut;
    private File baitBiasDetailsOut;
    private File errorSummaryFile;

    private IntervalListReferenceSequenceMask intervalMask;
    private DbSnpBitSetUtil dbSnpMask;
    private SamRecordFilter recordFilter;

    private String currentRefString = null;
    private int currentRefIndex = -1;

    private final Set<String> samples = new HashSet<>();
    private final Set<String> libraries = new HashSet<>();
    private final Map<String, ArtifactCounter> artifactCounters = new HashMap<>();

    private static final Log log = Log.getInstance(CollectSequencingArtifactMetrics.class);

    @Override
    protected boolean requiresReference() {
        return true;
    }

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> messages = new ArrayList<>();

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
        final String outext = (null != FILE_EXTENSION) ? FILE_EXTENSION : ""; // Add a file extension if desired
        preAdapterSummaryOut = new File(OUTPUT + SequencingArtifactMetrics.PRE_ADAPTER_SUMMARY_EXT + outext);
        preAdapterDetailsOut = new File(OUTPUT + SequencingArtifactMetrics.PRE_ADAPTER_DETAILS_EXT + outext);
        baitBiasSummaryOut   = new File(OUTPUT + SequencingArtifactMetrics.BAIT_BIAS_SUMMARY_EXT + outext);
        baitBiasDetailsOut   = new File(OUTPUT + SequencingArtifactMetrics.BAIT_BIAS_DETAILS_EXT + outext);
        errorSummaryFile     = new File(OUTPUT + SequencingArtifactMetrics.ERROR_SUMMARY_EXT + outext);
        IOUtil.assertFilesAreWritable(Arrays.asList(preAdapterSummaryOut, preAdapterDetailsOut, baitBiasSummaryOut, baitBiasDetailsOut, errorSummaryFile));

        for (final SAMReadGroupRecord rec : header.getReadGroups()) {
            samples.add(getOrElse(rec.getSample(), UNKNOWN_SAMPLE));
            libraries.add(getOrElse(rec.getLibrary(), UNKNOWN_LIBRARY));
        }

        if (INTERVALS != null) {
            IOUtil.assertFileIsReadable(INTERVALS);

            final IntervalList intervalList = IntervalList.fromFile(INTERVALS).uniqued();
            intervalMask = new IntervalListReferenceSequenceMask(intervalList);

            if (DB_SNP != null) {
                IOUtil.assertFileIsReadable(DB_SNP);
                dbSnpMask = new DbSnpBitSetUtil(DB_SNP, header.getSequenceDictionary(), EnumSet.noneOf(VariantType.class), intervalList, Optional.of(log));
            }
        }
        else if (DB_SNP != null) {
            IOUtil.assertFileIsReadable(DB_SNP);
            dbSnpMask = new DbSnpBitSetUtil(DB_SNP, header.getSequenceDictionary(), EnumSet.noneOf(VariantType.class), null, Optional.of(log));
        }

        // set record-level filters
        final List<SamRecordFilter> filters = new ArrayList<>();
        if (!INCLUDE_NON_PF_READS) filters.add(new FailsVendorReadQualityFilter());
        filters.add(new NotPrimaryAlignmentFilter());
        if (!INCLUDE_DUPLICATES) filters.add(new DuplicateReadFilter());
        filters.add(new AlignedFilter(true)); // discard unmapped reads
        filters.add(new MappingQualityFilter(MINIMUM_MAPPING_QUALITY));
        if (!INCLUDE_UNPAIRED) {
            final int effectiveMaxInsertSize = (MAXIMUM_INSERT_SIZE == 0) ? Integer.MAX_VALUE : MAXIMUM_INSERT_SIZE;
            filters.add(new InsertSizeFilter(MINIMUM_INSERT_SIZE, effectiveMaxInsertSize));
        }
        recordFilter = new AggregateFilter(filters);

        // set up the artifact counters
        final String sampleAlias = StringUtil.join(",", new ArrayList<>(samples));
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

                // skip non-ACGT bases
                if (!SequenceUtil.isUpperACGTN((byte)readBase))
                    continue;

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
        final MetricsFile<ErrorSummaryMetrics,?> errorSummaryMetricsFile = getMetricsFile();

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

        // Calculate the summary error rates - it's CRITICAL that the other files are written out
        // first as this code modifies the pre-adapter detail metrics!
        if (!preAdapterDetailMetricsFile.getMetrics().isEmpty()) {
            final List<PreAdapterDetailMetrics> in = preAdapterDetailMetricsFile.getMetrics();
            in.forEach(m -> {
                if (m.REF_BASE == 'G' || m.REF_BASE == 'T') {
                    m.REF_BASE = (char) SequenceUtil.complement((byte) m.REF_BASE);
                    m.ALT_BASE = (char) SequenceUtil.complement((byte) m.ALT_BASE);
                }
            });

            // Group the metrics by error type
            final Map<String,List<PreAdapterDetailMetrics>> byError =
                    in.stream().collect(Collectors.groupingBy(m -> m.REF_BASE + ">" + m.ALT_BASE));

            for (final String error : new TreeSet<>(byError.keySet())) {
                final List<PreAdapterDetailMetrics> ms = byError.get(error);
                final ErrorSummaryMetrics summary = new ErrorSummaryMetrics();
                summary.REF_BASE = ms.get(0).REF_BASE;
                summary.ALT_BASE = ms.get(0).ALT_BASE;
                summary.SUBSTITUTION = error;
                summary.REF_COUNT = ms.stream().mapToLong(m -> m.PRO_REF_BASES + m.CON_REF_BASES).sum();
                summary.ALT_COUNT = ms.stream().mapToLong(m -> m.PRO_ALT_BASES + m.CON_ALT_BASES).sum();
                summary.calculateDerivedFields();
                errorSummaryMetricsFile.addMetric(summary);
            }
        }

        errorSummaryMetricsFile.write(errorSummaryFile);
    }

    @Override
    protected boolean usesNoRefReads() { return false; }
}
