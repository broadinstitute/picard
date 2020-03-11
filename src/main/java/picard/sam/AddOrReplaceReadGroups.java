package picard.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Iso8601Date;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Assigns all the reads in a file to a single new read-group.
 *
 * <h3>Summary</h3>
 * Many tools (Picard and GATK for example) require or assume the presence of at least one <code>RG</code> tag, defining a "read-group"
 * to which each read can be assigned (as specified in the <code>RG</code> tag in the SAM record).
 * This tool enables the user to assign all the reads in the {@link #INPUT} to a single new read-group.
 * For more information about read-groups, see the <a href='https://www.broadinstitute.org/gatk/guide/article?id=6472'>
 * GATK Dictionary entry.</a>
 * <br />
 * This tool accepts as INPUT BAM and SAM files or URLs from the
 * <a href="http://ga4gh.org/#/documentation">Global Alliance for Genomics and Health (GA4GH)</a>.
 * <h3>Usage example:</h3>
 * <pre>
 * java -jar picard.jar AddOrReplaceReadGroups \
 *       I=input.bam \
 *       O=output.bam \
 *       RGID=4 \
 *       RGLB=lib1 \
 *       RGPL=ILLUMINA \
 *       RGPU=unit1 \
 *       RGSM=20
 * </pre>
 * <h3>Caveats</h3>
 * The value of the tags must adhere (according to the <a href="https://samtools.github.io/hts-specs/SAMv1.pdf">SAM-spec</a>)
 * with the regex <pre>{@value #READGROUP_ID_REGEX}</pre> (one or more characters from the ASCII range 32 through 126). In
 * particular <code>&lt;Space&gt;</code> is the only non-printing character allowed.
 * <br/>
 * The program enables only the wholesale assignment of all the reads in the {@link #INPUT} to a single read-group. If your file
 * already has reads assigned to multiple read-groups, the original <code>RG</code> value will be lost.
 *
 * @author mdepristo
 */
@CommandLineProgramProperties(
        summary = AddOrReplaceReadGroups.USAGE_SUMMARY + AddOrReplaceReadGroups.USAGE_DETAILS,
        oneLineSummary = AddOrReplaceReadGroups.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature
public class AddOrReplaceReadGroups extends CommandLineProgram {
    static final public String READGROUP_ID_REGEX="^[ -~]+$";

    static final String USAGE_SUMMARY = "Assigns all the reads in a file to a single new read-group.";
    static final String USAGE_DETAILS =
            "\n\nThis tool accepts INPUT BAM and SAM files or URLs from the <a href=\"http://ga4gh.org/#/documentation\">Global Alliance for Genomics and Health (GA4GH)</a>.\n" +
            "<h3>Usage example:</h3>" +
            "\n"+
            "java -jar picard.jar AddOrReplaceReadGroups \\\n" +
            "      I=input.bam \\\n" +
            "      O=output.bam \\\n" +
            "      RGID=4 \\\n" +
            "      RGLB=lib1 \\\n" +
            "      RGPL=ILLUMINA \\\n" +
            "      RGPU=unit1 \\\n" +
            "      RGSM=20\n " +
            "\n" +
            "<h3>Caveats</h3>\n" +
            "The value of the tags must adhere (according to the <a href=\"https://samtools.github.io/hts-specs/SAMv1.pdf\">SAM-spec</a>) " +
            "with the regex <code>'" + READGROUP_ID_REGEX + "'</code> (one or more characters from the ASCII range 32 through 126). " +
                    "In particular &lt;Space&gt; is the only non-printing character allowed.\n" +
            "\n" +
            "The program enables only the wholesale assignment of all the reads in the INPUT to a single read-group. If your file " +
            "already has reads assigned to multiple read-groups, the original RG value will be lost. \n\n" +
            "For more information about read-groups, see the <a href='https://www.broadinstitute.org/gatk/guide/article?id=6472'>" +
            "GATK Dictionary entry.</a>";

    @Argument(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input file (BAM or SAM or a GA4GH url).")
    public String INPUT = null;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (BAM or SAM).")
    public File OUTPUT = null;

    @Argument(shortName = StandardOptionDefinitions.SORT_ORDER_SHORT_NAME, optional = true,
            doc = "Optional sort order to output in. If not supplied OUTPUT is in the same order as INPUT.")
    public SortOrder SORT_ORDER;

    @Argument(shortName = "ID", doc = "Read-Group ID")
    public String RGID = "1";

    @Argument(shortName = "LB", doc = "Read-Group library")
    public String RGLB;

    @Argument(shortName = "PL", doc = "Read-Group platform (e.g. ILLUMINA, SOLID)")
    public String RGPL;

    @Argument(shortName = "PU", doc = "Read-Group platform unit (eg. run barcode)")
    public String RGPU;

    @Argument(shortName = "SM", doc = "Read-Group sample name")
    public String RGSM;

    @Argument(shortName = "CN", doc = "Read-Group sequencing center name", optional = true)
    public String RGCN;

    @Argument(shortName = "DS", doc = "Read-Group description", optional = true)
    public String RGDS;

    @Argument(shortName = "DT", doc = "Read-Group run date", optional = true)
    public Iso8601Date RGDT;

    @Argument(shortName = "KS", doc = "Read-Group key sequence", optional = true)
    public String RGKS;

    @Argument(shortName = "FO", doc = "Read-Group flow order", optional = true)
    public String RGFO;

    @Argument(shortName = "PI", doc = "Read-Group predicted insert size", optional = true)
    public Integer RGPI;

    @Argument(shortName = "PG", doc = "Read-Group program group", optional = true)
    public String RGPG;
    
    @Argument(shortName = "PM", doc = "Read-Group platform model", optional = true)
    public String RGPM;

    private final Log log = Log.getInstance(AddOrReplaceReadGroups.class);

    protected int doWork() {
        IOUtil.assertInputIsValid(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final SamReader in = SamReaderFactory.makeDefault()
            .referenceSequence(REFERENCE_SEQUENCE)
            .open(SamInputResource.of(INPUT));

        // create the read-group we'll be using
        final SAMReadGroupRecord rg = new SAMReadGroupRecord(RGID);
        rg.setLibrary(RGLB);
        rg.setPlatform(RGPL);
        rg.setSample(RGSM);
        rg.setPlatformUnit(RGPU);
        if (RGCN != null) rg.setSequencingCenter(RGCN);
        if (RGDS != null) rg.setDescription(RGDS);
        if (RGDT != null) rg.setRunDate(RGDT);
        if (RGPI != null) rg.setPredictedMedianInsertSize(RGPI);
        if (RGPG != null) rg.setProgramGroup(RGPG);
        if (RGPM != null) rg.setPlatformModel(RGPM);
        if (RGKS != null) rg.setKeySequence(RGKS);
        if (RGFO != null) rg.setFlowOrder(RGFO);

        log.info(String.format("Created read-group ID=%s PL=%s LB=%s SM=%s%n", rg.getId(), rg.getPlatform(), rg.getLibrary(), rg.getSample()));

        // create the new header and output file
        final SAMFileHeader inHeader = in.getFileHeader();
        final SAMFileHeader outHeader = inHeader.clone();
        outHeader.setReadGroups(Collections.singletonList(rg));
        if (SORT_ORDER != null) outHeader.setSortOrder(SORT_ORDER);

        final SAMFileWriter outWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(outHeader,
                outHeader.getSortOrder() == inHeader.getSortOrder(),
                OUTPUT);

        final ProgressLogger progress = new ProgressLogger(log);
        for (final SAMRecord read : in) {
            read.setAttribute(SAMTag.RG.name(), RGID);
            outWriter.addAlignment(read);
            progress.record(read);
        }

        // cleanup
        CloserUtil.close(in);
        outWriter.close();
        return 0;
    }

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> validationFailures = new ArrayList<>();

        checkTagValue("RGID", RGID).ifPresent(validationFailures::add);
        checkTagValue("RGLB", RGLB).ifPresent(validationFailures::add);
        checkTagValue("RGPL", RGPL).ifPresent(validationFailures::add);
        checkTagValue("RGPU", RGPU).ifPresent(validationFailures::add);
        checkTagValue("RGSM", RGSM).ifPresent(validationFailures::add);
        checkTagValue("RGCN", RGCN).ifPresent(validationFailures::add);
        checkTagValue("RGDS", RGDS).ifPresent(validationFailures::add);
        checkTagValue("RGKS", RGKS).ifPresent(validationFailures::add);
        checkTagValue("RGFO", RGFO).ifPresent(validationFailures::add);
        checkTagValue("RGPG", RGPG).ifPresent(validationFailures::add);
        checkTagValue("RGPM", RGPM).ifPresent(validationFailures::add);

        if (!validationFailures.isEmpty()) {
            return validationFailures.toArray(new String[validationFailures.size()]);
        }

        return super.customCommandLineValidation();
    }

    private final Pattern pattern = Pattern.compile(READGROUP_ID_REGEX);

    private Optional<String> checkTagValue(final String tagName, final String value) {
        if (value == null) {
            return Optional.empty();
        }

        final Matcher matcher = pattern.matcher(value);

        if (matcher.matches()) {
            return Optional.empty();
        } else {
            return Optional.of(String.format("The values of tags in a SAM header must adhere to the regular expression '%s'," +
                    "but the value provided for %s, '%s', doesn't.",READGROUP_ID_REGEX, tagName, value));
        }
    }
}
