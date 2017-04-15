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
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;
import java.util.Arrays;

/**
 * Replaces read groups in a BAM file
 *
 * @author mdepristo
 */
@CommandLineProgramProperties(
        summary = AddOrReplaceReadGroups.USAGE_SUMMARY + AddOrReplaceReadGroups.USAGE_DETAILS,
        oneLineSummary = AddOrReplaceReadGroups.USAGE_SUMMARY,
        programGroup = SamOrBam.class
)
public class AddOrReplaceReadGroups extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Replace read groups in a BAM file.";
    static final String USAGE_DETAILS = "This tool enables the user to replace all read groups in the INPUT file with a single new read " +
            "group and assign all reads to this read group in the OUTPUT BAM file.<br /><br />" +
            "For more information about read groups, see the <a href='https://www.broadinstitute.org/gatk/guide/article?id=6472'>" +
            "GATK Dictionary entry.</a> <br /><br /> " +
            "This tool accepts INPUT BAM and SAM files or URLs from the Global Alliance for Genomics and Health (GA4GH) (see http://ga4gh.org/#/documentation)." +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar AddOrReplaceReadGroups \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=output.bam \\<br />" +
            "      RGID=4 \\<br />" +
            "      RGLB=lib1 \\<br />" +
            "      RGPL=illumina \\<br />" +
            "      RGPU=unit1 \\<br />" +
            "      RGSM=20" +
            "</pre>" +
            "<hr />" ;
    @Argument(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input file (BAM or SAM or a GA4GH url).")
    public String INPUT = null;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (BAM or SAM).")
    public File OUTPUT = null;

    @Argument(shortName = StandardOptionDefinitions.SORT_ORDER_SHORT_NAME, optional = true,
            doc = "Optional sort order to output in. If not supplied OUTPUT is in the same order as INPUT.")
    public SortOrder SORT_ORDER;

    @Argument(shortName = "ID", doc = "Read Group ID")
    public String RGID = "1";

    @Argument(shortName = "LB", doc = "Read Group library")
    public String RGLB;

    @Argument(shortName = "PL", doc = "Read Group platform (e.g. illumina, solid)")
    public String RGPL;

    @Argument(shortName = "PU", doc = "Read Group platform unit (eg. run barcode)")
    public String RGPU;

    @Argument(shortName = "SM", doc = "Read Group sample name")
    public String RGSM;

    @Argument(shortName = "CN", doc = "Read Group sequencing center name", optional = true)
    public String RGCN;

    @Argument(shortName = "DS", doc = "Read Group description", optional = true)
    public String RGDS;

    @Argument(shortName = "DT", doc = "Read Group run date", optional = true)
    public Iso8601Date RGDT;

    @Argument(shortName = "PI", doc = "Read Group predicted insert size", optional = true)
    public Integer RGPI;
    
    @Argument(shortName = "PG", doc = "Read Group program group", optional = true)
    public String RGPG;
    
    @Argument(shortName = "PM", doc = "Read Group platform model", optional = true)
    public String RGPM;

    private final Log log = Log.getInstance(AddOrReplaceReadGroups.class);

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        new AddOrReplaceReadGroups().instanceMainWithExit(argv);
    }

    protected int doWork() {
        IOUtil.assertInputIsValid(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final SamReader in = SamReaderFactory.makeDefault()
            .referenceSequence(REFERENCE_SEQUENCE)
            .open(SamInputResource.of(INPUT));

        // create the read group we'll be using
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

        log.info(String.format("Created read group ID=%s PL=%s LB=%s SM=%s%n", rg.getId(), rg.getPlatform(), rg.getLibrary(), rg.getSample()));

        // create the new header and output file
        final SAMFileHeader inHeader = in.getFileHeader();
        final SAMFileHeader outHeader = inHeader.clone();
        outHeader.setReadGroups(Arrays.asList(rg));
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
}
