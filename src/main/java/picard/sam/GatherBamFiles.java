package picard.sam;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;
import java.util.List;

/**
 * Program to perform a rapid "gather" operation on BAM files after a scatter operations where
 * the same process has been performed on different regions of a BAM file creating many smaller
 * BAM files that now need to be concatenated back together.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = GatherBamFiles.USAGE_SUMMARY + GatherBamFiles.USAGE_DETAILS,
        oneLineSummary = GatherBamFiles.USAGE_SUMMARY,
        programGroup = SamOrBam.class
)
public class GatherBamFiles extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Concatenate one or more BAM files as efficiently as possible";
    static final String USAGE_DETAILS = "This tool performs a rapid \"gather\" operation on BAM files after scatter" +
            " operations where the same process has been performed on different regions of a BAM file creating many " +
            "smaller BAM files that now need to be concatenated (reassembled) back together." +
            "<br /><br />" +
            "Assumes that the list of BAM files provided as INPUT are in the order that they should be concatenated and" +
            " simply concatenates the bodies of the BAM files while retaining the header from the first file.  " +
            "Operates via copying of the gzip blocks directly for speed but also supports generation of an MD5 on the" +
            " output and indexing of the output BAM file. Only supports BAM files, does not support SAM files." +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar GatherBamFiles \\<br /> " +
            "     I=input1.bam \\ <br /> " +
            "     I=input2.bam \\ <br /> " +
            "     O=gathered_files.bam" +
            "</pre> " +
            "<hr />";
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "Two or more BAM files or text files containing lists of BAM files (one per line).")
    public List<File> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output BAM file to write.")
    public File OUTPUT;

    private static final Log log = Log.getInstance(GatherBamFiles.class);

    // Stock main method.
    public static void main(final String[] args) {
        final GatherBamFiles gatherer = new GatherBamFiles();
        gatherer.CREATE_INDEX = true;
        gatherer.instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        final List<File> inputs = IOUtil.unrollFiles(INPUT, BamFileIoUtils.BAM_FILE_EXTENSION, ".sam");
        for (final File f : inputs) IOUtil.assertFileIsReadable(f);
        IOUtil.assertFileIsWritable(OUTPUT);

        if (determineBlockCopyingStatus(inputs)) {
            BamFileIoUtils.gatherWithBlockCopying(inputs, OUTPUT, CREATE_INDEX, CREATE_MD5_FILE);
        } else {
            gatherNormally(inputs, OUTPUT, CREATE_INDEX, CREATE_MD5_FILE, REFERENCE_SEQUENCE);
        }

        return 0;
    }

    private boolean determineBlockCopyingStatus(final List<File> inputs) {
        boolean useBlockCopying = true;
        for (final File f : inputs) {
            if (!BamFileIoUtils.isBamFile(f)) {
                useBlockCopying = false;
            }
        }
        return useBlockCopying;
    }

    /**
     * Simple implementation of a gather operations that uses SAMFileReaders and Writers in order to concatenate
     * multiple BAM files.
     */
    private static void gatherNormally(final List<File> inputs, final File output, final boolean createIndex, final boolean createMd5,
                                       final File referenceFasta) {
        final SAMFileHeader header;
        {
            header = SamReaderFactory.makeDefault().referenceSequence(referenceFasta).getFileHeader(inputs.get(0));
        }

        final SAMFileWriter out = new SAMFileWriterFactory().setCreateIndex(createIndex).setCreateMd5File(createMd5).makeSAMOrBAMWriter(header, true, output);

        for (final File f : inputs) {
            log.info("Gathering " + f.getAbsolutePath());
            final SamReader in = SamReaderFactory.makeDefault().referenceSequence(referenceFasta).open(f);
            for (final SAMRecord rec : in) out.addAlignment(rec);
            CloserUtil.close(in);
        }

        out.close();
    }

}
