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
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.util.List;

/**
 * Concatenate efficiently BAM files that resulted from a scattered parallel analysis.
 * <p>
 * This tool performs a rapid "gather" or concatenation on BAM files. 
 * This is often needed in operations that have been run in parallel across genomics regions by scattering 
 * their execution across computing nodes and cores thus resulting in smaller BAM files.
 * </p>
 * <p>
 * This tool does not support SAM files.
 * </p>
 * <h3>Inputs</h3>
 * <p>
 * A list of BAM files to combine using the <code>INPUT</code> argument. 
 * These files must be provided in the order that they should be concatenated.
 * </p>
 *<h3>Output</h3>
 * <p>
 * A single BAM file. The header is copied from the first input file. 
 * </p>
 * <h3>Usage example:</h3>
 * <pre>
 * java -jar picard.jar GatherBamFiles \
 *      I=input1.bam \
 *      I=input2.bam \
 *      O=gathered_files.bam
 * </pre>
 * <h3>Notes</h3>
 * <p>
 * Operates via copying of the gzip blocks directly for speed but also supports generation of an MD5 on the
 * output and indexing of the output BAM file. 
 * </p>
 * 
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = "<p>" + GatherBamFiles.USAGE_SUMMARY + ".</p>" + GatherBamFiles.USAGE_DETAILS,
        oneLineSummary = GatherBamFiles.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature
public class GatherBamFiles extends CommandLineProgram {
    static final String USAGE_SUMMARY =
            "Concatenate efficiently BAM files that resulted from a scattered parallel analysis";
    static final String USAGE_DETAILS =
            "<p>This tool performs a rapid \"gather\" or concatenation on BAM files. " +
            "This is often needed in operations that have been run in parallel across genomics regions by scattering " +
            "their execution across computing nodes and cores thus resulting in smaller BAM files.</p><p>This tool does not support SAM files</p>" +
            "<h3>Inputs</h3>" +
            "<p>A list of BAM files to combine using the INPUT argument. " +
            "These files must be provided in the order that they should be concatenated.</p>" +
            "<h3>Output</h3>" +
            "<p>A single BAM file. The header is copied from the first input file.</p>" +
            "<h3>Usage example:</h3>" +
            "<pre>java -jar picard.jar GatherBamFiles \\\n" +
            "      I=input1.bam \\\n" +
            "      I=input2.bam \\\n" +
            "      O=gathered_files.bam</pre>" +
            "<h3>Notes</h3>" +
            "<p>Operates via copying of the gzip blocks directly for speed but also supports generation of an MD5 " +
            "on the output and indexing of the output BAM file.</p>" +
            "<hr/>";
    
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "Two or more BAM files or text files containing lists of BAM files (one per line).")
    public List<File> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output BAM file to write to.")
    public File OUTPUT;

    private static final Log log = Log.getInstance(GatherBamFiles.class);

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
