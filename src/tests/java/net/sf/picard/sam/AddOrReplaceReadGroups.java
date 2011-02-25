package net.sf.picard.sam;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.*;
import net.sf.samtools.SAMFileHeader.SortOrder;

import javax.xml.bind.SchemaOutputResolver;
import java.io.File;
import java.util.Arrays;

/**
 * Replaces read groups in a BAM file
 *
 * @author mdepristo
 */
public class AddOrReplaceReadGroups extends CommandLineProgram {
    @Usage(programVersion="1.0")
    public String USAGE = "Replaces all read groups in the INPUT file with a new read group and assigns " +
                          "all reads to this read group in the OUTPUT BAM";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input file (bam or sam).")
    public File INPUT = null;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file (bam or sam).")
    public File OUTPUT = null;

    @Option(shortName=StandardOptionDefinitions.SORT_ORDER_SHORT_NAME, optional=true,
            doc="Optional sort order to output in. If not supplied OUTPUT is in the same order as INPUT.")
    public SortOrder SORT_ORDER;

    @Option(shortName="ID",doc="Read Group ID")
    public String RGID = "1";

    @Option(shortName="LB",doc="Read Group Library")
    public String RGLB;

    @Option(shortName="PL",doc="Read Group platform (e.g. illumina, solid)")
    public String RGPL;

    @Option(shortName="PU",doc="Read Group platform unit (eg. run barcode)")
    public String RGPU;

    @Option(shortName="SM",doc="Read Group sample name")
    public String RGSM;

    @Option(shortName="CN", doc="Read Group sequencing center name", optional=true)
    public String RGCN;

    @Option(shortName="DS", doc="Read Group description", optional=true)
    public String RGDS;

    private final Log log = Log.getInstance(AddOrReplaceReadGroups.class);

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        new AddOrReplaceReadGroups().instanceMainWithExit(argv);
    }

    protected int doWork() {
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);

        SAMFileReader in = new SAMFileReader(INPUT);

        // create the read group we'll be using
        SAMReadGroupRecord rg = new SAMReadGroupRecord(RGID);
        rg.setLibrary(RGLB);
        rg.setPlatform(RGPL);
        rg.setSample(RGSM);
        rg.setPlatformUnit(RGPU);
        if (RGCN != null) rg.setSequencingCenter(RGCN);
        if (RGDS != null) rg.setDescription(RGDS);

        log.info(String.format("Created read group ID=%s PL=%s LB=%s SM=%s%n", rg.getId(), rg.getPlatform(), rg.getLibrary(), rg.getSample()));

        // create the new header and output file
        final SAMFileHeader inHeader = in.getFileHeader();
        final SAMFileHeader outHeader = inHeader.clone();
        outHeader.setReadGroups(Arrays.asList(rg));
        if (SORT_ORDER != null) outHeader.setSortOrder(SORT_ORDER);

        final SAMFileWriter outWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(outHeader,
                                                                                      outHeader.getSortOrder() == inHeader.getSortOrder(),
                                                                                      OUTPUT);

        for (final SAMRecord read : in) {
            read.setAttribute(SAMTag.RG.name(), RGID);
            outWriter.addAlignment(read);
        }

        // cleanup
        in.close();
        outWriter.close();
        return 0;
    }
}
