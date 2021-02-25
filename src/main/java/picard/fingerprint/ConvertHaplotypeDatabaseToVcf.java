package picard.fingerprint;

import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.io.FileNotFoundException;

@CommandLineProgramProperties(
        summary = "Convert Haplotype database file to vcf",
        oneLineSummary = "Convert Haplotype database file to vcf",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
public class ConvertHaplotypeDatabaseToVcf extends CommandLineProgram {

    @Argument(doc = "Haplotype database to be converted to VCF.", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "Where to write converted haplotype database VCF.", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Override
    protected boolean requiresReference() {
        return true;
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(INPUT);
        
        final HaplotypeMap haplotypeMap = new HaplotypeMap(INPUT);

        try {
            haplotypeMap.writeAsVcf(OUTPUT, REFERENCE_SEQUENCE);
        } catch (final FileNotFoundException ex) {
            throw new PicardException("uggh", ex);
        }

        return 0;
    }
}
