package picard.analysis;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.annotation.Gff3Codec;
import picard.annotation.GtfFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

@CommandLineProgramProperties(
        summary = ValidateGff.USAGE_SUMMARY + ValidateGff.USAGE_DETAILS,
        oneLineSummary = ValidateGff.USAGE_SUMMARY,
        programGroup = DiagnosticsAndQCProgramGroup.class
)

public class ValidateGff extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Validate Gff file";
    static final String USAGE_DETAILS = "Validate Gff file";

    @Argument(shortName = "G", doc="Gff file")
    public File GFF_FILE;

    @Override
    protected final int doWork() {
        try (AbstractFeatureReader<GtfFeature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(GFF_FILE.getAbsolutePath(), null, new Gff3Codec(), false)) {
            final Iterator<GtfFeature> iterator = reader.iterator();
            while (iterator.hasNext()) {
                iterator.next();
            }
        } catch (IOException e) {
            throw new PicardException("Error reading GFF3 file " + GFF_FILE.getAbsolutePath());
        }
        return 0;
    }
}
