package picard.cmdline.argumentcollections;

import htsjdk.samtools.Defaults;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;

/**
 * Picard default argument collection for an optional reference.
 */
public class OptionalReferenceArgumentCollection implements ReferenceArgumentCollection {

    @Argument(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence file.", common = true, optional = true)
    public File REFERENCE_SEQUENCE = Defaults.REFERENCE_FASTA;

    /**
     * @return The reference provided by the user, if any, or the default defined by {@code htsjdk.samtools.Defaults.REFERENCE_FASTA}
     */
    @Override
    public File getReferenceFile() {
        return REFERENCE_SEQUENCE;
    };
}
