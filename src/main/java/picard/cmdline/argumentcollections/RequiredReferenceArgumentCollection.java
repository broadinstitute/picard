package picard.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;

/**
 * Argument collection for references that are required (and not common).
 */
public class RequiredReferenceArgumentCollection implements ReferenceArgumentCollection {

    @Argument(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence file.", common = false)
    public File REFERENCE_SEQUENCE;

    /**
     * @return The reference provided by the user.
     */
    public File getReferenceFile() {
        return REFERENCE_SEQUENCE;
    };
}
