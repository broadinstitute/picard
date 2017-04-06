package picard.fingerprint;

import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.programgroups.Fingerprinting;

/**
 * Created by farjoun on 4/5/17.
 */
/**
 * Program to check that all fingerprints within the set of input files appear to come from the same
 * individual. Can be used to cross-check readgroups, libraries, samples, or files.
 *
 * @author Tim Fennell
 * @author Yossi Farjoun
 */
@CommandLineProgramProperties(
        usage = "Checks if all fingerprints within a set of files appear to come from the same individual.",
        usageShort = "Checks if all fingerprints appear to come from the same individual.",
        programGroup = Fingerprinting.class
)
@Deprecated // use CrosscheckFingerprints
public class CrosscheckReadGroupFingerprints extends CrosscheckFingerprints {
}
