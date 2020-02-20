package picard.fingerprint;

import org.testng.Assert;

import java.util.Arrays;
import java.util.Objects;

public class FingerprintingTestUtils {

    public static boolean areHaplotypeProbabilitiesEqual(final HaplotypeProbabilities lhs, final HaplotypeProbabilities rhs) {
        if (lhs == rhs) {
            return true;
        }
        if (lhs == null || rhs == null || lhs.getClass() != rhs.getClass()) {
            return false;
        }

        if (!Objects.equals(lhs.getHaplotype(), rhs.getHaplotype())) {
            return false;
        }

        return Arrays.equals(lhs.getLikelihoods(), rhs.getLikelihoods());
    }

    public static void assertHaplotypeProbabilitiesEqual(final HaplotypeProbabilities lhs, final HaplotypeProbabilities rhs) {
        Assert.assertTrue(areHaplotypeProbabilitiesEqual(lhs, rhs),
                "Expected HaplotypeProbabilities to be equal, but they differ: " +
                        lhs + ", " + rhs);
    }

    public static void assertFingerPrintHPsAreEqual(final Fingerprint lhs, final Fingerprint rhs) {
        Assert.assertEquals(lhs.keySet().size(), rhs.keySet().size());

        for (final HaplotypeBlock block : lhs.keySet()) {
            Assert.assertTrue(lhs.containsKey(block), "HaplotypeBlock was missing from lhs" + block);
            Assert.assertTrue(rhs.containsKey(block), "HaplotypeBlock was missing from rhs" + block);
            FingerprintingTestUtils.assertHaplotypeProbabilitiesEqual(lhs.get(block), rhs.get(block));
        }
    }
}