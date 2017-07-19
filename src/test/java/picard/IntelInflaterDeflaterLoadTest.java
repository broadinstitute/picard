package picard;

import com.intel.gkl.compression.IntelDeflater;
import com.intel.gkl.compression.IntelInflater;
import org.apache.commons.lang3.SystemUtils;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.Test;

/**
 * Test that the Intel Inflater and Deflater can be loaded successfully.
 */
public class IntelInflaterDeflaterLoadTest {
    @Test
    public void testIntelInflaterIsAvailable() {
        checkIntelSupported("IntelInflater");
        Assert.assertTrue(new IntelInflater().load(null),
                "Intel shared library was not loaded. This could be due to a configuration error, or your system might not support it.");
    }

    @Test
    public void testIntelDeflaterIsAvailable() {
        checkIntelSupported("IntelDeflater");
        Assert.assertTrue(new IntelDeflater().load(null),
                "Intel shared library was not loaded. This could be due to a configuration error, or your system might not support it.");
    }

    private void checkIntelSupported(final String componentName) {
        if (!SystemUtils.IS_OS_LINUX && !SystemUtils.IS_OS_MAC) {
            throw new SkipException(componentName + " is not available on this platform");
        }

        if (SystemUtils.OS_ARCH != null && SystemUtils.OS_ARCH.equals("ppc64le")) {
            throw new SkipException(componentName + " is not available for this architecture");
        }
    }
}
