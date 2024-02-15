package picard;

import com.intel.gkl.compression.IntelDeflater;
import com.intel.gkl.compression.IntelDeflaterFactory;
import com.intel.gkl.compression.IntelInflater;
import com.intel.gkl.compression.IntelInflaterFactory;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.zip.DeflaterFactory;
import htsjdk.samtools.util.zip.InflaterFactory;
import org.apache.commons.lang3.SystemUtils;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.programgroups.OtherProgramGroup;
import picard.cmdline.programgroups.Testing;

import java.lang.reflect.Field;

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

    @Test
    public void testIntelInflaterIsUsed(){
        final InflaterDeflaterTester cmd = new InflaterDeflaterTester();
        cmd.instanceMain(new String[]{});
        Assert.assertEquals(cmd.inflaterFactory.getClass(), IntelInflaterFactory.class);
    }

    @Test
    public void testDeflaterIsUsed(){
        final InflaterDeflaterTester cmd = new InflaterDeflaterTester();
        cmd.instanceMain(new String[]{});
        Assert.assertEquals(cmd.deflaterFactory.getClass(), IntelDeflaterFactory.class);
    }

    private void checkIntelSupported(final String componentName) {
        if (!SystemUtils.IS_OS_LINUX && !SystemUtils.IS_OS_MAC) {
            throw new SkipException(componentName + " is not available on this platform");
        }

        if (SystemUtils.OS_ARCH != null && SystemUtils.OS_ARCH.equals("ppc64le")) {
            throw new SkipException(componentName + " is not available for this architecture");
        }
    }


    @CommandLineProgramProperties(summary = "test program for checking if the intel optimized inflater/deflater are active",
            oneLineSummary = "test program please ignore",
            programGroup = Testing.class,
            omitFromCommandLine = true)
    public static class InflaterDeflaterTester extends CommandLineProgram {
        public InflaterFactory inflaterFactory;
        public DeflaterFactory deflaterFactory;

        @Override
        protected int doWork() {
            final SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
            inflaterFactory = getFieldValue(samReaderFactory, "inflaterFactory", InflaterFactory.class);

            final SAMFileWriterFactory samFileWriterFactory = new SAMFileWriterFactory();
            deflaterFactory = getFieldValue(samFileWriterFactory, "deflaterFactory", DeflaterFactory.class);

            return 0;
        }

        private <T,R> R getFieldValue(final T obj,final String fieldName,  Class<R> clazz) {
            try {
                final Field deflaterFactoryField = obj.getClass().getDeclaredField(fieldName);
                deflaterFactoryField.setAccessible(true);
                return clazz.cast(deflaterFactoryField.get(obj));
            } catch (NoSuchFieldException | IllegalAccessException e) {
                throw new RuntimeException(e);
            }
        }
    }
}
