package picard.vcf;

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.Set;

public class MakeVcfSampleNameMapTest extends CommandLineProgramTest {
    private static final Path TEST_RESOURCE_DIR = Paths.get("testdata/picard/vcf/MakeVcfSampleNameMap");

    public String getCommandLineProgramName() {
        return MakeVcfSampleNameMap.class.getSimpleName();
    }

    @Test
    public void testMakeMap() throws IOException {
        final File out = File.createTempFile("MakeMap", "test.sample_map");
        out.deleteOnExit();

        final String[] args = new String[]{
                "I=" + TEST_RESOURCE_DIR.resolve("NA12878.vcf").toString(),
                "I=" + TEST_RESOURCE_DIR.resolve("NA12891.vcf").toString(),
                "I=" + TEST_RESOURCE_DIR.resolve("NA12892.vcf").toString(),
                "O=" + out.getAbsolutePath()
        };

        final MakeVcfSampleNameMap clp = new MakeVcfSampleNameMap();
        Assert.assertEquals(clp.instanceMain(args), 0);

        final Set<String> expectedLines =
                new HashSet<>(Files.readAllLines(TEST_RESOURCE_DIR.resolve("expected.sample_map")));
        final Set<String> actualLines = new HashSet<>(Files.readAllLines(out.toPath()));

        Assert.assertEquals(actualLines, expectedLines);
    }

    @Test(expectedExceptions = PicardException.class)
    public void testFailMultiSampleVcf() {
        final String[] args = new String[]{
                "I=" + TEST_RESOURCE_DIR.getParent().resolve("mini.vcf"),
                "O=" + TEST_RESOURCE_DIR.resolve("not-used").toString()
        };

        final MakeVcfSampleNameMap clp = new MakeVcfSampleNameMap();
        clp.instanceMain(args);
    }
}
