package picard.nio;

import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.Files;

public class HttpNioIntegrationTest {
    @Test
    public void testCanReadFromHttpsPath() throws IOException {
        final String theWholeReadme = Files.readString(IOUtil.getPath("https://raw.githubusercontent.com/broadinstitute/picard/master/README.md"));
        Assert.assertTrue(theWholeReadme.contains("Picard"));
    }

}
