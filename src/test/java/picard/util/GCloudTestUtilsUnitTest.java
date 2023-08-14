package picard.util;

import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.stream.Stream;

public class GCloudTestUtilsUnitTest {

    @Test(groups = "bucket")
    public void testDownload() throws IOException {
        final String bigTextPath = GCloudTestUtils.getTestInputPath() + "nio/big.txt";
        try(final Stream<String> lines = Files.lines(IOUtil.getPath(bigTextPath))){
            String firstLine = lines.findFirst().orElseThrow();
            Assert.assertEquals(firstLine, "The Project Gutenberg EBook of The Adventures of Sherlock Holmes");
        }
    }

    @Test(groups = "bucket")
    public void testUpload() throws IOException {
        final Path uploadDir = Files.createTempDirectory(IOUtil.getPath(GCloudTestUtils.getTestStaging()), "picardTest");
        final Path txtUpload = uploadDir.resolve("tmp.txt");
        try {
            Assert.assertFalse(Files.exists(txtUpload));
            Files.writeString(txtUpload, "Hello there.");
            Assert.assertEquals(Files.readString(txtUpload), "Hello there.");
        } finally {
            Files.delete(txtUpload);
        }
        Assert.assertFalse(Files.exists(txtUpload));
        Assert.assertFalse(Files.exists(uploadDir));
    }
}