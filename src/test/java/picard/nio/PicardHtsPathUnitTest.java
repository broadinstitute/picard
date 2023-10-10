package picard.nio;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class PicardHtsPathUnitTest {
    @DataProvider(name="ReplaceExtensionTestInput")
    public Object[][] getReplaceExtensionTestInput() {
        return new Object[][] {
                {"gs://hellbender/test/resources/hg19mini.fasta", ".fai", false, "gs://hellbender/test/resources/hg19mini.fai"},
                {"gs://hellbender/test/resources/hg19mini.fasta", ".fai", true, "gs://hellbender/test/resources/hg19mini.fasta.fai"},
                {"testdata/picard/sam/test.bam", ".bai", false, new File("testdata/picard/sam/test.bai").getAbsolutePath()}, // a relative input will be converted to absolute input
                {"testdata/picard/sam/test.bam", ".bai", true, new File("testdata/picard/sam/test.bam.bai").getAbsolutePath()},
                {"/Users/jdoe/workspace/picard/testdata/picard/sam/test.bam", ".bai", false, "/Users/jdoe/workspace/picard/testdata/picard/sam/test.bai"},
                {"/Users/jdoe/workspace/picard/testdata/picard/sam/test.bam", ".md5", true, "/Users/jdoe/workspace/picard/testdata/picard/sam/test.bam.md5"}
        };
    }

    @Test(dataProvider = "ReplaceExtensionTestInput")
    public void testReplaceExtension(final String originalURI, final String newExtension, final boolean append,
                                     final String expectedURI){
        final PicardHtsPath originalPath = new PicardHtsPath(originalURI);
        final PicardHtsPath newPath = PicardHtsPath.replaceExtension(originalPath, newExtension, append);
        // We cannot directly compare the PicardHtsPath because currently replaceExtension() takes an PicardHtsPath with a
        // relative rawInputString and return an object with an absolute rawInputString. Instead, we check that the absolute URIs match.
        Assert.assertEquals(newPath.getURIString(), new PicardHtsPath(expectedURI).getURIString());
    }
}