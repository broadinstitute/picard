package picard.nio;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class PicardHtsPathUnitTest {
    @DataProvider(name="ReplaceExtensionTestInput")
    public Object[][] getReplaceExtensionTestInput() {
        return new Object[][] {
                {"gs://hellbender/test/resources/hg19mini.fasta", ".fai", false, "gs://hellbender/test/resources/hg19mini.fai"},
                {"gs://hellbender/test/resources/hg19mini.fasta", ".fai", true, "gs://hellbender/test/resources/hg19mini.fasta.fai"},
                {"testdata/picard/sam/test.bam", ".bai", false, "testdata/picard/sam/test.bai"}, // relative paths
                {"testdata/picard/sam/test.bam", ".bai", true, "testdata/picard/sam/test.bam.bai"},
                {"/Users/jsoto/workspace/picard/testdata/picard/sam/test.bam", ".bai", false, "/Users/jsoto/workspace/picard/testdata/picard/sam/test.bai"},
                {"/Users/jsoto/workspace/picard/testdata/picard/sam/test.bam", ".md5", true, "/Users/jsoto/workspace/picard/testdata/picard/sam/test.bam.md5"}
        };
    }

    @Test(dataProvider = "ReplaceExtensionTestInput")
    public void testReplaceExtension(final String originalURI, final String newExtension, final boolean append,
                                     final String expectedURI){
        final PicardHtsPath originalPath = new PicardHtsPath(originalURI);
        final PicardHtsPath newPath = PicardHtsPath.replaceExtension(originalPath, newExtension, append);
        // Assert.assertEquals(newPath, new PicardHtsPath(expectedURL)); // tsato: this would pass if we change HtsPath::equals() not to compare rawInputStrings.
        Assert.assertEquals(newPath.getURIString(), new PicardHtsPath(expectedURI).getURIString());
    }
}