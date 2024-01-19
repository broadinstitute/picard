package picard.nio;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

import java.io.File;

public class PicardHtsPathUnitTest {
    private static final boolean APPEND = true;
    private static final boolean REPLACE = false;


    @DataProvider(name="ReplaceExtensionTestInput")
    public Object[][] getReplaceExtensionTestInput() {
        return new Object[][] {
                {"gs://hellbender/test/resources/hg19mini.fasta", ".fai", REPLACE, "gs://hellbender/test/resources/hg19mini.fai"},
                {"gs://hellbender/test/resources/hg19mini.fasta", ".fai", APPEND, "gs://hellbender/test/resources/hg19mini.fasta.fai"},
                {"testdata/picard/sam/test.bam", ".bai", REPLACE, new File("testdata/picard/sam/test.bai").getAbsolutePath()}, // a relative input will be converted to absolute input
                {"testdata/picard/sam/test.bam", ".bai", APPEND, new File("testdata/picard/sam/test.bam.bai").getAbsolutePath()},
                {"/Users/jdoe/workspace/picard/testdata/picard/sam/test.bam", ".bai", REPLACE, "/Users/jdoe/workspace/picard/testdata/picard/sam/test.bai"},
                {"/Users/jdoe/workspace/picard/testdata/picard/sam/test.bam", ".md5", APPEND, "/Users/jdoe/workspace/picard/testdata/picard/sam/test.bam.md5"},
                {"/Users/jdoe/workspace/picard/testdata/picard/sam/my.fasta.gz", ".fai", REPLACE, "/Users/jdoe/workspace/picard/testdata/picard/sam/my.fasta.fai"},
                {"/Users/jdoe/workspace/picard/testdata/picard/sam/my.fasta.gz", ".md5", APPEND, "/Users/jdoe/workspace/picard/testdata/picard/sam/my.fasta.gz.md5"},
                {"gs://hellbender/test/resources.fasta/hg19mini.fasta", ".fai", REPLACE, "gs://hellbender/test/resources.fasta/hg19mini.fai"}, // "fasta" is part of the path
                {"gs://hellbender/test/resources.fasta/hg19mini.fasta", ".fai", APPEND, "gs://hellbender/test/resources.fasta/hg19mini.fasta.fai"},
                {"gs://hellbender/test/resources/hg19mini", ".fai", APPEND, "gs://hellbender/test/resources/hg19mini.fai"}, // no extension
        };
    }

    @Test(groups = "cloud", dataProvider = "ReplaceExtensionTestInput")
    public void testReplaceExtension(final String originalURI, final String newExtension, final boolean append,
                                     final String expectedURI){
        final PicardHtsPath originalPath = new PicardHtsPath(originalURI);
        final PicardHtsPath newPath = PicardHtsPath.replaceExtension(originalPath, newExtension, append);
        // We cannot directly compare PicardHtsPaths because replaceExtension() creates a new PicardHtsPath with an absolute path as the rawInputString,
        // even if the input has a rawInputString that is a relative path. Instead, we check that the absolute URIs match.
        Assert.assertEquals(newPath.getURIString(), new PicardHtsPath(expectedURI).getURIString());
    }

    @DataProvider(name="ReplaceExtensionFailureTestInput")
    public Object[][] getReplaceExtensionFailureTestInput() {
        return new Object[][]{
                {"gs://hellbender/test/resources/hg19mini", ".fai", REPLACE}
        };
    }

    @Test(groups = "cloud", dataProvider = "ReplaceExtensionFailureTestInput", expectedExceptions = PicardException.class)
    public void testReplaceExtensionFailure(final String originalURI, final String newExtension, final boolean append){
        final PicardHtsPath originalPath = new PicardHtsPath(originalURI);
        final PicardHtsPath newPath = PicardHtsPath.replaceExtension(originalPath, newExtension, append);
    }
}