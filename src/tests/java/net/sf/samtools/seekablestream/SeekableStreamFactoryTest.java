package net.sf.samtools.seekablestream;

import org.testng.Assert;
import org.testng.annotations.Test;

public class SeekableStreamFactoryTest {
    @Test
    public void testIsFilePath() throws Exception {
        Assert.assertEquals(SeekableStreamFactory.isFilePath("x"), true);
        Assert.assertEquals(SeekableStreamFactory.isFilePath(""), true);
        Assert.assertEquals(SeekableStreamFactory.isFilePath("http://broadinstitute.org"), false);
        Assert.assertEquals(SeekableStreamFactory.isFilePath("https://broadinstitute.org"), false);
        Assert.assertEquals(SeekableStreamFactory.isFilePath("ftp://broadinstitute.org"), false);
    }
}
