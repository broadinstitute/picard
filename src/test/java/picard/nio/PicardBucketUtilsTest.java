package picard.nio;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.util.GCloudTestUtils;

public class PicardBucketUtilsTest {

    @DataProvider(name="testGetTempFilePathDataProvider")
    public Object[][] testGetTempFilePathDataProvider() {
        return new Object[][] {
                {GCloudTestUtils.getTestInputPath(), "", ".txt"},
                {GCloudTestUtils.getTestInputPath(),  "test", ".txt"},
                {null, "", ".log"},
                {null, "test", ".log"}
        };
    }

    // Check that the extension scheme is consistent for cloud and local files
    @Test(dataProvider = "testGetTempFilePathDataProvider", groups = "cloud")
    public void testGetTempFilePath(final String directory, final String prefix, final String extension){
        PicardHtsPath path = PicardBucketUtils.getTempFilePath(directory, prefix, extension);
        Assert.assertTrue(path.hasExtension(extension));
        if (directory != null){
            Assert.assertTrue(path.getURIString().startsWith(directory));
        }

        if (directory == null){
            Assert.assertEquals(path.getScheme(), "file");
        }
    }

}