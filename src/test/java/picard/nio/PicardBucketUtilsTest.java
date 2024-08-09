package picard.nio;

import htsjdk.io.HtsPath;
import htsjdk.io.IOPath;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.util.GCloudTestUtils;

import java.net.URI;
import java.nio.file.Path;
import java.nio.file.Paths;

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
    public void testGetTempFilePath(final IOPath directory, final String prefix, final String extension){
        final IOPath path = PicardBucketUtils.getTempFilePath(directory, prefix, extension);
        Assert.assertTrue(path.hasExtension(extension));
        if (directory != null){
            Assert.assertTrue(path.getURIString().startsWith(directory.getURIString()));
        } else {
            Assert.assertEquals(path.getScheme(), PicardBucketUtils.FILE_SCHEME);
        }
    }

    @DataProvider
    public Object[][] getVariousPathsForPrefetching(){
        return new Object[][]{
                {new HtsPath("file:///local/file"), false},
                {new HtsPath("gs://abucket/bucket"), true},
                {new HtsPath("gs://abucket_with_underscores"), true},
        };
    }

    @Test(groups="bucket", dataProvider = "getVariousPathsForPrefetching")
    public void testIsEligibleForPrefetching(final IOPath path, boolean isPrefetchable){
        Assert.assertEquals(PicardBucketUtils.isEligibleForPrefetching(path), isPrefetchable);
    }

}