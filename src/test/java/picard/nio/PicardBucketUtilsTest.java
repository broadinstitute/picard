package picard.nio;

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.util.GCloudTestUtils;

public class PicardBucketUtilsTest {

    // Check that the extension scheme (e.g. "txt" vs ".txt" as the second argument) is consistent for cloud and local files
    @Test
    public void testExtensionConsistent(){
        final PicardHtsPath cloudPath = PicardBucketUtils.getTempFilePath(GCloudTestUtils.getTestInputPath(), ".txt");
        final PicardHtsPath localPath = PicardBucketUtils.getTempFilePath("test", ".txt");

        Assert.assertTrue(cloudPath.hasExtension(".txt"));
        Assert.assertTrue(localPath.hasExtension(".txt"));
    }

}