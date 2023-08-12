package picard.nio;

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.util.GCloudTestUtils;

import static org.testng.Assert.*;

public class GATKBucketUtilsTest {

    // Check that the extension scheme (e.g. "txt" vs ".txt" as the second argument) is consistent for cloud and local files
    @Test
    public void testExtensionConsistent(){
        final String cloudPath = GATKBucketUtils.getTempFilePath(GCloudTestUtils.getTestInputPath(), "txt");
        final String localPath = GATKBucketUtils.getTempFilePath("test", "txt");

        int d = 3;
        Assert.assertTrue(cloudPath.endsWith(".txt"));
        Assert.assertTrue(localPath.endsWith(".txt"));

    }

}