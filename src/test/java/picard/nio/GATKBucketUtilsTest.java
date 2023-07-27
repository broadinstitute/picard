package picard.nio;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import org.testng.Assert;

import java.net.URI;
import java.nio.file.Path;
import java.nio.file.Paths;

public class GATKBucketUtilsTest {

    @DataProvider
    public Object[][] getVariousPathsForPrefetching(){
        return new Object[][]{
                {"file:///local/file", false},
                {"gs://abucket/bucket", true},
                {"gs://abucket_with_underscores", true},
        };
    }

    @Test(groups="bucket", dataProvider = "getVariousPathsForPrefetching")
    public void testIsEligibleForPrefetching(String path, boolean isPrefetchable){
        final URI uri = URI.create(path);
        final Path uriPath = Paths.get(uri);
        Assert.assertEquals(GATKBucketUtils.isEligibleForPrefetching(uriPath), isPrefetchable);
    }
}
