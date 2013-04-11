package org.broad.tribble;

import org.broad.tribble.bed.BEDCodec;
import org.testng.annotations.Test;

import static org.testng.Assert.assertNotNull;

/**
 * @author jacob
 * @date 2013-Apr-10
 */
public class AbstractFeatureReaderTest {

    @Test(enabled = false, description = "Tribble doesn't have full FTP support yet")
    public void testLoadBEDFTP() throws Exception {
        String path = "ftp://ftp.broadinstitute.org/distribution/igv/TEST/cpgIslands with spaces.hg18.bed";
        FeatureCodec codec = new BEDCodec();
        AbstractFeatureReader<Feature> bfs = AbstractFeatureReader.getFeatureReader(path, codec, false);
        for(Feature feat: bfs.iterator()){
            assertNotNull(feat);
        }

    }

}
