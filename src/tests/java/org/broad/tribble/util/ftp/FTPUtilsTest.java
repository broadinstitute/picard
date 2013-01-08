package org.broad.tribble.util.ftp;

import net.sf.samtools.util.ftp.FTPUtils;
import org.testng.annotations.Test;

import java.net.URL;

import static org.testng.Assert.assertFalse;
import static org.testng.Assert.assertTrue;

/**
* @author Jim Robinson
* @since 10/4/11
*/
public class FTPUtilsTest {

    @Test
    public void testResourceAvailable() throws Exception {

        URL goodUrl = new URL("ftp://ftp.broadinstitute.org/pub/igv/TEST/test.txt");
        assertTrue(FTPUtils.resourceAvailable(goodUrl));

        URL nonExistentURL = new URL("ftp://ftp.broadinstitute.org/pub/igv/TEST/doesntExist");
        assertFalse(FTPUtils.resourceAvailable(nonExistentURL));

        URL nonExistentServer = new URL("ftp://noSuchServer/pub/igv/TEST/doesntExist");
        assertFalse(FTPUtils.resourceAvailable(nonExistentServer));


    }
}
