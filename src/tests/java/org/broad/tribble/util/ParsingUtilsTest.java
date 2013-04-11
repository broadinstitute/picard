package org.broad.tribble.util;


import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.IOException;
import java.io.InputStream;


/**
 * Parsing utils tests
 */
public class ParsingUtilsTest {

    static final String AVAILABLE_FTP_URL = "ftp://ftp.broadinstitute.org/pub/igv/TEST/test.txt";
    static final String UNAVAILABLE_FTP_URL = "ftp://www.example.com/file.txt";

    static final String AVAILABLE_HTTP_URL = "http://www.google.com";
    static final String UNAVAILABLE_HTTP_URL = "http://www.unknownhostwhichshouldntexist.com";

    @Test
    public void testSplit1() {
        String[] tokens = new String[10];
        String blankColumnLine = "a\tb\t\td";
        int nTokens = ParsingUtils.split(blankColumnLine, tokens, '\t');
        Assert.assertEquals(nTokens,4);
        Assert.assertEquals(tokens[0],"a");
        Assert.assertEquals(tokens[1],"b");
        Assert.assertEquals(tokens[2],"");
        Assert.assertEquals(tokens[3],"d");
    }

    @Test
    public void testSplit2() {
        String[] tokens = new String[10];
        String blankColumnLine = "a\tb\t\td\t";
        int nTokens = ParsingUtils.split(blankColumnLine, tokens, '\t');
        Assert.assertEquals(nTokens,5);
        Assert.assertEquals(tokens[0],"a");
        Assert.assertEquals(tokens[1],"b");
        Assert.assertEquals(tokens[2],"");
        Assert.assertEquals(tokens[3],"d");
        Assert.assertEquals(tokens[4],"");
    }

    @Test
    public void testSplitWhitespace1() {
        String[] tokens = new String[10];
        String blankColumnLine = "a b\t\td";
        int nTokens = ParsingUtils.splitWhitespace(blankColumnLine, tokens);
        Assert.assertEquals(nTokens,4);
        Assert.assertEquals(tokens[0],"a");
        Assert.assertEquals(tokens[1],"b");
        Assert.assertEquals(tokens[2],"");
        Assert.assertEquals(tokens[3],"d");
    }

    @Test
    public void testSplitWhitespace2() {
        String[] tokens = new String[10];
        String blankColumnLine = "a b\t\td\t";
        int nTokens = ParsingUtils.splitWhitespace(blankColumnLine, tokens);
        Assert.assertEquals(nTokens,5);
        Assert.assertEquals(tokens[0],"a");
        Assert.assertEquals(tokens[1],"b");
        Assert.assertEquals(tokens[2],"");
        Assert.assertEquals(tokens[3],"d");
    }

    @Test
    public void testFTPDoesExist() throws IOException{
        tstExists(AVAILABLE_FTP_URL, true);
    }

    @Test
    public void testFTPNotExist() throws IOException{
        tstExists(UNAVAILABLE_FTP_URL, false);
    }

    @Test
    public void testHTTPDoesExist() throws IOException{
        tstExists(AVAILABLE_HTTP_URL, true);
    }

    @Test
    public void testHTTPNotExist() throws IOException{
        tstExists(UNAVAILABLE_HTTP_URL, false);
    }

    private void tstExists(String path, boolean expectExists) throws IOException{
        boolean exists = ParsingUtils.resourceExists(path);
        Assert.assertEquals(exists, expectExists);
    }

    @Test
    public void testFTPOpenInputStream() throws IOException{
        tstStream(AVAILABLE_FTP_URL);
    }

    @Test
    public void testHTTPOpenInputStream() throws IOException{
        tstStream(AVAILABLE_HTTP_URL);
    }

    private void tstStream(String path) throws IOException{
        InputStream is = ParsingUtils.openInputStream(path);
        Assert.assertNotNull(is, "InputStream is null for " + path);
        int b = is.read();
        Assert.assertNotSame(b, -1);
    }


}
