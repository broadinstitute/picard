package org.broad.tribble.util;


import org.testng.Assert;
import org.testng.annotations.Test;


/**
 * Parsing utils tests
 */
public class ParsingUtilsTest {

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
        
}
