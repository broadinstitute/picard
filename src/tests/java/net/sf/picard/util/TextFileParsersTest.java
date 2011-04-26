/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package net.sf.picard.util;

import net.sf.picard.io.IoUtil;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;
import org.testng.Assert;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;

import net.sf.picard.PicardException;

public class TextFileParsersTest {

    private static final String testFile1 = "testdata/net/sf/picard/util/whitepace_text_file.txt";
    private static final String testFile2 = "testdata/net/sf/picard/util/all_ones_text_file.txt";
    private static final String testFile3 = "testdata/net/sf/picard/util/no_grouping_file.txt";
    private static final String testFile4 = "testdata/net/sf/picard/util/tabbed_text_file.txt";
    private static final Object[][] testFile1Data = {
        { "Now", "is", "the", "time" },
        { "for", "all", "good", "men" },
        { "to", "come", "to", "the" },
        { "aid", "of", "their", "country." },
        { 15.0d, 23, 55, 67.88888}
    };

    @Test(dataProvider = "basicInputParserData")
    public void testTextFileParser(Object fileOrStream) {
        FormatUtil format = new FormatUtil();
        BasicInputParser parser = fileOrStream instanceof File
            ? new BasicInputParser(true, (File)fileOrStream )
            : new BasicInputParser(true, (InputStream)fileOrStream);
        int index = 0;
        while (parser.hasNext())
        {
            String parts[] = parser.next();
            Assert.assertEquals(parts.length, 4);
            if (index < 4) {
                for (int i = 0; i < parts.length; i++) {
                    Assert.assertEquals(parts, testFile1Data[index]);
                }
            }
            else {
                Assert.assertEquals(testFile1Data[index][0], format.parseDouble(parts[0]));
                Assert.assertEquals(testFile1Data[index][1], format.parseInt(parts[1]));
                Assert.assertEquals(testFile1Data[index][2], format.parseInt(parts[2]));
                Assert.assertEquals(testFile1Data[index][3], format.parseDouble(parts[3]));
            }
            index++;
        }
    }

    @DataProvider(name = "basicInputParserData")
    private Object[][] getBasicInputParserData()
    {
        return new Object[][] {
                {new File(testFile1)},
                {IoUtil.openFileForReading(new File(testFile1))}
        };
    }

    @Test(dataProvider = "noGroupingData")
    public void testTextFileParserNoGrouping(Object fileOrStream) {
        BasicInputParser parser = fileOrStream instanceof File
            ? new BasicInputParser(true, (File)fileOrStream)
            : new BasicInputParser(true, (InputStream)fileOrStream);
        parser.setTreatGroupedDelimitersAsOne(false);
        while (parser.hasNext()) {
            String parts[] = parser.next();
            for (int i = 0; i < parts.length; i++) {
                if (parts[i] != null) {
                    Assert.assertEquals(Integer.parseInt(parts[i]), i+1);
                }
            }
        }
    }

    @DataProvider(name = "noGroupingData")
    private Object[][] getNoGroupingData()
    {
        return new Object[][] {
                {new File(testFile3)},
                {IoUtil.openFileForReading(new File(testFile3))}
        };
    }


    @Test(dataProvider = "leadingWhiteSpaceData")
    public void testTextFileParserLeadingWhitespace(Object fileOrStream) {
        BasicInputParser parser = fileOrStream instanceof File
            ? new BasicInputParser(true, (File)fileOrStream)
            : new BasicInputParser(true, (InputStream)fileOrStream);
        while (parser.hasNext())
        {
            String parts[] = parser.next();
            Assert.assertEquals(parts.length, 1);
            Assert.assertEquals("1", parts[0]);
        }
    }

    @DataProvider(name = "leadingWhiteSpaceData")
    private Object[][] getLeadingWhiteSpaceData()
    {
        return new Object[][] {
                {new File(testFile2)},
                {IoUtil.openFileForReading(new File(testFile2))}
        };
    }


    @Test(expectedExceptions= PicardException.class, dataProvider = "tooManyWordsData")
    public void testTooManyWords(Object fileOrStream) {
        BasicInputParser parser = fileOrStream instanceof File
            ? new BasicInputParser(true, 3, (File)fileOrStream)
            : new BasicInputParser(true, 3, (InputStream)fileOrStream);
        if (parser.hasNext()) {
            String parts[] = parser.next();
        }
        Assert.fail("Attempt to parse extra-long file should have failed but didn't.");
    }

    @DataProvider(name = "tooManyWordsData")
    private Object[][] getTooManyWordsData()
    {
        return new Object[][] {
                {new File(testFile1)},
                {IoUtil.openFileForReading(new File(testFile1))}
        };
    }

    @Test(dataProvider = "tabbedData")
    public void testTabbedFileParser(Object fileOrStream) {
        TabbedInputParser parser = fileOrStream instanceof File
            ? new TabbedInputParser(false, (File)fileOrStream)
            : new TabbedInputParser(false, (InputStream)fileOrStream);
        while (parser.hasNext()) {
            String parts[] = parser.next();
            for (int i = 0; i < parts.length; i++) {
                if (parts[i] != null && !parts[i].equals("")) {
                    Assert.assertEquals(parts[i].trim(), String.valueOf(i+1));
                }
            }
        }
    }

    @DataProvider(name = "tabbedData")
    private Object[][] getTabbedData()
    {
        return new Object[][] {
                {new File(testFile4)},
                {IoUtil.openFileForReading(new File(testFile4))}
        };
    }

    @Test(dataProvider="data")
    public void testWordCountCalculation(String line, boolean groupDelimiters, String name) {

        WordCountTestParser parser = new WordCountTestParser();
        parser.setDelimiter("\t ");
        parser.setTreatGroupedDelimitersAsOne(groupDelimiters);
        parser.calculateWordCount(line.getBytes());
        Assert.assertEquals(parser.getWordCount(), 3, name);
    }

    @DataProvider(name = "data")
    private Object[][] getWordCountCalculationData()
    {
        return new Object[][]{
                {"1\t2\t3", false, "Tabs with all fields filled."},
                {"1\t2\t", false, "Tabs with no final field."},
                {"\t2\t3", false, "Tabs with no first field."},
                {"\t2\t", false, "Tabs with no first or final field."},
                {"1  2  3", true, "Spaces with all fields filled  (grouping on)."},
                {"1  2  3  ", true, "Spaces with no final field (grouping on)."},
                {"   2   3   4", true, "Spaces with no first field (grouping on)."},
                {" 2 ", false, "Spaces with no first or final field."}
        };
    }

    /**
     * Toy class for testing the word count functionality
     */
    private static class WordCountTestParser extends AbstractInputParser {

        private char delimiters[] = null;
        
        public WordCountTestParser() {
        }

        public void setDelimiter(String delim) {
            delimiters = delim.toCharArray();
        }

        protected boolean isDelimiter(final byte b) {
            for (int i = 0; i < delimiters.length; i++) {
                if (b == delimiters[i]) {
                    return true;
                }
            }
            return false;
        }

        protected byte[] readNextLine() {  return new byte[0]; }
        public String getFileName() { return null; }
        public void close() {}
    }
}
