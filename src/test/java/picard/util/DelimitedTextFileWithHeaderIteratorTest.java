package picard.util;

import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.StringUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class DelimitedTextFileWithHeaderIteratorTest {
    @Test
    public void basicParsingTest() throws Exception {
        final String[][] data = new String[][] {
                new String[] {"FOO", "BAR", "SPLAT"},
                new String[] {"1", "2", "3"},
                new String[] {"a", "b", "c"},
                new String[] {"foo", "bar", "splat"},
        };

        final File tmp = File.createTempFile("tabbedTextTest.", ".txt");
        tmp.deleteOnExit();
        final BufferedWriter out = new BufferedWriter(new FileWriter(tmp));

        for (final String[] fields : data) {
            out.write(StringUtil.join("\t", fields));
            out.newLine();
        }

        out.close();

        final TabbedInputParser parser = new TabbedInputParser(false, tmp);
        final DelimitedTextFileWithHeaderIterator fileIterator = new DelimitedTextFileWithHeaderIterator(parser);
        for (final String col : data[0]) Assert.assertTrue(fileIterator.hasColumn(col));

        int i=1;
        for (final DelimitedTextFileWithHeaderIterator.Row row : new IterableAdapter<DelimitedTextFileWithHeaderIterator.Row>(fileIterator)) {
            final String[] expected = data[i++];
            Assert.assertEquals(row.getFields(), expected);
            Assert.assertEquals(row.getCurrentLine(), StringUtil.join("\t", expected));
        }

        Assert.assertEquals(i, data.length);
    }

    @Test
    public void parsingWithColumnHeadersTest() throws Exception {
        final String[] headers = {"STRING", "STRING2", "NUMBER"};

        final String[][] data = new String[][] {
                headers,
                new String[] {"1", "2", "3"},
                new String[] {"a", "b", "2"},
                new String[] {"foo", "bar", ""},
        };


        final File tmp = File.createTempFile("tabbedTextTest.", ".txt");
        tmp.deleteOnExit();
        final BufferedWriter out = new BufferedWriter(new FileWriter(tmp));

        for (final String[] fields : data) {
            out.write(StringUtil.join("\t", fields));
            out.newLine();
        }

        out.close();

        final TabbedInputParser parser = new TabbedInputParser(false, tmp);
        final DelimitedTextFileWithHeaderIterator fileIterator = new DelimitedTextFileWithHeaderIterator(parser);
        for (final String col : headers) Assert.assertTrue(fileIterator.hasColumn(col));

        int i=1;
        for (final DelimitedTextFileWithHeaderIterator.Row row : new IterableAdapter<DelimitedTextFileWithHeaderIterator.Row>(fileIterator)) {
            final String[] expected = data[i++];
            final String[] actual = row.getFields();
            for (int j = 0; j < expected.length; j++) {
                Assert.assertTrue((expected[j].equals("") && actual[j] == null) || expected[j].equals(actual[j]));
            }
            Assert.assertEquals(row.getCurrentLine(), StringUtil.join("\t", expected));
            try {
                row.getField(headers[0]);
                row.getField(headers[1]);
                row.getIntegerField(headers[2]);
            }
            catch(Exception e) {
                Assert.fail("Failed to parse one of the fields in " + row.getCurrentLine() + ": " + e.getMessage());
            }
        }

        Assert.assertEquals(i, data.length);
    }


}
