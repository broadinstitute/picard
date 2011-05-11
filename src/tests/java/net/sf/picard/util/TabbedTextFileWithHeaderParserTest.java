package net.sf.picard.util;

import net.sf.picard.io.IoUtil;
import net.sf.picard.util.TabbedTextFileWithHeaderParser.Row;
import net.sf.samtools.util.StringUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

/**
 *
 */
public class TabbedTextFileWithHeaderParserTest {
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

        final TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(tmp);
        for (final String col : data[0]) Assert.assertTrue(parser.hasColumn(col));

        int i=1;
        for (final Row row : parser) {
            final String[] expected = data[i++];
            Assert.assertEquals(row.getFields(), expected);
            Assert.assertEquals(row.getCurrentLine(), StringUtil.join("\t", expected));
        }

        Assert.assertEquals(i, data.length);
    }
}
