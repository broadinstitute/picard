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
package net.sf.samtools.util;

import net.sf.samtools.*;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class IupacTest {
    @Test(dataProvider = "basicDataProvider")
    public void basic(final String tempFileExtension) throws Exception {
        final File outputFile = File.createTempFile("iupacTest.", tempFileExtension);
        outputFile.deleteOnExit();
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(new SAMFileHeader(), false, outputFile);
        final String bases1 = "=ACMGRSVTWYHKDBNA";
        final String bases2 = "A=ACMGRSVTWYHKDBN"; // Test both high-order and low-order base encoding.
        final byte[] quals = new byte[bases1.length()];
        Arrays.fill(quals, (byte)20);
        final String[] reads = {bases1,  bases1.toLowerCase(), bases2, bases2.toLowerCase()};
        for (int i = 0; i < reads.length; ++i) {
            final SAMRecord rec = new SAMRecord(writer.getFileHeader());
            rec.setReadName("read" + i);
            rec.setReadUnmappedFlag(true);
            rec.setReadString(reads[i]);
            rec.setBaseQualities(quals);
            writer.addAlignment(rec);
        }
        writer.close();
        final SAMFileReader reader = new SAMFileReader(outputFile);
        final SAMRecordIterator it = reader.iterator();
        for (int i = 0; i < reads.length; ++i) {
            final SAMRecord rec = it.next();
            Assert.assertEquals(rec.getReadString(), reads[i].toUpperCase());
        }
        reader.close();
    }

    @DataProvider(name = "basidDataProvider")
    public Object[][] basicDataProvider() {
        return new Object[][] {
                {".bam"},
                {".sam"}
        };
    }
}
