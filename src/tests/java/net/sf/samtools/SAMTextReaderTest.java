/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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
package net.sf.samtools;

import net.sf.samtools.util.CloseableIterator;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;

public class SAMTextReaderTest {
    // Simple input, spot check that parsed correctly, and make sure nothing blows up.
    @Test
    public void testBasic() throws Exception {
        final String seq1 = "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG";
        final String seq2 = "ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA";
        final String qual1 = "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<";
        final String qual2 = "<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<";
        final String fileFormatVersion = "1.0";
        final String sequence = "chr20";
        final int sequenceLength = 62435964;
        final String charTag = "XC";
        final char charValue = 'q';
        final String intTag = "XI";
        final int intValue = 12345;
        final String floatTag = "XF";
        final float floatValue = 1.2345f;
        final String stringTag = "XS";
        final String stringValue = "Hi,Mom!";
        final String samExample = "@HD\tVN:" + fileFormatVersion + "\t" + charTag + ":A:" + charValue + "\n" +
                "@SQ\tSN:" + sequence + "\tAS:HG18\tLN:" + sequenceLength + "\t" + intTag + ":i:" + intValue + "\n" +
                "@RG\tID:L1\tPU:SC_1_10\tLB:SC_1\tSM:NA12891" + "\t" + floatTag + ":f:" + floatValue + "\n" +
                "@RG\tID:L2\tPU:SC_2_12\tLB:SC_2\tSM:NA12891\n" +
                "@PG\tID:0\tVN:1.0\tCL:yo baby\t" + stringTag + ":Z:" + stringValue + "\n" +
                "@PG\tID:2\tVN:1.1\tCL:whassup? ? ? ?\n" +
                "read_28833_29006_6945\t99\tchr20\t28833\t20\t10M1D25M\t=\t28993\t195\t" +
                seq1.toLowerCase() + "\t" + qual1 + "\t" +
                "MF:i:130\tNm:i:1\tH0:i:0\tH1:i:0\tRG:Z:L1\n" +
                "read_28701_28881_323b\t147\tchr20\t28834\t30\t35M\t=\t28701\t-168\t" +
                seq2 + "\t" + qual2 + "\t" +
                "MF:i:18\tNm:i:0\tH0:i:1\tH1:i:0\tRG:Z:L2\n";

        final String[] samResults =
                {"read_28833_29006_6945\t99\tchr20\t28833\t20\t10M1D25M\tchr20\t28993\t195\t" + seq1 + "\t" + qual1 +
                        "\tMF:i:130\tNm:i:1\tH0:i:0\tH1:i:0\tRG:Z:L1",
                "read_28701_28881_323b\t147\tchr20\t28834\t30\t35M\tchr20\t28701\t-168\t" + seq2 + "\t" + qual2 +
                        "\tMF:i:18\tNm:i:0\tH0:i:1\tH1:i:0\tRG:Z:L2"
        };

        final SAMFileReader samReader = createSamFileReader(samExample);
        final SAMFileHeader fileHeader = samReader.getFileHeader();

        Assert.assertEquals(fileHeader.getVersion(), fileFormatVersion);
        Assert.assertEquals(fileHeader.getAttribute(charTag), charValue);
        final SAMSequenceRecord sequenceRecord = fileHeader.getSequence(sequence);
        Assert.assertNotNull(sequenceRecord);
        Assert.assertEquals(sequenceRecord.getSequenceLength(), sequenceLength);
        Assert.assertEquals(sequenceRecord.getAttribute(intTag), intValue);
        Assert.assertEquals(fileHeader.getReadGroup("L1").getAttribute(floatTag), floatValue);
        Assert.assertEquals(fileHeader.getProgramRecord("0").getAttribute(stringTag), stringValue);

        final CloseableIterator<SAMRecord> iterator = samReader.iterator();
        int i = 0;
        while (iterator.hasNext()) {
            final SAMRecord rec = iterator.next();
            Assert.assertEquals(rec.format(), samResults[i++]);
        }
        iterator.close();
        iterator.close();
        samReader.close();        
    }

    private SAMFileReader createSamFileReader(final String samExample) {
        final ByteArrayInputStream inputStream = new ByteArrayInputStream(samExample.getBytes());
        return new SAMFileReader(inputStream);
    }

    @Test
    public void testUnmapped() {
        final String alignmentFromKris =
          "0\t4\t*\t0\t0\t*\t*\t0\t0\tGCCTCGTAGTGCGCCATCAGTCTATCGATGTCGTTG\t44\"44===;;;;;;;;;::::88844\"4\"\"\"\"\"\"\"\"\n";
        final SAMFileReader samReader = createSamFileReader(alignmentFromKris);
        final CloseableIterator<SAMRecord> iterator = samReader.iterator();
        while (iterator.hasNext()) {
            iterator.next();
        }
        iterator.close();

    }

    /**
     * Colon separates fields of a text tag, but colon is also valid in a tag value, so assert that works properly.
     */
    @Test
    public void testTagWithColon() {
        // Create a SAMRecord with a String tag containing a colon
        final SAMRecordSetBuilder samBuilder = new SAMRecordSetBuilder();
        samBuilder.addUnmappedFragment("Hi,Mom!");
        final SAMRecord rec = samBuilder.iterator().next();
        final String valueWithColons = "A:B::C:::";
        rec.setAttribute(SAMTag.CQ.name(),  valueWithColons);
        // Write the record as SAM Text
        final ByteArrayOutputStream os = new ByteArrayOutputStream();
        final SAMFileWriter textWriter = new SAMFileWriterFactory().makeSAMWriter(samBuilder.getHeader(),
                true, os);
        textWriter.addAlignment(rec);
        textWriter.close();

        final SAMFileReader reader = new SAMFileReader(new ByteArrayInputStream(os.toByteArray()));
        final SAMRecord recFromText = reader.iterator().next();
        Assert.assertEquals(recFromText.getAttribute(SAMTag.CQ.name()), valueWithColons);
    }
}
