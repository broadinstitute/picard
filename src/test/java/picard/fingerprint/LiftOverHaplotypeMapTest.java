/*
 * The MIT License
 *
 * Copyright (c) 2020 The Broad Institute
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
package picard.fingerprint;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SequenceUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

/**
 * Basic test for HaplotypeMap liftover
 */
public class LiftOverHaplotypeMapTest {

    private static final String TEST_DATA_DIR = "testdata/picard/fingerprint/";
    private static final String HG18_HAPLOTYPE_MAP = TEST_DATA_DIR + "haplotypeMapHg18.txt";
    private static final String CHAIN = TEST_DATA_DIR + "hg18ToHg19.clipped.over.chain";
    private static final String HG_19_SD = TEST_DATA_DIR + "Homo_sapiens_assembly19.dict";

    @Test
    public void testHaplotypeMapLiftover() throws Exception {

        Map<String, Interval> hg19 = new HashMap<>();
        hg19.put("rs1210110", new Interval("1", 14096821, 14096821));
        hg19.put("rs7555566", new Interval("1", 14804874, 14804874));
        hg19.put("rs1364054", new Interval("2", 8178735, 8178735));
        hg19.put("rs6734275", new Interval("2", 67241174, 67241174));
        hg19.put("rs7584993", new Interval("2", 223845942, 223845942));

        final File tmp = File.createTempFile("hg19Map", ".txt");
        tmp.deleteOnExit();

        new LiftOverHaplotypeMap().instanceMain(new String[]{"INPUT=" + HG18_HAPLOTYPE_MAP,
                "OUTPUT=" + tmp.getAbsolutePath(), "SEQUENCE_DICTIONARY=" + HG_19_SD,
                "CHAIN=" + CHAIN});
        final HaplotypeMap results = new HaplotypeMap(tmp);
        try(final SamReader samReader = SamReaderFactory.makeDefault().open(new File(HG_19_SD))) {
            SequenceUtil.assertSequenceDictionariesEqual(
                    results.getHeader().getSequenceDictionary(),
                    samReader.getFileHeader().getSequenceDictionary());

            for (Snp snp : results.getAllSnps()) {
                Interval i = hg19.remove(snp.getName());
                Assert.assertNotNull(i, "No lifted-over interval found for Snp: " + snp.getName());
                Assert.assertEquals(snp.getChrom(), i.getContig(), "Sequence doesn't match for " + snp.getName());
                Assert.assertEquals(snp.getPos(), i.getStart(), "Start position doesn't match for " + snp.getName());
            }
            Assert.assertEquals(hg19.size(), 0, "Extra intervals remain after liftover: " + hg19);
        }
    }
}
