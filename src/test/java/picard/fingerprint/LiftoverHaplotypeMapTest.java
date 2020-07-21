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
import htsjdk.samtools.util.CloserUtil;
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
public class LiftoverHaplotypeMapTest {

    private static final String TEST_DATA_DIR="testdata/picard/fingerprint/";
    private static final String HG18_HAPLOTYPE_MAP = TEST_DATA_DIR + "haplotypeMapHg18.txt";
    private static final String CHAIN = TEST_DATA_DIR + "hg18ToHg19.broad.over.chain";
    private static final String HG_19_SD = TEST_DATA_DIR + "Homo_sapiens_assembly19.dict";

    @Test
    public void testHaplotypeMapLiftover() throws Exception {

        Map<String, Interval> hg19 = new HashMap<String, Interval>();
        hg19.put("rs1210110", new Interval("1", 14096821, 14096821));
        hg19.put("rs7555566", new Interval("1", 14804874, 14804874));
        hg19.put("rs1364054", new Interval("2", 8178735, 8178735));
        hg19.put("rs6734275", new Interval("2", 67241174, 67241174));
        hg19.put("rs7584993", new Interval("2", 223845942, 223845942));
        hg19.put("rs17272796", new Interval("3", 17077268, 17077268));
        hg19.put("rs1155741", new Interval("3", 37627112, 37627112));
        hg19.put("rs1155744", new Interval("3", 37627598, 37627598));
        hg19.put("rs161792", new Interval("3", 151899704, 151899704));
        hg19.put("rs11940551", new Interval("4", 27162478, 27162478));
        hg19.put("rs9293511", new Interval("5", 88416354, 88416354));
        hg19.put("rs10943567", new Interval("6", 79402451, 79402451));
        hg19.put("rs9343786", new Interval("6", 79414728, 79414728));
        hg19.put("rs9352613", new Interval("6", 79424433, 79424433));
        hg19.put("rs685449", new Interval("6", 153344531, 153344531));
        hg19.put("rs7808249", new Interval("7", 86983715, 86983715));
        hg19.put("rs1106334", new Interval("8", 71012811, 71012811));
        hg19.put("rs11017876", new Interval("10", 129200964, 129200964));
        hg19.put("rs9572094", new Interval("13", 35252882, 35252882));
        hg19.put("rs4905366", new Interval("14", 96103099, 96103099));
        hg19.put("rs4775699", new Interval("15", 47873549, 47873549));
        hg19.put("rs1528601", new Interval("16", 51098427, 51098427));
        hg19.put("rs11655512", new Interval("17", 20851735, 20851735));
        hg19.put("rs4793172", new Interval("17", 43131480, 43131480));
        hg19.put("rs242076", new Interval("22", 33229830, 33229830));
        hg19.put("rs6603251", new Interval("X", 320580, 320580));


        final File tmp = File.createTempFile("hg19Map", ".txt");
        tmp.deleteOnExit();

        new LiftOverHaplotypeMap().instanceMain(new String[] {"INPUT=" + HG18_HAPLOTYPE_MAP,
                "OUTPUT=" + tmp.getAbsolutePath(), "SEQUENCE_DICTIONARY=" + HG_19_SD,
                "CHAIN=" + CHAIN });
        final HaplotypeMap results = new HaplotypeMap(tmp);
        final SamReader samReader = SamReaderFactory.makeDefault().open(new File(HG_19_SD));
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
        CloserUtil.close(samReader);
    }
}
