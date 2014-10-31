/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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

package picard.illumina;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

import static picard.illumina.CollectIlluminaSummaryMetrics.*;

/**
 * @author farjoun@broadinstitute.org
 */

public class CollectIlluminaSummaryMetricsTest {
    public static final String TEST_DATA_DIR = "testdata/picard/quality/";
    public static final String CHR_M_REFERENCE = TEST_DATA_DIR+"chrM.reference.fasta";

    private File rootTestDir;

    @BeforeTest
    private void setUp() throws Exception {
        rootTestDir = File.createTempFile("CollectIlluminaSummaryMetricsTest.", ".tmp");
        Assert.assertTrue(rootTestDir.delete());
        Assert.assertTrue(rootTestDir.mkdir());
    }

    @AfterTest
    private void tearDown() {
        IOUtil.deleteDirectoryTree(rootTestDir);
    }

    @Test(expectedExceptions = PicardException.class)
    public void testIlluminaSummaryReadLengthFail() throws Exception {
        getSummaryFile("null",TEST_DATA_DIR+"chrMReadsDiffereingLengths.sam", CHR_M_REFERENCE,"/dev/null");
    }

    @DataProvider(name="IlluminaSummaryData")
    public Object[][] testIlluminaSummaryData() {

        return new Object[][]{
                {"",TEST_DATA_DIR+"chrMReads.sam",                       5, 2, 2, 5*36L,   2*36L, 0L,   0D,            (5*36L)/5, 16299L}, // only dups and failing reads
                {"",TEST_DATA_DIR+"chrMReadsMated.sam",                  2, 2, 0, 2*36L,   2*36L,66L,   66.0D/16299,   36L,         16299L}, // nice paired reads
                {"",TEST_DATA_DIR+"chrMReadsMatedRev.sam",               2, 2, 0, 2*36L,   2*36L,66L,   66.0D/16299,   36L,         16299L}, // nice paired reads (reversed order)
                {"",TEST_DATA_DIR+"chrMReadsMatedSameStrand.sam",        2, 2, 0, 2*36L,   2*36L,46L,   46.0D/16299,   36L,         16299L}, // nice paired reads (same strand)
                {"",TEST_DATA_DIR+"chrMReadsMateUnmapped.sam",           2, 2, 0, 2*36L,   2*36L,2*36L, 2*36.0D/16299, 36L,         16299L}, // one mate is unmapped
                {"",TEST_DATA_DIR+"chrMReadsUnMated.sam",                1, 1, 0, 1*36L,   1*36L,1*36L, 1*36.0D/16299, 36L,         16299L}, // one mate is unmapped
                {"",TEST_DATA_DIR+"chrMReadsSecondaryAlignSoft.sam",     2, 2, 0, 2*36L,   2*36L,66L,   66.0D/16299,   36L,         16299L}, // secondary read with soft clips
                {"",TEST_DATA_DIR+"chrMReadsSupplementalAlignSoft.sam",  2, 2, 0, 2*36L,   2*36L,66L,   66.0D/16299,   36L,         16299L}, // secondary read with soft clips
                {TEST_DATA_DIR+"chrM.part.intervallist",
                TEST_DATA_DIR+"chrMReadsMated.sam",                      2, 2, 0, 2*36L,   2*36L,66L,   66.0D/2000,    36L,         2000L},  // nice paired reads  with an interval list
        };
    }

    @Test(dataProvider = "IlluminaSummaryData")
    public void testIlluminaSummary(final String intervalFile, final String fileBAM, final Integer TOTAL_READS, final Integer PF_READS, final Integer DUPLICATE_READS,
                                    final Long TOTAL_BASES, final Long PF_BASES, final Long TOTAL_ILMN_BASES, final Double AVERAGE_ILMN_DEPTH, final Long READ_LENGTH, final Long TARGET_REGION) throws Exception {

        final IlluminaSummaryMetrics expectedMetric=new IlluminaSummaryMetrics( TOTAL_READS,  PF_READS,  DUPLICATE_READS, TOTAL_BASES,  PF_BASES,  TOTAL_ILMN_BASES,  AVERAGE_ILMN_DEPTH,  READ_LENGTH, TARGET_REGION);

            final MetricsFile<IlluminaSummaryMetrics, ?> metricsFile = getSummaryFile(intervalFile,fileBAM, CHR_M_REFERENCE, rootTestDir + "/READ_TEST");
        final IlluminaSummaryMetrics metrics = metricsFile.getMetrics().get(0);

        Assert.assertEquals(metrics.toString(), expectedMetric.toString());
    }

    private MetricsFile<IlluminaSummaryMetrics, ?> getSummaryFile(final String intervalFile, final String input, final String reference, final String outputFile) throws Exception {
        final List<String> argList = new ArrayList<String>();
        argList.add("I=" + input);
        argList.add("O=" + outputFile);
        argList.add("R=" + reference);
        if(!intervalFile.equals("")) argList.add("TARGET_REGION="+intervalFile);

        final String[] args = new String[argList.size()];
        argList.toArray(args);

        Assert.assertEquals(new CollectIlluminaSummaryMetrics().instanceMain(args), 0);

        final MetricsFile<IlluminaSummaryMetrics, Integer> retVal = new MetricsFile<IlluminaSummaryMetrics, Integer>();
        retVal.read(new FileReader(outputFile));
        return retVal;
    }


}
