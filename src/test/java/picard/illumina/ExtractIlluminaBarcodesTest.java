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
package picard.illumina;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.illumina.parser.BaseIlluminaDataProvider;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.IlluminaDataProviderFactory;
import picard.illumina.parser.IlluminaDataType;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.util.BasicInputParser;

import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author alecw@broadinstitute.org
 */
public class ExtractIlluminaBarcodesTest extends CommandLineProgramTest {
    private static final File SINGLE_DATA_DIR = new File("testdata/picard/illumina/25T8B25T/Data/Intensities/BaseCalls");
    private static final File DUAL_DATA_DIR = new File("testdata/picard/illumina/25T8B8B25T/Data/Intensities/BaseCalls");
    private static final String[] BARCODES = {
            "CAACTCTC",
            "CAACTCTG", // This one is artificial -- one edit away from the first one
            "ACAGGTAT",
            "GACCGTTG",
            "ATTATCAA",
            "TGCTGCTG",
            "AACAATGG",
            "TGTAATCA",
            "GCCGTCGA",
            "GTCCACAG",
            "TTGTCTAT",
            "GTGGAGAC",
            "TTGCAAAT"
    };

    private File basecallsDir;
    private File dual;
    private File qual;

    public String getCommandLineProgramName() {
        return ExtractIlluminaBarcodes.class.getSimpleName();
    }

    @BeforeTest
    private void setUp() throws Exception {
        basecallsDir = File.createTempFile("eib.", ".tmp");
        Assert.assertTrue(basecallsDir.delete());
        Assert.assertTrue(basecallsDir.mkdir());
        IOUtil.copyDirectoryTree(SINGLE_DATA_DIR, basecallsDir);
        dual = File.createTempFile("eib_dual", ".tmp");
        Assert.assertTrue(dual.delete());
        Assert.assertTrue(dual.mkdir());
        IOUtil.copyDirectoryTree(DUAL_DATA_DIR, dual);
        qual = File.createTempFile("eib_qual", ".tmp");
        Assert.assertTrue(qual.delete());
        Assert.assertTrue(qual.mkdir());
        IOUtil.copyDirectoryTree(DUAL_DATA_DIR, qual);
    }

    @AfterTest
    private void tearDown() {
        IOUtil.deleteDirectoryTree(basecallsDir);
        IOUtil.deleteDirectoryTree(dual);
        IOUtil.deleteDirectoryTree(qual);
    }

    @Test
    public void testSingleEndWithBarcodeAtStart() throws Exception {
        final MetricsFile<BarcodeMetric, Integer> metricsFile = runIt(1, "8B25T");
        Assert.assertEquals(metricsFile.getMetrics().get(11).PERFECT_MATCHES, 1);
    }

    @Test
    public void testSingleEndWithBarcodeAtStartAndMolecularIndicies() throws Exception {
        final MetricsFile<BarcodeMetric, Integer> metricsFile = runIt(1, "8B4M21T");
        Assert.assertEquals(metricsFile.getMetrics().get(11).PERFECT_MATCHES, 1);
    }

    @Test
    public void testSingleEndWithBarcodeAtEnd() throws Exception {
        final MetricsFile<BarcodeMetric, Integer> metricsFile = runIt(1, "25T8B");
        Assert.assertEquals(metricsFile.getMetrics().get(0).PERFECT_MATCHES, 5);
    }

    @Test
    public void testSingleEndWithBarcodeAtEndAndMolecularIndicies() throws Exception {
        final MetricsFile<BarcodeMetric, Integer> metricsFile = runIt(1, "4M21T8B");
        Assert.assertEquals(metricsFile.getMetrics().get(0).PERFECT_MATCHES, 5);
    }

    @Test
    public void testPairedEndWithBarcodeOnFirstEnd() throws Exception {
        final MetricsFile<BarcodeMetric, Integer> metricsFile = runIt(1, "25T8B25T");
        Assert.assertEquals(metricsFile.getMetrics().get(0).PERFECT_MATCHES, 5);
    }

    @Test
    public void testPairedEndWithBarcodeAndMolecularIndicies() throws Exception {
        final MetricsFile<BarcodeMetric, Integer> metricsFile = runIt(1, "4M21T8B21T4M");
        Assert.assertEquals(metricsFile.getMetrics().get(0).PERFECT_MATCHES, 5);
    }

    @Test
    public void testPairedEndWithBarcodeOnSecondEnd() throws Exception {
        final MetricsFile<BarcodeMetric, Integer> metricsFile = runIt(1, "25T25T8B");
        Assert.assertEquals(metricsFile.getMetrics().get(12).PERFECT_MATCHES, 1);
    }

    @Test
    public void testNonWritableOutputFile() throws Exception {
        final File existingFile = new File(basecallsDir, "s_1_1101_barcode.txt.gz");
        try {
            existingFile.setReadOnly();
            final String readStructure = "25T8B25T";
            final int lane = 1;

            final File metricsFile = File.createTempFile("eib.", ".metrics");
            metricsFile.deleteOnExit();

            final List<String> args = new ArrayList<String>(Arrays.asList(
                    "BASECALLS_DIR=" + basecallsDir.getPath(),
                    "LANE=" + lane,
                    "READ_STRUCTURE=" + readStructure,
                    "METRICS_FILE=" + metricsFile.getPath(),
                    "COMPRESS_OUTPUTS=true"
            ));
            for (final String barcode : BARCODES) {
                args.add("BARCODE=" + barcode);
            }
            Assert.assertEquals(runPicardCommandLine(args), 4);
        } finally {
            existingFile.setWritable(true);
        }

    }

    /**
     * 4 cases tested:
     * * exact match to ACAGTG
     * * inexact match within threshold to TGACCA
     * * inexact match not within threshold to TGACCA
     * * inexact match where the next match is too close to ACAGTG
     *
     * @throws Exception
     */
    @Test
    public void testBarcodeMatching() throws Exception {
        final int lane = 1;
        final int barcodePosition = 26;
        final MetricsFile<BarcodeMetric, Integer> metricsFile = runIt(lane, "25T8B25T");

        BarcodeMetric metricOne = null;
        BarcodeMetric metricTwo = null;
        BarcodeMetric metricNoMatch = null;
        for (final BarcodeMetric metric : metricsFile.getMetrics()) {
            if (metric.BARCODE.equals(BARCODES[0])) {
                metricOne = metric;
            } else if (metric.BARCODE.equals(BARCODES[2])) {
                metricTwo = metric;
            } else if (metric.BARCODE.equals("NNNNNNNN")) {
                metricNoMatch = metric;
            }
        }
        Assert.assertEquals(metricOne.PERFECT_MATCHES, 5);
        Assert.assertEquals(metricOne.ONE_MISMATCH_MATCHES, 0);
        Assert.assertEquals(metricOne.PF_READS, 3);
        Assert.assertEquals(metricOne.READS, 5);

        // one inexact match
        Assert.assertEquals(metricTwo.READS, 4);
        Assert.assertEquals(metricTwo.ONE_MISMATCH_MATCHES, 0);

        Assert.assertEquals(metricNoMatch.READS, 140);
        Assert.assertEquals(metricNoMatch.PF_READS, 112);

        // Check the barcode files themselves
        final File[] barcodeFiles = IOUtil.getFilesMatchingRegexp(basecallsDir, "s_" + lane + "_\\d{4}_barcode.txt");
        Arrays.sort(barcodeFiles);

        final BasicInputParser barcodeParser = new BasicInputParser(true, barcodeFiles);

        // Exact match
        String[] illuminaFields = barcodeParser.next();
        Assert.assertEquals(illuminaFields[1], "Y");
        Assert.assertEquals(illuminaFields[2], "CAACTCTC");

        // Inexact match
        illuminaFields = barcodeParser.next();
        Assert.assertEquals(illuminaFields[1], "Y");
        Assert.assertEquals(illuminaFields[2], "ACAGGTAT");

        // Too many mismatches
        illuminaFields = barcodeParser.next();
        Assert.assertEquals(illuminaFields[1], "N");

        barcodeParser.close();

        // Tack on test of barcode-informed Illumina Basecall parsing
        final ReadStructure rs = new ReadStructure("25T8B25T");
        final IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(basecallsDir, lane, rs,
                new BclQualityEvaluationStrategy(BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY),
                IlluminaDataType.BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.Barcodes);
        testParsing(factory, rs, metricOne, barcodePosition);
    }

    @Test
    public void testDualBarcodes() throws Exception {
        final File metricsFile = File.createTempFile("dual.", ".metrics");
        metricsFile.deleteOnExit();

        final String[] args = new String[]{
                "BASECALLS_DIR=" + dual.getAbsolutePath(),
                "LANE=" + 1,
                "METRICS_FILE=" + metricsFile.getPath(),
                "READ_STRUCTURE=" + "25T8B8B25T",
                "BARCODE=" + "CAATAGTCCGACTCTC"
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);
        final MetricsFile<BarcodeMetric, Integer> result = new MetricsFile<>();
        result.read(new FileReader(metricsFile));
        Assert.assertEquals(result.getMetrics().get(0).PERFECT_MATCHES, 1, "Got wrong number of perfect matches");
        Assert.assertEquals(result.getMetrics().get(0).ONE_MISMATCH_MATCHES, 0, "Got wrong number of one-mismatch matches");
    }

    /**
     * Testing the quality thresholding. Looking at a single barcode (ACAGTG) with a min quality of 25 and no mismatches
     */
    @Test(dataProvider = "qualityBarcodeData")
    public void testQualityBarcodes(final int quality,
                                    final int maxMismatches, final int perfectMatches, final int oneMismatch,
                                    final String testName) throws Exception {
        final File metricsFile = File.createTempFile("qual.", ".metrics");
        metricsFile.deleteOnExit();

        final String[] args = new String[]{
                "BASECALLS_DIR=" + qual.getPath(),
                "LANE=" + 1,
                "READ_STRUCTURE=25T8B25T",
                "METRICS_FILE=" + metricsFile.getPath(),
                "MINIMUM_BASE_QUALITY=" + quality,
                "MAX_MISMATCHES=" + maxMismatches,
                "BARCODE=CAATAGTC"
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);
        final MetricsFile<BarcodeMetric, Integer> result = new MetricsFile<>();
        result.read(new FileReader(metricsFile));
        Assert.assertEquals(result.getMetrics().get(0).PERFECT_MATCHES, perfectMatches, "Got wrong number of perfect matches for test: '" + testName + "'");
        Assert.assertEquals(result.getMetrics().get(0).ONE_MISMATCH_MATCHES, oneMismatch, "Got wrong number of one-mismatch matches for test: '" + testName + "'");
    }

    @DataProvider(name = "qualityBarcodeData")
    public Object[][] getQualityTestData() {
        return new Object[][]{
                {16, 0, 1, 0, "Barcode has good quality, 1 match"},
                {25, 0, 0, 0, "Barcode has quality failures, no matches"}
        };
    }

    private void testParsing(final IlluminaDataProviderFactory factory, final ReadStructure readStructure, final BarcodeMetric metricACAGTG, final int barcodePosition) {

        int numReads = 0;

        final BaseIlluminaDataProvider dataProvider = factory.makeDataProvider();
        while (dataProvider.hasNext()) {
            final ClusterData cluster = dataProvider.next();

            if (metricACAGTG.BARCODE.equals(cluster.getMatchedBarcode())) {
                ++numReads;
            }

            Assert.assertEquals(cluster.getRead(readStructure.templates.getIndices()[0]).getQualities().length, barcodePosition - 1);
            Assert.assertEquals(cluster.getRead(readStructure.templates.getIndices()[0]).getBases().length, barcodePosition - 1);
        }
        Assert.assertEquals(numReads, metricACAGTG.READS);
        dataProvider.close();
    }

    private MetricsFile<BarcodeMetric, Integer> runIt(final int lane, final String readStructure)
            throws Exception {
        final File metricsFile = File.createTempFile("eib.", ".metrics");
        metricsFile.deleteOnExit();

        final List<String> args = new ArrayList<String>(Arrays.asList(
                "BASECALLS_DIR=" + basecallsDir.getPath(),
                "LANE=" + lane,
                "READ_STRUCTURE=" + readStructure,
                "METRICS_FILE=" + metricsFile.getPath()
        ));
        for (final String barcode : BARCODES) {
            args.add("BARCODE=" + barcode);
        }
        return runIt(args, metricsFile);
    }

    private MetricsFile<BarcodeMetric, Integer> runIt(final List<String> args, final File metricsFile) throws Exception {
        // Generate _barcode.txt files and metrics file.
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<BarcodeMetric, Integer> retval = new MetricsFile<>();
        retval.read(new FileReader(metricsFile));
        return retval;
    }
}
