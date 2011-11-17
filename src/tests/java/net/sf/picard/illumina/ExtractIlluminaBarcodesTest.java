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
package net.sf.picard.illumina;

import net.sf.picard.util.BasicInputParser;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.AfterTest;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.io.File;
import java.io.FileReader;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

import net.sf.picard.io.IoUtil;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.illumina.parser.*;

/**
 * @author alecw@broadinstitute.org
 */
public class ExtractIlluminaBarcodesTest {
    private static final File TEST_DATA_DIR = new File("testdata/net/sf/picard/illumina/ExtractIlluminaBarcodes");
    private static final String[] BARCODES = {
            "ACAGTG",
            "ACAGTT", // This one is artificial -- one edit away from the first one
            "ACTTGA",
            "ATCACG",
            "CAGATC",
            "CGATGT",
            "CTTGTA",
            "GATCAG",
            "GCCAAT",
            "TAGCTT",
            "TGACCA"
    };

    private File basecallsDir;

    @BeforeTest
    private void setUp() throws Exception {
        basecallsDir = File.createTempFile("eib.", ".tmp");
        Assert.assertTrue(basecallsDir.delete());
        Assert.assertTrue(basecallsDir.mkdir());
        for (final File source : TEST_DATA_DIR.listFiles()) {
            if (!source.isFile()) {
                continue;
            }
            final File dest = new File(basecallsDir, source.getName());
            IoUtil.copyFile(source, dest);
        }
    }

    @AfterTest
    private void tearDown() {
        IoUtil.deleteDirectoryTree(basecallsDir);
    }

    @Test
    public void testSingleEndWithBarcodeAtStart() throws Exception {
        final MetricsFile<ExtractIlluminaBarcodes.BarcodeMetric, Integer> metricsFile = runIt(1, 1);
        Assert.assertEquals(metricsFile.getMetrics().get(0).PERFECT_MATCHES, 1);
    }

    @Test
    public void testSingleEndWithBarcodeAtEnd() throws Exception {
        final MetricsFile<ExtractIlluminaBarcodes.BarcodeMetric, Integer> metricsFile = runIt(2, 37);
        Assert.assertEquals(metricsFile.getMetrics().get(0).PERFECT_MATCHES, 1);
    }

    @Test
    public void testPairedEndWithBarcodeOnFirstEnd() throws Exception {
        final MetricsFile<ExtractIlluminaBarcodes.BarcodeMetric, Integer> metricsFile = runIt(3, 37);
        Assert.assertEquals(metricsFile.getMetrics().get(0).PERFECT_MATCHES, 1);
    }

    @Test
    public void testPairedEndWithBarcodeOnSecondEnd() throws Exception {
        final MetricsFile<ExtractIlluminaBarcodes.BarcodeMetric, Integer> metricsFile = runIt(4, 73);
        Assert.assertEquals(metricsFile.getMetrics().get(0).PERFECT_MATCHES, 1);
    }

    /**
     * 4 cases tested:
     * * exact match to ACAGTG
     * * inexact match within threshold to TGACCA
     * * inexact match not within threshold to TGACCA
     * * inexact match where the next match is too close to ACAGTG
     * @throws Exception
     */
    @Test
    public void testBarcodeMatching() throws Exception {
        final int lane = 5;
        final int barcodePosition = 37;
        final MetricsFile<ExtractIlluminaBarcodes.BarcodeMetric, Integer> metricsFile = runIt(lane, barcodePosition);

        ExtractIlluminaBarcodes.BarcodeMetric metricACAGTG = null;
        ExtractIlluminaBarcodes.BarcodeMetric metricTGACCA = null;
        ExtractIlluminaBarcodes.BarcodeMetric metricNoMatch = null;
        for (final ExtractIlluminaBarcodes.BarcodeMetric metric : metricsFile.getMetrics()) {
            if (metric.BARCODE.equals("ACAGTG")) {
                metricACAGTG = metric;
            } else if (metric.BARCODE.equals("TGACCA")) {
                metricTGACCA = metric;
            } else if (metric.BARCODE.equals("NNNNNN")) {
                metricNoMatch = metric;
            }
        }
        Assert.assertEquals(metricACAGTG.PERFECT_MATCHES, 1);
        Assert.assertEquals(metricACAGTG.ONE_MISMATCH_MATCHES, 0);
        Assert.assertEquals(metricACAGTG.PF_READS, 1);
        Assert.assertEquals(metricACAGTG.READS, 1);

        for (final ExtractIlluminaBarcodes.BarcodeMetric metric : metricsFile.getMetrics()) {
            if (metric == metricACAGTG || metric == metricTGACCA || metric == metricNoMatch) {
                continue;
            }
            Assert.assertEquals(metric.READS, 0);
        }

        // one inexact match
        Assert.assertEquals(metricTGACCA.READS, 1);
        Assert.assertEquals(metricTGACCA.ONE_MISMATCH_MATCHES, 1);

        Assert.assertEquals(metricNoMatch.READS, 2);
        Assert.assertEquals(metricNoMatch.PF_READS, 1);

        // Check the barcode files themselves
        final File[] barcodeFiles = IoUtil.getFilesMatchingRegexp(basecallsDir, "s_" + lane + "_\\d{4}_barcode.txt");
        Arrays.sort(barcodeFiles);

        final BasicInputParser barcodeParser = new BasicInputParser(true, barcodeFiles);

        // Exact match
        String[] illuminaFields = barcodeParser.next();
        Assert.assertEquals(illuminaFields[1], "Y");
        Assert.assertEquals(illuminaFields[2], "ACAGTG");

        // Inexact match
        illuminaFields = barcodeParser.next();
        Assert.assertEquals(illuminaFields[1], "Y");
        Assert.assertEquals(illuminaFields[2], "TGACCA");

        // Too many mismatches
        illuminaFields = barcodeParser.next();
        Assert.assertEquals(illuminaFields[1], "N");

        // Next match too close
        illuminaFields = barcodeParser.next();
        Assert.assertEquals(illuminaFields[1], "N");

        Assert.assertFalse(barcodeParser.hasNext());
        barcodeParser.close();

        // Tack on test of barcode-informed Illumina Basecall parsing
        final IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(basecallsDir, lane, barcodePosition,
                metricACAGTG.BARCODE.length(), IlluminaDataType.BaseCalls, IlluminaDataType.QualityScores,
                IlluminaDataType.Barcodes);
        int numReads = 0;

        final IlluminaDataProvider dataProvider = factory.makeDataProvider();
        while (dataProvider.hasNext()) {
            final ClusterData cluster = dataProvider.next();

            if(metricACAGTG.BARCODE.equals(cluster.getMatchedBarcode())) {
                ++numReads;
            }
            
            Assert.assertEquals(getFirstTemplate(cluster).getQualities().length, barcodePosition - 1);
            Assert.assertEquals(getFirstTemplate(cluster).getBases().length, barcodePosition - 1);
        }
        Assert.assertEquals(numReads, metricACAGTG.READS);
    }

    //TODO: When we extract out the IlluminaRunConfig this should be eliminated
    private static ReadData getFirstTemplate(ClusterData cd) {
        for(int i = 0; i < cd.getNumReads(); i++) {
            final ReadData rd = cd.getRead(i);
            if(rd.getReadType() == ReadType.Template)
                return rd;
        }

        return null;
    }

    private MetricsFile<ExtractIlluminaBarcodes.BarcodeMetric, Integer> runIt(final int lane, final int position)
            throws Exception {
        final File metricsFile = File.createTempFile("eib.", ".metrics");
        metricsFile.deleteOnExit();
        final List<String> args = new ArrayList<String>(Arrays.asList(
                "BASECALLS_DIR=" + basecallsDir.getPath(),
                "LANE=" + lane,
                "BARCODE_POSITION=" + position,
                "METRICS_FILE=" + metricsFile.getPath()
                ));
        for (final String barcode : BARCODES) {
            args.add("BARCODE=" + barcode);
        }
        // Generate _barcode.txt files and metrics file.
        Assert.assertEquals(new ExtractIlluminaBarcodes().instanceMain(args.toArray(new String[args.size()])), 0);

        final MetricsFile<ExtractIlluminaBarcodes.BarcodeMetric,Integer> retval =  new MetricsFile<ExtractIlluminaBarcodes.BarcodeMetric,Integer>();
        retval.read(new FileReader(metricsFile));
        return retval;
    }
}
