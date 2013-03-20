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

import net.sf.picard.io.IoUtil;
import net.sf.samtools.util.BufferedLineReader;
import net.sf.samtools.util.LineReader;
import net.sf.samtools.util.StringUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Run IlluminaBasecallsToSam in various barcode & non-barcode modes
 *
 * @author alecw@broadinstitute.org
 */
public class IlluminaBasecallsToSamTest {

    private static final File BASECALLS_DIR = new File("testdata/net/sf/picard/illumina/IlluminaTests/BasecallsDir");
    private static final File TEST_DATA_DIR = new File("testdata/net/sf/picard/illumina/IlluminaBasecallsToSamTest");

    @Test
    public void testTileNumberComparator() {
        Assert.assertTrue(IlluminaBasecallsConverter.TILE_NUMBER_COMPARATOR.compare(100, 10) < 0, "");
        Assert.assertTrue(IlluminaBasecallsConverter.TILE_NUMBER_COMPARATOR.compare(20, 200) > 0, "");
        Assert.assertTrue(IlluminaBasecallsConverter.TILE_NUMBER_COMPARATOR.compare(10, 10) == 0, "");
    }

    @Test
    public void testNonBarcoded() throws Exception {
        final File outputBam = File.createTempFile("nonBarcoded.", ".sam");
        outputBam.deleteOnExit();
        final int lane = 1;
        new IlluminaBasecallsToSam().instanceMain(new String[]{
                "BASECALLS_DIR=" + BASECALLS_DIR,
                "LANE=" + lane,
                "READ_STRUCTURE=76T76T",
                "OUTPUT=" + outputBam,
                "RUN_BARCODE=HiMom",
                "SAMPLE_ALIAS=HiDad",
                "LIBRARY_NAME=Hello, World"
        });
        IoUtil.assertFilesEqual(outputBam, new File(TEST_DATA_DIR, "nonBarcoded.sam"));
    }

    @Test
    public void testMultiplexed() throws Exception {
        runStandardTest(7, "multiplexedBarcode.", "barcode.params", 1, "30T8B");
    }

    //Same as testMultiplexed except we use BARCODE_1 instead of BARCODE
    @Test
    public void testMultiplexedWithAlternateBarcodeName() throws Exception {
        runStandardTest(7, "singleBarcodeAltName.", "multiplexed_positive_rgtags.params", 1, "30T8B");
    }

    @Test
    public void testDualBarcodes() throws Exception {
        runStandardTest(9, "dualBarcode.", "barcode_double.params", 2, "30T8B8B");
    }

    /**
     * Ensures that a run missing a barcode from the parameters file throws an error.
     * 
     * TODO: This testcase isn't broken, but can spawn an issue with FileChannelJDKBugWorkAround since it expects
     * an exception to be thrown.
     */
    @Test(groups={"broken"})
    public void testCorruptDataReturnCode() throws Exception {
        boolean exceptionThrown = false;
        try {
            runStandardTest(9, "dualBarcode.", "negative_test.params", 2, "30T8B8B");
        } catch (Throwable e) {
            exceptionThrown = true;
        } finally {
            Assert.assertTrue(exceptionThrown);
        }
    }

    /**
     * This test utility takes a libraryParamsFile and generates output sam files through IlluminaBasecallsToSam to compare against
     * preloaded test data
     *
     * @param jobName
     * @param libraryParamsFile
     * @param concatNColumnFields
     * @param readStructure
     * @throws Exception
     */
    private void runStandardTest(final int lane, final String jobName, final String libraryParamsFile, final int concatNColumnFields, final String readStructure) throws Exception {
        final File outputDir = File.createTempFile(jobName, ".dir");
        outputDir.delete();
        outputDir.mkdir();
        outputDir.deleteOnExit();
        // Create barcode.params with output files in the temp directory
        final File libraryParams = new File(outputDir, libraryParamsFile);
        libraryParams.deleteOnExit();
        final List<File> samFiles = new ArrayList<File>();
        final LineReader reader = new BufferedLineReader(new FileInputStream(new File(TEST_DATA_DIR, libraryParamsFile)));
        final PrintWriter writer = new PrintWriter(libraryParams);
        final String header = reader.readLine();
        writer.println(header + "\tOUTPUT");
        while (true) {
            final String line = reader.readLine();
            if (line == null) {
                break;
            }
            final String[] fields = line.split("\t");
            final File outputSam = new File(outputDir, StringUtil.join("", Arrays.copyOfRange(fields, 0, concatNColumnFields)) + ".sam");
            outputSam.deleteOnExit();
            samFiles.add(outputSam);
            writer.println(line + "\t" + outputSam);
        }
        writer.close();

        new IlluminaBasecallsToSam().instanceMain(new String[]{
                "BASECALLS_DIR=" + BASECALLS_DIR,
                "LANE=" + lane,
                "RUN_BARCODE=HiMom",
                "READ_STRUCTURE=" + readStructure,
                "LIBRARY_PARAMS=" + libraryParams
        });

        for (final File outputSam : samFiles) {
            IoUtil.assertFilesEqual(outputSam, new File(TEST_DATA_DIR, outputSam.getName()));
        }
    }
}
