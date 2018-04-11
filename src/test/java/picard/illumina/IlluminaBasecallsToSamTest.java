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

import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.LineReader;
import htsjdk.samtools.util.StringUtil;
import htsjdk.samtools.util.TestUtil;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Run IlluminaBasecallsToSam in various barcode & non-barcode modes
 *
 * @author alecw@broadinstitute.org
 */
public class IlluminaBasecallsToSamTest extends CommandLineProgramTest {

    private static final File BASECALLS_DIR = new File("testdata/picard/illumina/25T8B25T/Data/Intensities/BaseCalls");
    private static final File DUAL_BASECALLS_DIR = new File("testdata/picard/illumina/25T8B8B25T/Data/Intensities/BaseCalls");
    private static final File TEST_DATA_DIR = new File("testdata/picard/illumina/25T8B25T/sams");
    private static final File DUAL_TEST_DATA_DIR = new File("testdata/picard/illumina/25T8B8B25T/sams");
    private static final File TEST_DATA_DIR_WITH_4M_INDEX = new File("testdata/picard/illumina/25T8B25T/sams_with_4M");
    private static final File TEST_DATA_DIR_WITH_4M4M_INDEX = new File("testdata/picard/illumina/25T8B25T/sams_with_4M4M");
    private static final File TEST_DATA_DIR_WITH_CBCLS = new File("testdata/picard/illumina/151T8B8B151T_cbcl/Data/Intensities/BaseCalls");
    private static final File DUAL_CBCL_TEST_DATA_DIR = new File("testdata/picard/illumina/151T8B8B151T_cbcl/sams");
    private static final File TEST_DATA_HISEQX_SINGLE_LOCS = new File("testdata/picard/illumina/25T8B8B25T_hiseqx/Data/Intensities/BaseCalls");
    private static final File HISEQX_TEST_DATA_DIR = new File("testdata/picard/illumina/25T8B8B25T_hiseqx/sams");

    public String getCommandLineProgramName() {
        return IlluminaBasecallsToSam.class.getSimpleName();
    }

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

        Assert.assertEquals(runPicardCommandLine(new String[]{
                "BASECALLS_DIR=" + BASECALLS_DIR,
                "LANE=" + lane,
                "READ_STRUCTURE=25S8S25T",
                "OUTPUT=" + outputBam,
                "RUN_BARCODE=HiMom",
                "SAMPLE_ALIAS=HiDad",
                "LIBRARY_NAME=Hello, World"
        }), 0);
        IOUtil.assertFilesEqual(outputBam, new File(TEST_DATA_DIR, "nonBarcoded.sam"));
    }

    @Test
    public void testNonBarcodedWithMolecularIndex() throws Exception {
        final File outputBam = File.createTempFile("nonBarcodedWithMI.", ".sam");
        outputBam.deleteOnExit();
        final int lane = 1;

        Assert.assertEquals(runPicardCommandLine(new String[]{
                "BASECALLS_DIR=" + BASECALLS_DIR,
                "LANE=" + lane,
                "READ_STRUCTURE=25S8M25T",
                "OUTPUT=" + outputBam,
                "RUN_BARCODE=HiMom",
                "SAMPLE_ALIAS=HiDad",
                "LIBRARY_NAME=Hello, World"
        }), 0);
        IOUtil.assertFilesEqual(outputBam, new File(TEST_DATA_DIR, "nonBarcodedWithMolecularIndex8M.sam"));
    }

    @Test
    public void testNonBarcodedWithTagPerMolecularIndexIsNUll() throws Exception {
        final File outputBam = File.createTempFile("nonBarcodedWithMI.", ".sam");
        outputBam.deleteOnExit();
        final int lane = 1;

        Assert.assertEquals(runPicardCommandLine(new String[]{
                "BASECALLS_DIR=" + BASECALLS_DIR,
                "LANE=" + lane,
                "READ_STRUCTURE=25S8M25T",
                "OUTPUT=" + outputBam,
                "RUN_BARCODE=HiMom",
                "SAMPLE_ALIAS=HiDad",
                "LIBRARY_NAME=Hello, World",
                "TAG_PER_MOLECULAR_INDEX=null"
        }), 0);
        IOUtil.assertFilesEqual(outputBam, new File(TEST_DATA_DIR, "nonBarcodedWithMolecularIndex8M.sam"));
    }

    @Test
    public void testNonBarcodedWithDualMoleclarIndex() throws Exception {
        final File outputBam = File.createTempFile("nonBarcodedWithDualMI.", ".sam");
        outputBam.deleteOnExit();
        final int lane = 1;

        Assert.assertEquals(runPicardCommandLine(new String[]{
                "BASECALLS_DIR=" + BASECALLS_DIR,
                "LANE=" + lane,
                "READ_STRUCTURE=25S4M4M25T",
                "OUTPUT=" + outputBam,
                "RUN_BARCODE=HiMom",
                "SAMPLE_ALIAS=HiDad",
                "LIBRARY_NAME=Hello, World"
        }), 0);
        IOUtil.assertFilesEqual(outputBam, new File(TEST_DATA_DIR, "nonBarcodedWithMolecularIndex4M4M.sam"));
    }

    // This *should* store molecular indexes individually in ZA and ZB
    @Test
    public void testNonBarcodedWithTagPerMolecularIndexDual() throws Exception {
        final File outputBam = File.createTempFile("nonBarcodedWithDualMI.", ".sam");
        outputBam.deleteOnExit();
        final int lane = 1;

        Assert.assertEquals(runPicardCommandLine(new String[]{
                "BASECALLS_DIR=" + BASECALLS_DIR,
                "LANE=" + lane,
                "READ_STRUCTURE=25S4M4M25T",
                "OUTPUT=" + outputBam,
                "RUN_BARCODE=HiMom",
                "SAMPLE_ALIAS=HiDad",
                "LIBRARY_NAME=Hello, World",
                "TAG_PER_MOLECULAR_INDEX=ZA",
                "TAG_PER_MOLECULAR_INDEX=ZB"

        }), 0);
        IOUtil.assertFilesEqual(outputBam, new File(TEST_DATA_DIR, "nonBarcodedWithTagPerMolecularIndex4M4M.sam"));
    }

    // Too many tags
    @Test
    public void testNonBarcodedWithTagPerMolecularIndexDualTooManyTags() throws Exception {
        final File outputBam = File.createTempFile("nonBarcodedWithDualMI.", ".sam");
        outputBam.deleteOnExit();
        final int lane = 1;

        Assert.assertEquals(runPicardCommandLine(new String[]{
                "BASECALLS_DIR=" + BASECALLS_DIR,
                "LANE=" + lane,
                "READ_STRUCTURE=25S4M4M25T",
                "OUTPUT=" + outputBam,
                "RUN_BARCODE=HiMom",
                "SAMPLE_ALIAS=HiDad",
                "LIBRARY_NAME=Hello, World",
                "TAG_PER_MOLECULAR_INDEX=ZA",
                "TAG_PER_MOLECULAR_INDEX=ZB",
                "TAG_PER_MOLECULAR_INDEX=ZC"

        }), 1);
    }

    // Too few tags
    @Test
    public void testNonBarcodedWithTagPerMolecularIndexDualTooFewTags() throws Exception {
        final File outputBam = File.createTempFile("nonBarcodedWithDualMI.", ".sam");
        outputBam.deleteOnExit();
        final int lane = 1;

        Assert.assertEquals(runPicardCommandLine(new String[]{
                "BASECALLS_DIR=" + BASECALLS_DIR,
                "LANE=" + lane,
                "READ_STRUCTURE=25S2M2M2M2M25T",
                "OUTPUT=" + outputBam,
                "RUN_BARCODE=HiMom",
                "SAMPLE_ALIAS=HiDad",
                "LIBRARY_NAME=Hello, World",
                "TAG_PER_MOLECULAR_INDEX=ZA",
                "TAG_PER_MOLECULAR_INDEX=ZB",
                "TAG_PER_MOLECULAR_INDEX=ZC"

        }), 1);
    }

    // Just the right number of tags
    @Test
    public void testNonBarcodedWithTagPerMolecularIndexDualFourMolecularIndexes() throws Exception {
        final File outputBam = File.createTempFile("nonBarcodedWithDualMI.", ".sam");
        outputBam.deleteOnExit();
        final int lane = 1;

        Assert.assertEquals(runPicardCommandLine(new String[]{
                "BASECALLS_DIR=" + BASECALLS_DIR,
                "LANE=" + lane,
                "READ_STRUCTURE=25S2M2M2M2M25T",
                "OUTPUT=" + outputBam,
                "RUN_BARCODE=HiMom",
                "SAMPLE_ALIAS=HiDad",
                "LIBRARY_NAME=Hello, World",
                "TAG_PER_MOLECULAR_INDEX=ZA",
                "TAG_PER_MOLECULAR_INDEX=ZB",
                "TAG_PER_MOLECULAR_INDEX=ZC",
                "TAG_PER_MOLECULAR_INDEX=ZD"
        }), 0);
        IOUtil.assertFilesEqual(outputBam, new File(TEST_DATA_DIR, "nonBarcodedWithTagPerMolecularIndex2M2M2M2M.sam"));
    }

    @Test
    public void testNonBarcodedTileSubset() throws Exception {
        final File outputBam = File.createTempFile("nonBarcoded.", ".sam");
        outputBam.deleteOnExit();
        final int lane = 1;

        Assert.assertEquals(runPicardCommandLine(new String[]{
                "BASECALLS_DIR=" + BASECALLS_DIR,
                "LANE=" + lane,
                "READ_STRUCTURE=25S8S25T",
                "OUTPUT=" + outputBam,
                "RUN_BARCODE=HiMom",
                "SAMPLE_ALIAS=HiDad",
                "LIBRARY_NAME=Hello, World",
                "FIRST_TILE=1201",
                "TILE_LIMIT=1"
        }), 0);
        IOUtil.assertFilesEqual(outputBam, new File(TEST_DATA_DIR, "nonBarcodedTileSubset.sam"));
    }

    @Test
    public void testMultiplexed() throws Exception {
        runStandardTest(1, "multiplexedBarcode.", "barcode.params", 1, "25T8B25T", BASECALLS_DIR, TEST_DATA_DIR);
    }

    @Test
    public void testMultiplexedWith4MIndex() throws Exception {
        runStandardTest(1, "multiplexedBarcode.", "barcode.params", 1, "25T8B4M21T", BASECALLS_DIR, TEST_DATA_DIR_WITH_4M_INDEX);
    }

    @Test
    public void testMultiplexedWith4M4MIndex() throws Exception {
        runStandardTest(1, "multiplexedBarcode2.", "barcode.params", 1, "25T8B4M4M17T", BASECALLS_DIR, TEST_DATA_DIR_WITH_4M4M_INDEX);
    }

    //Same as testMultiplexed except we use BARCODE_1 instead of BARCODE
    @Test
    public void testMultiplexedWithAlternateBarcodeName() throws Exception {
        runStandardTest(1, "singleBarcodeAltName.", "multiplexed_positive_rgtags.params", 1, "25T8B25T", BASECALLS_DIR, TEST_DATA_DIR);
    }

    @Test
    public void testDualBarcodes() throws Exception {
        runStandardTest(1, "dualBarcode.", "barcode_double.params", 2, "25T8B8B25T", DUAL_BASECALLS_DIR, DUAL_TEST_DATA_DIR);
    }

    @Test
    public void testCbclConvert() throws Exception {
        runStandardTest(1, "cbclConvert.", "barcode_double.params", 2, "151T8B8B151T", TEST_DATA_DIR_WITH_CBCLS, DUAL_CBCL_TEST_DATA_DIR);
    }

    @Test
    public void testHiseqxSingleLocs() throws Exception {
        runStandardTest(1, "hiseqxSingleLocs.", "barcode_double.params", 2, "25T8B8B25T", TEST_DATA_HISEQX_SINGLE_LOCS, HISEQX_TEST_DATA_DIR);
    }

    @Test
    public void testDualBarcodesSingleTile() throws Exception {
        runStandardTest(1, "dualBarcode.", "barcode_double.params", 2, "25T8B8B25T", DUAL_BASECALLS_DIR, DUAL_TEST_DATA_DIR, 1101);
    }

    @Test
    public void testCbclConvertSingleTile() throws Exception {
        runStandardTest(1, "cbclConvert.", "barcode_double.params", 2, "151T8B8B151T", TEST_DATA_DIR_WITH_CBCLS, DUAL_CBCL_TEST_DATA_DIR, 1102);
    }

    /**
     * Ensures that a run missing a barcode from the parameters file throws an error.
     * <p>
     * TODO: This testcase isn't broken, but can spawn an issue with FileChannelJDKBugWorkAround since it expects
     * an exception to be thrown.
     */
    @Test(groups = {"broken"})
    public void testCorruptDataReturnCode() throws Exception {
        boolean exceptionThrown = false;
        try {
            runStandardTest(9, "dualBarcode.", "negative_test.params", 2, "30T8B8B", BASECALLS_DIR, TEST_DATA_DIR);
        } catch (Throwable e) {
            exceptionThrown = true;
        } finally {
            Assert.assertTrue(exceptionThrown);
        }
    }

    private void runStandardTest(final int lane, final String jobName, final String libraryParamsFile,
                                 final int concatNColumnFields, final String readStructure,
                                 final File baseCallsDir, final File testDataDir) throws Exception {
        runStandardTest(lane, jobName, libraryParamsFile, concatNColumnFields, readStructure, baseCallsDir, testDataDir, null);
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
    private void runStandardTest(final int lane, final String jobName, final String libraryParamsFile,
                                 final int concatNColumnFields, final String readStructure,
                                 final File baseCallsDir, final File testDataDir, final Integer tile) throws Exception {
        final Path outputDir = Files.createTempDirectory(jobName);
        try {
            final String tilePrefix = (tile != null) ? tile + "." : "";

            // Create barcode.params with output files in the temp directory
            final File libraryParams = new File(outputDir.toFile(), libraryParamsFile);
            libraryParams.deleteOnExit();
            final List<File> samFiles = new ArrayList<File>();
            final LineReader reader = new BufferedLineReader(new FileInputStream(new File(testDataDir, libraryParamsFile)));
            final PrintWriter writer = new PrintWriter(libraryParams);
            final String header = reader.readLine();
            writer.println(header + "\tOUTPUT");
            while (true) {
                final String line = reader.readLine();
                if (line == null) {
                    break;
                }
                final String[] fields = line.split("\t");
                final File outputSam = new File(outputDir.toFile(), StringUtil.join("", Arrays.copyOfRange(fields, 0, concatNColumnFields)) + ".sam");
                outputSam.deleteOnExit();
                samFiles.add(new File(outputSam.getParentFile(), tilePrefix + outputSam.getName()));
                writer.println(line + "\t" + outputSam);
            }
            writer.close();
            reader.close();

            List<String> args = new ArrayList<>();
            args.add("BASECALLS_DIR=" + baseCallsDir);
            args.add("LANE=" + lane);
            args.add("RUN_BARCODE=HiMom");
            args.add("READ_STRUCTURE=" + readStructure);
            args.add("LIBRARY_PARAMS=" + libraryParams);
            if (tile != null) {
                args.add("PROCESS_SINGLE_TILE=" + tile);
            }

            Assert.assertEquals(runPicardCommandLine(args), 0);

            for (final File outputSam : samFiles) {
                IOUtil.assertFilesEqual(outputSam, new File(testDataDir, outputSam.getName()));
            }
        } finally {
            TestUtil.recursiveDelete(outputDir.toFile());
            Files.delete(outputDir);
        }
    }
}
