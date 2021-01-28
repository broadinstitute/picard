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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamFileValidator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.*;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.sam.util.SamComparison;

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

    private static final File ILLUMINA_TEST_DIR = new File("testdata/picard/illumina/");
    private static final File BASECALLS_DIR = new File(ILLUMINA_TEST_DIR, "25T8B25T/Data/Intensities/BaseCalls");
    private static final File DUAL_BASECALLS_DIR = new File(ILLUMINA_TEST_DIR, "25T8B8B25T/Data/Intensities/BaseCalls");
    private static final File TEST_DATA_DIR = new File(ILLUMINA_TEST_DIR, "25T8B25T/sams");
    private static final File DUAL_TEST_DATA_DIR = new File(ILLUMINA_TEST_DIR, "25T8B8B25T/sams");
    private static final File TEST_DATA_DIR_WITH_4M_INDEX = new File(ILLUMINA_TEST_DIR, "25T8B25T/sams_with_4M");
    private static final File TEST_DATA_DIR_WITH_4M4M_INDEX = new File(ILLUMINA_TEST_DIR, "25T8B25T/sams_with_4M4M");
    private static final File TEST_DATA_DIR_WITH_CBCLS = new File(ILLUMINA_TEST_DIR, "151T8B8B151T_cbcl/Data/Intensities/BaseCalls");
    private static final File DUAL_CBCL_TEST_DATA_DIR = new File(ILLUMINA_TEST_DIR, "151T8B8B151T_cbcl/sams");
    private static final File TEST_DATA_HISEQX_SINGLE_LOCS = new File(ILLUMINA_TEST_DIR, "25T8B8B25T_hiseqx/Data/Intensities/BaseCalls");
    private static final File HISEQX_TEST_DATA_DIR = new File(ILLUMINA_TEST_DIR, "25T8B8B25T_hiseqx/sams");

    public String getCommandLineProgramName() {
        return IlluminaBasecallsToSam.class.getSimpleName();
    }

    @Test
    public void testTileNumberComparator() {
        Assert.assertTrue(BasecallsConverter.TILE_NUMBER_COMPARATOR.compare(100, 10) < 0, "");
        Assert.assertTrue(BasecallsConverter.TILE_NUMBER_COMPARATOR.compare(20, 200) > 0, "");
        Assert.assertTrue(BasecallsConverter.TILE_NUMBER_COMPARATOR.compare(10, 10) == 0, "");
    }

    @Test
    public void testNonBarcodedWithCenter() throws Exception {
        for (final boolean sort : new boolean[]{false, true}) {
            final File outputBam = File.createTempFile("nonBarcodedDescriptionNonBI.", ".sam");
            outputBam.deleteOnExit();
            final int lane = 1;

            Assert.assertEquals(runPicardCommandLine(new String[]{
                    "BASECALLS_DIR=" + BASECALLS_DIR,
                    "LANE=" + lane,
                    "READ_STRUCTURE=25S8S25T",
                    "OUTPUT=" + outputBam,
                    "RUN_BARCODE=HiMom",
                    "SAMPLE_ALIAS=HiDad",
                    "SEQUENCING_CENTER=TEST_CENTER123",
                    "LIBRARY_NAME=Hello, World",
                    "SORT=" + sort
            }), 0);
            final SamFileValidator validator = new SamFileValidator(new PrintWriter(System.out), 100);
            validator.validateSamFileSummary(SamReaderFactory.makeDefault().open(outputBam), null);
            String fileName = sort ? "nonBarcodedDescriptionNonBI.sam" : "nonBarcodedDescriptionNonBI.unsorted.sam"  ;
            final File expectedSamFile = new File(TEST_DATA_DIR, fileName);

            if (sort) {
                IOUtil.assertFilesEqual(outputBam, expectedSamFile);
            } else {
                compareSams(outputBam, expectedSamFile);
            }
        }
    }

    @DataProvider
    public Object[][] molecularBarcodeData() {
        return new Object[][]{
                {"25S8S25T", null, null, "nonBarcoded", 0},
                {"25S8M25T", null, null, "nonBarcodedWithMolecularIndex8M", 0},
                {"25S8M25T", null, new String[]{"TAG_PER_MOLECULAR_INDEX=null"}, "nonBarcodedWithMolecularIndex8M", 0},
                {"25S4M4M25T", null, null, "nonBarcodedWithMolecularIndex4M4M", 0},
                {"25S4M4M25T", new String[]{"ZA", "ZB"}, null, "nonBarcodedWithTagPerMolecularIndex4M4M", 0},
                {"25S4M4M25T", new String[]{"ZA", "ZB", "ZC"}, null, null, 1},
                {"25S2M2M2M2M25T", new String[]{"ZA", "ZB", "ZC"}, null, null, 1},
                {"25S2M2M2M2M25T", new String[]{"ZA", "ZB", "ZC", "ZD"}, null, "nonBarcodedWithTagPerMolecularIndex2M2M2M2M", 0},
                {"25S8S25T", null, new String[]{"FIRST_TILE=1201", "TILE_LIMIT=1"}, "nonBarcodedTileSubset", 0}
        };
    }

    @Test(dataProvider = "molecularBarcodeData")
    public void testMolecularBarcodes(final String readStructure, final String[] umiTags, final String[] extraArgs, final String expectedSam, final int expectedReturn) throws Exception {
        for (final boolean sort : new boolean[]{false, true}) {
            final File outputBam = File.createTempFile("molecularBarcodeTest.", ".sam");
            outputBam.deleteOnExit();
            final int lane = 1;

            List<String> args = new ArrayList<>(CollectionUtil.makeList(
                    "BASECALLS_DIR=" + BASECALLS_DIR,
                    "LANE=" + lane,
                    "READ_STRUCTURE=" + readStructure,
                    "OUTPUT=" + outputBam,
                    "RUN_BARCODE=HiMom",
                    "SAMPLE_ALIAS=HiDad",
                    "SEQUENCING_CENTER=BI",
                    "LIBRARY_NAME=Hello, World",
                    "SORT=" + sort));

            if (umiTags != null) {
                for (final String umiTag : umiTags) {
                    args.add("TAG_PER_MOLECULAR_INDEX=" + umiTag);
                }
            }

            if (extraArgs != null) {
                args.addAll(Arrays.asList(extraArgs));
            }

            Assert.assertEquals(runPicardCommandLine(args), expectedReturn);
            if (expectedSam != null) {
                final SamFileValidator validator = new SamFileValidator(new PrintWriter(System.out), 100);
                validator.validateSamFileSummary(SamReaderFactory.makeDefault().open(outputBam), null);
                String expectedSamFilename = sort ? expectedSam + ".sam" : expectedSam + ".unsorted.sam";
                final File expectedSamFile = new File(TEST_DATA_DIR, expectedSamFilename);
                if (sort) {
                    IOUtil.assertFilesEqual(outputBam, expectedSamFile);
                } else {
                    compareSams(outputBam, expectedSamFile);
                }
            }
        }
    }

    @DataProvider
    public Object[][] multiplexedData() {
        final File testDirParent = TEST_DATA_DIR.getParentFile();
        return new Object[][]{
                {true, ClusterDataToSamConverter.PopulateBarcode.ALWAYS, false, new File(testDirParent, "sams_with_BC_ALWAYS")},
                {false, ClusterDataToSamConverter.PopulateBarcode.INEXACT_MATCH, true, new File(testDirParent, "sams_with_INEXACT_MATCH_QUALITY")},
                {true, ClusterDataToSamConverter.PopulateBarcode.ALWAYS, true, new File(testDirParent, "sams_with_BC_ALWAYS_QUALITY")},
                {false, ClusterDataToSamConverter.PopulateBarcode.ORPHANS_ONLY, false, TEST_DATA_DIR}
        };
    }

    @Test(dataProvider = "multiplexedData")
    public void testMultiplexed(final boolean includeBcInHeader, final ClusterDataToSamConverter.PopulateBarcode populateBarcode,
                                final boolean includeBarcodeQuality, final File testDataDir) throws Exception {
        runStandardTest(1, "multiplexedBarcode.", "library.params", 1, "25T8B25T", BASECALLS_DIR, testDataDir, null, includeBcInHeader, populateBarcode, includeBarcodeQuality);
    }

    @DataProvider
    public Object[][] variousConfigurationsData() {
        return new Object[][]{
                {"multiplexedBarcode.", "library.params", 1, "25T8B25T", BASECALLS_DIR, new File(TEST_DATA_DIR.getParentFile(),"sams_with_DS"), null},
                {"multiplexedBarcode.", "library.params", 1, "25T8B25T", BASECALLS_DIR, TEST_DATA_DIR, null},
                {"multiplexedBarcode.", "library.params", 1, "25T8B4M21T", BASECALLS_DIR, TEST_DATA_DIR_WITH_4M_INDEX, null},
                {"multiplexedBarcode2.", "library.params", 1, "25T8B4M4M17T", BASECALLS_DIR, TEST_DATA_DIR_WITH_4M4M_INDEX, null},
                {"singleBarcodeAltName.", "multiplexed_positive_rgtags.params", 1, "25T8B25T", BASECALLS_DIR, TEST_DATA_DIR, null},
                {"dualBarcode.", "library_double.params", 2, "25T8B8B25T", DUAL_BASECALLS_DIR, DUAL_TEST_DATA_DIR, null},
                {"cbclConvert.", "library_double.params", 2, "151T8B8B151T", TEST_DATA_DIR_WITH_CBCLS, DUAL_CBCL_TEST_DATA_DIR, null},
                {"hiseqxSingleLocs.", "library_double.params", 2, "25T8B8B25T", TEST_DATA_HISEQX_SINGLE_LOCS, HISEQX_TEST_DATA_DIR, null},
                {"hiseqxSingleLocs.", "library_double.params", 2, "25T8B8B25T", TEST_DATA_HISEQX_SINGLE_LOCS, HISEQX_TEST_DATA_DIR, null},
                {"dualBarcode.", "library_double.params", 2, "25T8B8B25T", DUAL_BASECALLS_DIR, DUAL_TEST_DATA_DIR, 1101},
                {"cbclConvert.", "library_double.params", 2, "151T8B8B151T", TEST_DATA_DIR_WITH_CBCLS, DUAL_CBCL_TEST_DATA_DIR, 1102}
        };
    }

    @Test(dataProvider = "variousConfigurationsData")
    public void testVariousConfigurations(final String jobName, final String libraryParamsFile, final int nColumnFields, final String cigar, final File baseCallingDir, final File samDir, final Integer tile) throws Exception {
        runStandardTest(1, jobName, libraryParamsFile, nColumnFields, cigar, baseCallingDir, samDir, tile, false, ClusterDataToSamConverter.PopulateBarcode.ORPHANS_ONLY, false);
    }

    /**
     * Ensures that a run missing a barcode from the parameters file throws an error.
     */
    @Test
    public void testCorruptDataReturnCode() throws Exception {
        boolean exceptionThrown = false;
        try {
            runStandardTest(9, "dualBarcode.", "negative_test.params", 2, "30T8B8B", BASECALLS_DIR, TEST_DATA_DIR, null, false, ClusterDataToSamConverter.PopulateBarcode.ORPHANS_ONLY, false);
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
     * @param includeBcInHeader
     * @param populateBarcode
     * @param includeBarcodeQuality @throws Exception
     */
    private void runStandardTest(final int lane, final String jobName, final String libraryParamsFile,
                                 final int concatNColumnFields, final String readStructure,
                                 final File baseCallsDir, final File testDataDir, final Integer tile, final boolean includeBcInHeader, final ClusterDataToSamConverter.PopulateBarcode populateBarcode,
                                 final boolean includeBarcodeQuality) throws Exception {
        for (final boolean sort : new boolean[]{false, true}) {
            final Path outputDir = Files.createTempDirectory(jobName + sort);
            try {
                final String tilePrefix = (tile != null) ? tile + "." : "";

                // Create library.params with output files in the temp directory
                final File libraryParams = new File(outputDir.toFile(), libraryParamsFile);
                libraryParams.deleteOnExit();
                final List<File> samFiles = new ArrayList<>();
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
                args.add("SEQUENCING_CENTER=BI");
                args.add("LIBRARY_PARAMS=" + libraryParams);
                args.add("INCLUDE_BC_IN_RG_TAG=" + includeBcInHeader);
                args.add("BARCODE_POPULATION_STRATEGY=" + populateBarcode.name());
                args.add("INCLUDE_BARCODE_QUALITY=" + includeBarcodeQuality);
                args.add("SORT=" + sort);

                if (tile != null) {
                    args.add("PROCESS_SINGLE_TILE=" + tile);
                }

                Assert.assertEquals(runPicardCommandLine(args), 0);

                for (final File outputSam : samFiles) {
                    if (sort) {
                        IOUtil.assertFilesEqual(outputSam, new File(testDataDir, outputSam.getName()));
                    } else {
                        compareSams(outputSam, new File(testDataDir, outputSam.getName()));
                    }
                }
            } finally {
                IOUtil.recursiveDelete(outputDir);
            }
        }
    }

    private void compareSams(File testSam, File sam) {
        SamReader testReader = SamReaderFactory.makeDefault().open(testSam);
        SamReader samReader = SamReaderFactory.makeDefault().open(sam);
        samReader.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.unsorted);

        SamComparison comparison = new SamComparison(testReader, samReader);
        Assert.assertTrue(comparison.areEqual());
    }
}
