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
package picard.sam;

import htsjdk.io.IOPath;
import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.utils.ValidationUtils;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.http.nio.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramTest;
import picard.nio.PicardBucketUtils;
import picard.nio.PicardHtsPath;
import picard.nio.PicardIOUtils;
import picard.util.GCloudTestUtils;
import picard.util.TabbedInputParser;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Created by IntelliJ IDEA.
 * User: ktibbett
 * Date: Jul 20, 2010
 * Time: 10:27:58 AM
 * To change this template use File | Settings | File Templates.
 */
public class RevertSamTest extends CommandLineProgramTest {
    private final static String DEFAULT_CLOUD_TEST_OUTPUT_RELATIVE_PATH = "test/RevertSam/";
    private final static PicardHtsPath DEFAULT_CLOUD_TEST_OUTPUT_DIR = PicardHtsPath.resolve(GCloudTestUtils.TEST_OUTPUT_DEFAULT_GCLOUD, DEFAULT_CLOUD_TEST_OUTPUT_RELATIVE_PATH);
    public static final String REVERT_SAM_LOCAL_TEST_DATA_DIR = "testdata/picard/sam/RevertSam/";
    private static final String basicSamToRevert = REVERT_SAM_LOCAL_TEST_DATA_DIR + "revert_sam_basic.sam";
    private static final String sampleLibraryOverrideSam = REVERT_SAM_LOCAL_TEST_DATA_DIR + "revert_sam_sample_library_override.sam";
    private static final File validOutputMap = new File(REVERT_SAM_LOCAL_TEST_DATA_DIR + "revert_sam_valid_output_map.txt");
    private static final File nonExistentOutputMap = new File(REVERT_SAM_LOCAL_TEST_DATA_DIR + "revert_sam_does_not_exist.txt");
    private static final File badHeaderOutputMap = new File(REVERT_SAM_LOCAL_TEST_DATA_DIR + "revert_sam_bad_header_output_map.txt");
    private static final File samTestData = new File("testdata/picard/sam");
    private static final File writablePath = new File(REVERT_SAM_LOCAL_TEST_DATA_DIR + "revert_sam_writable.bam");
    private static final File referenceFasta = new File("testdata/picard/reference/test.fasta");
    private static final String singleEndSamToRevert = REVERT_SAM_LOCAL_TEST_DATA_DIR + "revert_sam_single_end.sam";
    private static final File hardClipFasta = new File("testdata/picard/sam/MergeBamAlignment/cliptest.fasta");
    private static final File hardClippedAlignedSam = new File("testdata/picard/sam/MergeBamAlignment/hardclip.aligned.sam");
    private static final File hardClippedUnmappedSam = new File("testdata/picard/sam/MergeBamAlignment/hardclip.unmapped.sam");

    private static final String revertedQualities =
            "11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111";

    private static final String unmappedRead = "both_reads_present_only_first_aligns/2";

    public String getCommandLineProgramName() {
        return RevertSam.class.getSimpleName();
    }

    @Test(dataProvider="positiveTestData")
    public void basicPositiveTests(final SAMFileHeader.SortOrder so, final boolean removeDuplicates, final boolean removeAlignmentInfo, final boolean restoreHardClips,
                                   final boolean restoreOriginalQualities, final boolean outputByReadGroup, final String sample, final String library,
                                   final List<String> attributesToClear) throws Exception {

        final File output;
        File output0 = null;
        File output1 = null;
        File output2 = null;
        if (outputByReadGroup) {
            output = Files.createTempDirectory("picardRevertSamTest").toFile();
            output0 = Paths.get(output.toString(), "0.sam").toFile();
            output1 = Paths.get(output.toString(), "1.sam").toFile();
            output2 = Paths.get(output.toString(), "2.sam").toFile();
        } else {
            output = File.createTempFile("reverted", ".sam");
        }
        output.deleteOnExit();
        final RevertSam reverter = new RevertSam();
        int argSize = 5 + (so != null ? 1 : 0) + attributesToClear.size() + (sample != null ? 1 : 0) + (library != null ? 1 : 0);
        if (outputByReadGroup) {
            argSize++;
        }
        if (!restoreHardClips) {
            argSize++;
        }
        final String args[] = new String[argSize];
        int index = 0;
        args[index++] = "INPUT=" + basicSamToRevert;
        args[index++] = "OUTPUT=" + output.getAbsolutePath();
        if (outputByReadGroup) {
            args[index++] = "OUTPUT_BY_READGROUP=" + outputByReadGroup;
        }
        if (so != null) {
            args[index++] = "SORT_ORDER=" + so.name();
        }
        args[index++] = "REMOVE_DUPLICATE_INFORMATION=" + removeDuplicates;
        args[index++] = "REMOVE_ALIGNMENT_INFORMATION=" + removeAlignmentInfo;
        if (!restoreHardClips) {
            args[index++] = "RESTORE_HARDCLIPS=" + restoreHardClips;
        }
        args[index++] = "RESTORE_ORIGINAL_QUALITIES=" + restoreOriginalQualities;
        if (sample != null) {
            args[index++] = "SAMPLE_ALIAS=" + sample;
        }
        if (library != null) {
            args[index++] = "LIBRARY_NAME=" + library;
        }
        for (final String attr : attributesToClear) {
            args[index++] = "ATTRIBUTE_TO_CLEAR=" + attr;
        }

        runPicardCommandLine(args);

        if (outputByReadGroup) {
            verifyPositiveResults(output0, reverter, removeDuplicates, removeAlignmentInfo, restoreOriginalQualities, outputByReadGroup, "0", 2, sample, library);
            verifyPositiveResults(output1, reverter, removeDuplicates, removeAlignmentInfo, restoreOriginalQualities, outputByReadGroup, "1", 4, sample, library);
            verifyPositiveResults(output2, reverter, removeDuplicates, removeAlignmentInfo, restoreOriginalQualities, outputByReadGroup, "2", 2, sample, library);
        } else {
            verifyPositiveResults(output, reverter, removeDuplicates, removeAlignmentInfo, restoreOriginalQualities, outputByReadGroup, null, 8, sample, library);
        }
    }

    @Test
    public void testOutputByReadGroupWithOutputMap() throws Exception {
        final File outputDir = Files.createTempDirectory("tmpPicardTest").toFile();
        //outputDir.deleteOnExit();
        // Create the output map
        final File outputMapFile = Files.createTempFile("picardRevertSamTestOutputMap", ".txt").toFile();
        final PrintWriter mapWriter = new PrintWriter(outputMapFile);
        final String outputPath0 = outputDir + "/my_rg0.sam";
        final String outputPath1 = outputDir + "/rg1.cram";
        final String outputPath2 = outputDir + "/my_rg2.bam";
        final String outputPath3 = outputDir + "/my_rg3.sam";
        mapWriter.println("READ_GROUP_ID\tOUTPUT");
        mapWriter.println("0\t" + outputPath0);
        mapWriter.println("2\t" + outputPath2);
        mapWriter.println("1\t" + outputPath1);
        mapWriter.println("3\t" + outputPath3);
        System.out.println("outputFile: " + outputPath0);
        System.out.println("outputFile: " + outputPath1);
        System.out.println("outputFile: " + outputPath2);
        System.out.println("outputFile: " + outputPath3);
        mapWriter.close();
        outputMapFile.deleteOnExit();

        final RevertSam reverter = new RevertSam();

        final String args[] = new String[11];
        int index = 0;
        args[index++] = "INPUT=" + basicSamToRevert;
        args[index++] = "OUTPUT_BY_READGROUP=true";
        args[index++] = "OUTPUT_MAP=" + outputMapFile;
        args[index++] = "REFERENCE_SEQUENCE=" + referenceFasta;
        args[index++] = "SORT_ORDER=" + SAMFileHeader.SortOrder.queryname.name();
        args[index++] = "REMOVE_DUPLICATE_INFORMATION=" + true;
        args[index++] = "REMOVE_ALIGNMENT_INFORMATION=" + true;
        args[index++] = "RESTORE_ORIGINAL_QUALITIES=" + true;
        args[index++] = "SAMPLE_ALIAS=" + "test_sample_1";
        args[index++] = "LIBRARY_NAME=" + "test_library_1";
        args[index++] = "ATTRIBUTE_TO_CLEAR=" + SAMTag.NM.name();

        runPicardCommandLine(args);

        final File output0 = new File(outputPath0);
        final File output1 = new File(outputPath1);
        final File output2 = new File(outputPath2);
        verifyPositiveResults(output0, reverter, true, true, true, true, "0", 2, "test_sample_1", "test_library_1");
        verifyPositiveResults(output1, reverter, true, true, true, true, "1", 4, "test_sample_1", "test_library_1");
        verifyPositiveResults(output2, reverter, true, true, true, true, "2", 2, "test_sample_1", "test_library_1");
    }

    @Test
    public void testSingleEnd() throws Exception {
        final File output = File.createTempFile("single_end_reverted", ".sam");
        output.deleteOnExit();
        final String args[] = { "INPUT=" + singleEndSamToRevert, "OUTPUT=" + output.getAbsolutePath() };
        runPicardCommandLine(args);
        final ValidateSamFile validator = new ValidateSamFile();
        validator.INPUT = new PicardHtsPath(output);
        validator.VALIDATION_STRINGENCY = ValidationStringency.STRICT;
        validator.MODE = ValidateSamFile.Mode.VERBOSE;
        final int result = validator.doWork();
        Assert.assertEquals(result, 0, "Validation of reverted single-end sample failed.");
    }

    @Test
    public void testSingleEndSanitize() throws Exception {
        final File output = File.createTempFile("single_end_reverted", ".sam");
        output.deleteOnExit();
        final String args[] = { "INPUT=" + singleEndSamToRevert, "OUTPUT=" + output.getAbsolutePath(), "SANITIZE=true"};
        Assert.assertEquals(runPicardCommandLine(args), 0, "Sanitation of single-end sample failed.");
    }

    private void verifyPositiveResults(
            final File outputFile,
            final RevertSam reverter,
            final boolean removeDuplicates,
            final boolean removeAlignmentInfo,
            final boolean restoreOriginalQualities,
            final boolean outputByReadGroup,
            final String readGroupId,
            final int numReadsExpected,
            final String sample,
            final String library) {

        outputFile.deleteOnExit();
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(referenceFasta).open(outputFile);
        final SAMFileHeader header = reader.getFileHeader();
        Assert.assertEquals(header.getSortOrder(), SAMFileHeader.SortOrder.queryname);
        Assert.assertEquals(header.getProgramRecords().size(), removeAlignmentInfo ? 0 : 1);
        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        if (outputByReadGroup) {
            Assert.assertEquals(readGroups.size(), 1);
            Assert.assertEquals(readGroups.get(0).getId(), readGroupId);
        }
        for (final SAMReadGroupRecord rg : header.getReadGroups()) {
            if (sample != null) {
                Assert.assertEquals(rg.getSample(), sample);
            } else {
                Assert.assertEquals(rg.getSample(), "Hi,Mom!");
            }
            if (library != null) {
                Assert.assertEquals(rg.getLibrary(), library);
            } else {
                Assert.assertEquals(rg.getLibrary(), "my-library");
            }
        }
        int numReads = 0;
        for (final SAMRecord rec : reader) {
            numReads++;
            if (removeDuplicates) {
                Assert.assertFalse(rec.getDuplicateReadFlag(),
                        "Duplicates should have been removed: " + rec.getReadName());
            }

            if (removeAlignmentInfo) {
                Assert.assertTrue(rec.getReadUnmappedFlag(),
                        "Alignment info should have been removed: " + rec.getReadName());
            }

            if (restoreOriginalQualities && !unmappedRead.equals(
                    rec.getReadName() + "/" + (rec.getFirstOfPairFlag() ? "1" : "2"))) {

                Assert.assertEquals(rec.getBaseQualityString(), revertedQualities);
            } else {
                Assert.assertNotSame(rec.getBaseQualityString(), revertedQualities);
            }

            for (final SAMRecord.SAMTagAndValue attr : rec.getAttributes()) {
                if (removeAlignmentInfo || (!attr.tag.equals("PG") && !attr.tag.equals("NM")
                        && !attr.tag.equals(SAMTag.MQ.toString()))) {
                    Assert.assertFalse(reverter.ATTRIBUTE_TO_CLEAR.contains(attr.tag),
                            attr.tag + " should have been cleared.");
                }
            }
        }
        Assert.assertEquals(numReads, numReadsExpected);
        CloserUtil.close(reader);
    }

    @DataProvider(name = "positiveTestData")
    public Object[][] getPostitiveTestData() {
        return new Object[][]{
                {null, true, true, true, true, true, null, null, Collections.EMPTY_LIST},
                {SAMFileHeader.SortOrder.queryname, true, true, true, true, false, "Hey,Dad!", null, Arrays.asList("XT")},
                {null, false, true, true, false, false, "Hey,Dad!", "NewLibraryName", Arrays.asList("XT")},
                {null, false, false, false, false, false, null, null, Collections.EMPTY_LIST}
        };
    }

    @Test(dataProvider = "overrideTestData", expectedExceptions = {PicardException.class})
    public void testSampleLibraryOverride(final String sample, final String library) throws Exception {

        final File output = File.createTempFile("bad", ".sam");
        output.deleteOnExit();
        final String args[] = new String[2 + (sample != null ? 1 : 0) + (library != null ? 1 : 0)];
        int index = 0;
        args[index++] = "INPUT=" + sampleLibraryOverrideSam;
        args[index++] = "OUTPUT=" + output.getAbsolutePath();
        if (sample != null) {
            args[index++] = "SAMPLE_ALIAS=" + sample;
        }
        if (library != null) {
            args[index++] = "LIBRARY_NAME=" + library;
        }
        runPicardCommandLine(args);
        Assert.fail("Negative test should have thrown an exception and didn't");
    }

    @DataProvider(name = "overrideTestData")
    public Object[][] getNegativeTestData() {
        return new Object[][]{
                {"NewSample", null},
                {null, "NewLibrary"},
                {"NewSample", "NewLibrary"}
        };
    }

    @Test
    public void testMutexOutputMapVsOutput() throws Exception {
        final File outputDir = Files.createTempDirectory("picardRevertSamTest").toFile();
        outputDir.deleteOnExit();

        final String args[] = new String[4];
        int index = 0;
        args[index++] = "INPUT=" + basicSamToRevert;
        args[index++] = "OUTPUT_BY_READGROUP=true";
        args[index++] = "OUTPUT=" + outputDir;
        args[index++] = "OUTPUT_MAP=" + validOutputMap;

        try {
            final int returnCode = runPicardCommandLine(args);
            Assert.assertEquals(returnCode, 1);
        } catch (CommandLineException e) {
            // Barclay parser throws on mutex violation
            Assert.assertFalse(CommandLineProgram.useLegacyParser());
        }
    }

    @Test
    public void testNoInput() throws Exception {
        final File outputDir = Files.createTempDirectory("picardRevertSamTest").toFile();
        outputDir.deleteOnExit();

        final String args[] = new String[0];
        try {
            final int returnCode = runPicardCommandLine(args);
            Assert.assertEquals(returnCode, 1);
        } catch (CommandLineException e) {
            // Barclay parser throws on command line errors
            Assert.assertFalse(CommandLineProgram.useLegacyParser());
        }
    }

    @Test
    public void testCommandLineHelp() throws Exception {
        final File outputDir = Files.createTempDirectory("picardRevertSamTest").toFile();
        outputDir.deleteOnExit();

        final String args[] = new String[1];
        args[0] = "--help";
        final int returnCode = runPicardCommandLine(args);
        Assert.assertEquals(returnCode, 1);
    }

    @Test
    public void testValidateOutputParamsByReadGroupMapValid() {
        final List<String> errors = new ArrayList<>();
        RevertSam.ValidationUtil.validateOutputParamsByReadGroup(null, new PicardHtsPath(validOutputMap), errors);
        Assert.assertEquals(errors.size(), 0);
    }

    @Test
    public void testValidateOutputParamsByReadGroupMissingMap() {
        final List<String> errors = new ArrayList<>();
        RevertSam.ValidationUtil.validateOutputParamsByReadGroup(null, new PicardHtsPath(nonExistentOutputMap), errors);
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("Cannot read"), true);
    }

    @Test
    public void testValidateOutputParamsByReadGroupBadHeaderMap() {
        final List<String> errors = new ArrayList<>();
        RevertSam.ValidationUtil.validateOutputParamsByReadGroup(null, new PicardHtsPath(badHeaderOutputMap), errors);
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("Invalid header"), true);
    }

    @Test
    public void testValidateOutputParamsByReadGroupNoMapOrDir() {
        final List<String> errors = new ArrayList<>();
        RevertSam.ValidationUtil.validateOutputParamsByReadGroup(null, null, errors);
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("Must provide either"), true);
    }

    @Test
    public void testValidateOutputParamsByReadGroupDirValid() {
        final List<String> errors = new ArrayList<>();
        RevertSam.ValidationUtil.validateOutputParamsByReadGroup(new PicardHtsPath(samTestData), null, errors);
        Assert.assertEquals(errors.size(), 0);
    }

    @Test
    public void testValidateOutputParamsNotByReadGroupValid() {
        final List<String> errors = new ArrayList<>();
        RevertSam.ValidationUtil.validateOutputParamsNotByReadGroup(new PicardHtsPath(writablePath), null, errors);
        Assert.assertEquals(errors.size(), 0);
    }

    @Test
    public void testValidateOutputParamsNotByReadGroupNoOutput() {
        final List<String> errors = new ArrayList<>();
        RevertSam.ValidationUtil.validateOutputParamsNotByReadGroup(null, null, errors);
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("OUTPUT is required"), true);
    }

    @Test
    public void testValidateOutputParamsNotByReadGroupMap() {
        final List<String> errors = new ArrayList<>();
        RevertSam.ValidationUtil.validateOutputParamsNotByReadGroup(null, new PicardHtsPath(validOutputMap), errors);
        Assert.assertEquals(errors.size(), 2);
        Assert.assertEquals(errors.get(0).contains("Cannot provide OUTPUT_MAP"), true);
        Assert.assertEquals(errors.get(1).contains("OUTPUT is required"), true);
    }

    @Test
    public void testValidateOutputParamsNotByReadGroupDir() {
        final List<String> errors = new ArrayList<>();
        RevertSam.ValidationUtil.validateOutputParamsNotByReadGroup(new PicardHtsPath(samTestData), null, errors);
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("should not be a directory"), true);
    }

    @Test
    public void testAssertAllReadGroupsMappedSuccess() {
        final SAMReadGroupRecord rg1 = new SAMReadGroupRecord("rg1");
        final SAMReadGroupRecord rg2 = new SAMReadGroupRecord("rg2");

        final Map<String, Path> outputMap = new HashMap<>();
        outputMap.put("rg1", Paths.get("/foo/bar/rg1.bam"));
        outputMap.put("rg2", Paths.get("/foo/bar/rg2.bam"));
        RevertSam.ValidationUtil.assertAllReadGroupsMapped(outputMap, Arrays.asList(rg1, rg2));
        RevertSam.ValidationUtil.assertAllReadGroupsMapped(outputMap, Arrays.asList(rg1));
        RevertSam.ValidationUtil.assertAllReadGroupsMapped(outputMap, Arrays.asList(rg2));
    }

    @Test(expectedExceptions = {PicardException.class})
    public void testAssertAllReadGroupsMappedFailure() {
        final SAMReadGroupRecord rg1 = new SAMReadGroupRecord("rg1");
        final SAMReadGroupRecord rg2 = new SAMReadGroupRecord("rg2");
        final SAMReadGroupRecord rg3 = new SAMReadGroupRecord("rg3");

        final Map<String, Path> outputMap = new HashMap<>();
        outputMap.put("rg1", Paths.get("/foo/bar/rg1.bam"));
        outputMap.put("rg2", Paths.get("/foo/bar/rg2.bam"));
        RevertSam.ValidationUtil.assertAllReadGroupsMapped(outputMap, Arrays.asList(rg1, rg2, rg3));
    }

    @Test
    public void testIsOutputMapHeaderValid() {
        boolean isValid = RevertSam.ValidationUtil.isOutputMapHeaderValid(Arrays.asList("READ_GROUP_ID", "OUTPUT"));
        Assert.assertEquals(isValid, true);

        isValid = RevertSam.ValidationUtil.isOutputMapHeaderValid(Arrays.asList("OUTPUT"));
        Assert.assertEquals(isValid, false);

        isValid = RevertSam.ValidationUtil.isOutputMapHeaderValid(Collections.EMPTY_LIST);
        Assert.assertEquals(isValid, false);
    }

    @Test
    public void testFilePathsWithoutMapFile() {
        final SAMReadGroupRecord rg1 = new SAMReadGroupRecord("rg1");
        final SAMReadGroupRecord rg2 = new SAMReadGroupRecord("rg2");

        final Map<String, Path> outputMap = RevertSam.createOutputMapFromReadGroups(Arrays.asList(rg1, rg2), Paths.get("/foo/bar"), ".bam");
        Assert.assertEquals(outputMap.get("rg1"), Paths.get("/foo/bar/rg1.bam"));
        Assert.assertEquals(outputMap.get("rg2"), Paths.get("/foo/bar/rg2.bam"));
    }

    @Test
    public void testFilePathsWithMapFile() {
        final Map<String, Path> outputMap = RevertSam.readOutputMap(validOutputMap.toPath());
        Assert.assertEquals(outputMap.get("rg1"), Paths.get("/path/to/my_rg_1.ubam"));
        Assert.assertEquals(outputMap.get("rg2"), Paths.get("/path/to/my_rg_2.ubam"));
    }

    @Test
    public void testGetDefaultExtension() {
        Assert.assertEquals(RevertSam.getDefaultExtension("this.is.a.sam"), ".sam");
        Assert.assertEquals(RevertSam.getDefaultExtension("this.is.a.cram"), ".cram");
        Assert.assertEquals(RevertSam.getDefaultExtension("this.is.a.bam"), ".bam");
        Assert.assertEquals(RevertSam.getDefaultExtension("foo"), ".bam");
    }

    @Test(expectedExceptions = PicardException.class)
    public void testNoRgInfoOutputByRg() {
        final String[] args = new String[]{
                "I=testdata/picard/sam/bam2fastq/paired/bad/missing-rg-info.sam",
                "OUTPUT_BY_READGROUP=true",
                "O=."
        };
        runPicardCommandLine(args);
    }

    @Test
    public void testNoRgInfoSanitize() throws Exception {
        final File output = File.createTempFile("no-rg-reverted", ".sam");
        output.deleteOnExit();
        final String[] args = new String[]{
                "I=testdata/picard/sam/bam2fastq/paired/bad/missing-rg-info.sam",
                "SANITIZE=true",
                "O=" + output.getAbsolutePath()
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);
        verifyPositiveResults(output, new RevertSam(), true, true, false, false, null, 240, null, null);
    }

    @Test
    public void testSanitizeAndDeduplicateRecords() throws Exception {
        final File input = File.createTempFile("test-input-santize-and-deduplicate-records", ".sam");
        final File output = File.createTempFile("test-output-santize-and-deduplicate-records", ".sam");

        // Create a SAM file that has duplicate records
        final SamReader reader = SamReaderFactory.makeDefault().open(Paths.get(basicSamToRevert));
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), false, input);
        int numDuplicated = 0;
        for (final SAMRecord rec : reader) {
            writer.addAlignment(rec);
            if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) {
                writer.addAlignment(rec);
                numDuplicated++;
            }
        }
        reader.close();
        writer.close();

        // Make sure some records are duplicated
        Assert.assertTrue(numDuplicated > 0);

        final String[] args = new String[]{
                "I=" + input.getAbsolutePath(),
                "SANITIZE=true",
                "KEEP_FIRST_DUPLICATE=true",
                "MAX_DISCARD_FRACTION=1",
                "O=" + output.getAbsolutePath()
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);
        verifyPositiveResults(output, new RevertSam(), true, true, false, false, null, 8, null, null);
    }

    @Test
    public void testHardClipRoundTrip() throws Exception {
        // Runs sam files through MergeBamAlignment using hard clipping on overlapping reads.
        // Tests to ensure that RevertSam can reconstruct the reads and base qualities from reads that have been hard clipped.

        final File outputMBA = File.createTempFile("test-output-hard-clipped-round-trip-mba", ".sam");
        outputMBA.deleteOnExit();
        final String[] mergeBamAlignmentsArgs = new String[]{
            "UNMAPPED_BAM=" + hardClippedUnmappedSam.getAbsolutePath(),
            "ALIGNED_BAM=" + hardClippedAlignedSam.getAbsolutePath(),
            "OUTPUT=" + outputMBA.getAbsolutePath(),
            "REFERENCE_SEQUENCE=" + hardClipFasta.getAbsolutePath(),
            "HARD_CLIP_OVERLAPPING_READS=true"
        };
        Assert.assertEquals(runPicardCommandLine("MergeBamAlignment", mergeBamAlignmentsArgs), 0);

        final File outputRevert = File.createTempFile("test-output-hard-clipped-round-trip-reverted", ".sam");
        outputRevert.deleteOnExit();
        final String[] revertSamArgs = new String[]{
                "I=" + outputMBA.getAbsolutePath(),
                "RESTORE_HARDCLIPS=true",
                "O=" + outputRevert.getAbsolutePath()
        };
        Assert.assertEquals(runPicardCommandLine("RevertSam", revertSamArgs), 0);

        try (SamReader revertedReader = SamReaderFactory.makeDefault().referenceSequence(referenceFasta).open(outputRevert);
             SamReader unmappedReader = SamReaderFactory.makeDefault().referenceSequence(referenceFasta).open(hardClippedUnmappedSam)) {
            final SAMRecordIterator revertedIterator = revertedReader.iterator();
            final SAMRecordIterator unmappedIterator = unmappedReader.iterator();

            while (revertedIterator.hasNext() && unmappedIterator.hasNext()) {
                final SAMRecord reverted = revertedIterator.next();
                final SAMRecord unmapped = unmappedIterator.next();

                Assert.assertEquals(reverted.getReadString(), unmapped.getReadString());
                Assert.assertEquals(SAMUtils.phredToFastq(reverted.getBaseQualities()), SAMUtils.phredToFastq(unmapped.getBaseQualities()));
            }
            if (revertedIterator.hasNext() || unmappedIterator.hasNext()) {
                Assert.fail("Reverted sam file should be identical in length to unmapped sam file in test, but was not.");
            }
        }
    }

    private final static boolean OUTPUT_BY_READ_GROUP = true;
    @DataProvider(name = "cloudTestData")
    public Object[][] getCloudTestData() {
        final PicardHtsPath tempCloudCram = PicardBucketUtils.getTempFilePath(DEFAULT_CLOUD_TEST_OUTPUT_DIR, "RevertSam", ".cram");
        final PicardHtsPath tempLocalCram = PicardBucketUtils.getLocalTempFilePath("RevertSam", ".cram");
        return new Object[][]{
                // Output by read group without the output map, write output bams in the cloud.
                {NA12878_MEDIUM_GCLOUD, PicardBucketUtils.getRandomGCSDirectory(DEFAULT_CLOUD_TEST_OUTPUT_RELATIVE_PATH), OUTPUT_BY_READ_GROUP, null, null },
                // Output by read group using the output map, write output bams in the cloud
                {NA12878_MEDIUM_GCLOUD, null, OUTPUT_BY_READ_GROUP, DEFAULT_CLOUD_TEST_OUTPUT_DIR.getURIString(), null },
                // Output by read group using the local output map, write output bams in the cloud
                {NA12878_MEDIUM_GCLOUD, null, OUTPUT_BY_READ_GROUP, REVERT_SAM_LOCAL_TEST_DATA_DIR, null },
                // Cram input, output a CRAM for each read group
                {NA12878_MEDIUM_CRAM_GCLOUD, null, OUTPUT_BY_READ_GROUP, DEFAULT_CLOUD_TEST_OUTPUT_DIR.getURIString(), HG19_CHR2021 },
                // Cram input in the cloud, single cloud cram output
                {NA12878_MEDIUM_CRAM_GCLOUD, tempCloudCram, !OUTPUT_BY_READ_GROUP, null, HG19_CHR2021 },
                // Cram input in the cloud, single local cram output
                {NA12878_MEDIUM_CRAM_GCLOUD, tempLocalCram, !OUTPUT_BY_READ_GROUP, null, HG19_CHR2021 }};
    }

    /**
     *
     * Write the output map for RevertSam based on the read groups present in the header file of the given SAM file.
     *
     * @param sam        The read names in the header of this sam file will be used as rows
     * @param perReadGroupPathBase This path will be the base name of the SAM files listed in the outputMap table.
     *                   Not to be confused with the location of the table file (the return value of this method).
     * @return The PicardHtsPath to the tsv with each row a pair (read name, output), where output points to a temp file
     * created in the directory specified in outputBase.
     *
     * Also see RevertSam::createOutputMapFromReadGroups(), which does part of what this method does (creates a dictionary rather than two lists).
     */
    private PicardHtsPath getTmpTestReadGroupMapTable(final IOPath sam, final String perReadGroupPathBase) {
        ValidationUtils.validateArg(sam.getExtension().isPresent(), "The variable sam must be a SAM file but is missing an extension: " + sam.getURIString());

        final String samExtension = sam.getExtension().get();
        final PicardHtsPath perReadGroupTableFile = PicardBucketUtils.getTempFilePath(DEFAULT_CLOUD_TEST_OUTPUT_DIR.getURIString(), ".txt");
        final List<String> readGroups = getReadGroups(sam.toPath()).stream().map(rg -> rg.getId()).collect(Collectors.toList());
        final List<String> readGroupOutput = new ArrayList<>();

        for (String readGroup : readGroups){
            final PicardHtsPath readGroupOutputPath = PicardBucketUtils.getTempFilePath(perReadGroupPathBase, readGroup, samExtension);
            readGroupOutput.add(readGroupOutputPath.getURIString());
        }

        writeOutputMapToFile(perReadGroupTableFile.toPath(), readGroups, readGroupOutput);
        return perReadGroupTableFile;
    }

    /**
     *
     * Write an output map i.e. tab-delimited pairs (read name, file name) to outputTablePath
     *
     * @param outputTablePath The Path to the file where the table is to be written.
     * @param readGroupNames The list of ordered read names, to be matched with the output files by index.
     * @param outputFiles The list of ordered names of the output files.
     */
    private static void writeOutputMapToFile(final Path outputTablePath, final List<String> readGroupNames, final List<String> outputFiles) {
        ValidationUtils.validateArg(readGroupNames.size() == outputFiles.size(), "readGroupNames and output files must be the same size but got: " + readGroupNames.size() + " and " + outputFiles.size());
        try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outputTablePath))) {
            // Write header
            writeRow(writer, RevertSam.READ_GROUP_ID_COLUMN_NAME, RevertSam.OUTPUT_COLUMN_NAME);

            for (int i = 0; i < readGroupNames.size(); i++){
                writeRow(writer, readGroupNames.get(i), outputFiles.get(i));
            }
        } catch (IOException e) {
            throw new PicardException("Error writing a tsv file to " + outputTablePath, e);
        }
    }

    // Write to a file (writer) the given array of strings as a tab-delimited row
    private static void writeRow(final PrintWriter writer, final String... row) {
        final int delimiter = '\t';
        for (int j = 0; j < row.length; j++){
            writer.write(row[j]);
            if (j < row.length - 1) {
                writer.write(delimiter);
            } else {
                writer.println();
            }
        }
    }

    private List<SAMReadGroupRecord> getReadGroups(final Path sam){
        final List<SAMReadGroupRecord> readGroupsInInput = SamReaderFactory.makeDefault().open(sam).getFileHeader().getReadGroups();
        return readGroupsInInput;
    }

    /***
     *
     * @param inputSAM The path to the input SAM file
     * @param outputPath The path to a directory, when output by read group is true and the map isn't provided.
     *                   The path to a single output SAM file, when output by read group is turned off. May be null.
     * @param outputByReadGroup Whether to output by read group
     * @param outputMapDir The location for a dynamically created output map file. May be null.
     */
    @Test(dataProvider = "cloudTestData", groups = "cloud")
    public void testCloud(final PicardHtsPath inputSAM, final PicardHtsPath outputPath, final boolean outputByReadGroup,
                          final String outputMapDir, final PicardHtsPath reference) {
        final List<String> args = new ArrayList<>(Arrays.asList(
                "INPUT=" + inputSAM.getURIString(),
                "OUTPUT_BY_READGROUP=" + outputByReadGroup));

        if (outputPath != null){
            args.add("OUTPUT=" + outputPath.getURIString());
        }

        PicardHtsPath outputMap = null;
        if (outputMapDir != null){
            // Dynamically create a temporary output map file based on the read groups in inputSAM.
            outputMap = getTmpTestReadGroupMapTable(inputSAM, outputMapDir);
            args.add("OUTPUT_MAP=" + outputMap.getURIString());
        }

        if (inputSAM.isCram()){
            Utils.nonNull(reference, () -> "A reference must also be supplied with a CRAM input");
            args.add("REFERENCE_SEQUENCE=" + reference.getURIString());
        }

        Assert.assertEquals(runPicardCommandLine(args), 0);

        //                     |Output by Read Group|
        //                    /                    \
        //                  YES                    NO
        //                  /                        \
        //             |OUTPUT_MAP provided?|    |OUTPUT={file}|
        //               /                \       *** Case 1 ***
        //             YES                NO
        //             /                    \
        //     |OUTPUT_MAP={map}|    |OUTPUT={directory}|
        //       *** Case 3 ***         *** Case 2 ***

        // Collect read groups in the input SAM file to compare with the output read groups.
        final List<SAMReadGroupRecord> readGroupsInInput = SamReaderFactory.makeDefault().open(inputSAM.toPath()).getFileHeader().getReadGroups();

        if (!outputByReadGroup){
            // *** Case 1: Single output file ***
            ValidationUtils.validateArg(outputPath != null, "outputPath must be provided and point to a SAM file to be created.");
            final List<SAMReadGroupRecord> readGroupsInOutput = SamReaderFactory.makeDefault().open(outputPath.toPath()).getFileHeader().getReadGroups();
            Assert.assertEquals(readGroupsInOutput, readGroupsInInput);
        } else if (outputMapDir == null){
            // *** Case 2: OUTPUT_BY_READ_GROUP=true, but OUTPUT_MAP was not provided. ***
            ValidationUtils.validateArg(outputPath != null, "outputPath must be provided.");

            // Confirm that the expected output files have been created in the directory specified by outputPath
            final Map<String, Path> expectedOutputPaths = RevertSam.createOutputMapFromReadGroups(readGroupsInInput, outputPath.toPath(), RevertSam.getDefaultExtension(inputSAM.getURIString()));

            for (Map.Entry<String, Path> entry : expectedOutputPaths.entrySet()){
                final Path outputFile = entry.getValue();
                Assert.assertTrue(Files.exists(outputFile));
                final List<SAMReadGroupRecord> readGroups = SamReaderFactory.makeDefault().open(outputFile).getFileHeader().getReadGroups();
                Assert.assertEquals(readGroups.size(), 1);
                Assert.assertEquals(readGroups.get(0).getId(), entry.getKey());

                // We must mark these output files to be deleted manually, because we cannot designate them to be temp files
                if (Files.exists(entry.getValue())){
                    PicardIOUtils.deleteOnExit(outputFile);
                }
            }
        } else {
            // *** Case 3: OUTPUT_BY_READ_GROUP=true and outputMap was provided ***
            try (final TabbedInputParser intermediateParser = new TabbedInputParser(false, Files.newInputStream(outputMap.toPath()));
                 final TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(intermediateParser)){
                final Iterator<TabbedTextFileWithHeaderParser.Row> rowIterator = parser.iterator();
                while (rowIterator.hasNext()){
                    final TabbedTextFileWithHeaderParser.Row row = rowIterator.next();
                    final String readGroup = row.getField(RevertSam.READ_GROUP_ID_COLUMN_NAME);
                    final String outputForReadGroup = row.getField(RevertSam.OUTPUT_COLUMN_NAME);
                    final PicardHtsPath outputForReadGroupPath = new PicardHtsPath(outputForReadGroup);
                    Assert.assertTrue(Files.exists(outputForReadGroupPath.toPath()), "The expected output file does not exist: " +  outputForReadGroup);

                    final List<SAMReadGroupRecord> readGroupsInOutput = SamReaderFactory.makeDefault().open(outputForReadGroupPath.toPath()).getFileHeader().getReadGroups();
                    Assert.assertEquals(readGroupsInOutput.size(), 1);
                    Assert.assertEquals(readGroupsInOutput.get(0).getId(), readGroup);
                }
            } catch (IOException e){
                throw new PicardException("Encountered an exception while parsing the output map", e);
            }
        }
    }
}
