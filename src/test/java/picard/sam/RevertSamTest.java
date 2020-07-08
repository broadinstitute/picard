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

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramTest;
import picard.cmdline.PicardCommandLine;

import java.io.File;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: ktibbett
 * Date: Jul 20, 2010
 * Time: 10:27:58 AM
 * To change this template use File | Settings | File Templates.
 */
public class RevertSamTest extends CommandLineProgramTest {
    private static final String basicSamToRevert = "testdata/picard/sam/revert_sam_basic.sam";
    private static final String sampleLibraryOverrideSam = "testdata/picard/sam/revert_sam_sample_library_override.sam";
    private static final File validOutputMap = new File("testdata/picard/sam/revert_sam_valid_output_map.txt");
    private static final File nonExistentOutputMap = new File("testdata/picard/sam/revert_sam_does_not_exist.txt");
    private static final File badHeaderOutputMap = new File("testdata/picard/sam/revert_sam_bad_header_output_map.txt");
    private static final File samTestData = new File("testdata/picard/sam");
    private static final File writablePath = new File("testdata/picard/sam/revert_sam_writable.bam");
    private static final File referenceFasta = new File("testdata/picard/reference/test.fasta");
    private static final String singleEndSamToRevert = "testdata/picard/sam/revert_sam_single_end.sam";
    private static final File hardClipFasta = new File("testdata/picard/sam/MergeBamAlignment/cliptest.fasta");
    private static final File hardClippedAlignedSam = new File("testdata/picard/sam/MergeBamAlignment/hardclip.aligned.sam");
    private static final File hardClippedUnmappedSam = new File("testdata/picard/sam/MergeBamAlignment/hardclip.unmapped.sam");

    private static final String revertedQualities  =
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
        validator.INPUT = output;
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

    @DataProvider(name="positiveTestData")
    public Object[][] getPostitiveTestData() {
        return new Object[][] {
                {null, true, true, true, true, true, null, null, Collections.EMPTY_LIST},
                {SAMFileHeader.SortOrder.queryname, true, true, true, true, false, "Hey,Dad!", null, Arrays.asList("XT")},
                {null, false, true, true, false, false, "Hey,Dad!", "NewLibraryName", Arrays.asList("XT")},
                {null, false, false, false, false, false, null, null, Collections.EMPTY_LIST}
        };
    }

    @Test(dataProvider="overrideTestData", expectedExceptions = {PicardException.class})
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

    @DataProvider(name="overrideTestData")
    public Object[][] getNegativeTestData() {
        return new Object[][] {
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
        final List<String> errors = new ArrayList<String>();
        RevertSam.ValidationUtil.validateOutputParamsByReadGroup(null, validOutputMap, errors);
        Assert.assertEquals(errors.size(), 0);
    }

    @Test
    public void testValidateOutputParamsByReadGroupMissingMap() {
        final List<String> errors = new ArrayList<String>();
        RevertSam.ValidationUtil.validateOutputParamsByReadGroup(null, nonExistentOutputMap, errors);
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("Cannot read"), true);
    }

    @Test
    public void testValidateOutputParamsByReadGroupBadHeaderMap() {
        final List<String> errors = new ArrayList<String>();
        RevertSam.ValidationUtil.validateOutputParamsByReadGroup(null, badHeaderOutputMap, errors);
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("Invalid header"), true);
    }

    @Test
    public void testValidateOutputParamsByReadGroupNoMapOrDir() {
        final List<String> errors = new ArrayList<String>();
        RevertSam.ValidationUtil.validateOutputParamsByReadGroup(null, null, errors);
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("Must provide either"), true);
    }

    @Test
    public void testValidateOutputParamsByReadGroupDirValid() {
        final List<String> errors = new ArrayList<String>();
        RevertSam.ValidationUtil.validateOutputParamsByReadGroup(samTestData, null, errors);
        Assert.assertEquals(errors.size(), 0);
    }

    @Test
    public void testValidateOutputParamsNotByReadGroupValid() {
        final List<String> errors = new ArrayList<String>();
        RevertSam.ValidationUtil.validateOutputParamsNotByReadGroup(writablePath, null, errors);
        Assert.assertEquals(errors.size(), 0);
    }

    @Test
    public void testValidateOutputParamsNotByReadGroupNoOutput() {
        final List<String> errors = new ArrayList<String>();
        RevertSam.ValidationUtil.validateOutputParamsNotByReadGroup(null, null, errors);
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("OUTPUT is required"), true);
    }

    @Test
    public void testValidateOutputParamsNotByReadGroupMap() {
        final List<String> errors = new ArrayList<String>();
        RevertSam.ValidationUtil.validateOutputParamsNotByReadGroup(null, validOutputMap, errors);
        Assert.assertEquals(errors.size(), 2);
        Assert.assertEquals(errors.get(0).contains("Cannot provide OUTPUT_MAP"), true);
        Assert.assertEquals(errors.get(1).contains("OUTPUT is required"), true);
    }

    @Test
    public void testValidateOutputParamsNotByReadGroupDir() {
        final List<String> errors = new ArrayList<String>();
        RevertSam.ValidationUtil.validateOutputParamsNotByReadGroup(samTestData, null, errors);
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("should not be a directory"), true);
    }

    @Test
    public void testAssertAllReadGroupsMappedSuccess() {
        final SAMReadGroupRecord rg1 = new SAMReadGroupRecord("rg1");
        final SAMReadGroupRecord rg2 = new SAMReadGroupRecord("rg2");

        final Map<String, File> outputMap = new HashMap<String, File>();
        outputMap.put("rg1", new File("/foo/bar/rg1.bam"));
        outputMap.put("rg2", new File("/foo/bar/rg2.bam"));
        RevertSam.ValidationUtil.assertAllReadGroupsMapped(outputMap, Arrays.asList(rg1, rg2));
        RevertSam.ValidationUtil.assertAllReadGroupsMapped(outputMap, Arrays.asList(rg1));
        RevertSam.ValidationUtil.assertAllReadGroupsMapped(outputMap, Arrays.asList(rg2));
    }

    @Test(expectedExceptions = {PicardException.class})
    public void testAssertAllReadGroupsMappedFailure() {
        final SAMReadGroupRecord rg1 = new SAMReadGroupRecord("rg1");
        final SAMReadGroupRecord rg2 = new SAMReadGroupRecord("rg2");
        final SAMReadGroupRecord rg3 = new SAMReadGroupRecord("rg3");

        final Map<String, File> outputMap = new HashMap<String, File>();
        outputMap.put("rg1", new File("/foo/bar/rg1.bam"));
        outputMap.put("rg2", new File("/foo/bar/rg2.bam"));
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

        final Map<String, File> outputMap = RevertSam.createOutputMap(null, new File("/foo/bar"), ".bam", Arrays.asList(rg1, rg2));
        Assert.assertEquals(outputMap.get("rg1"), new File("/foo/bar/rg1.bam"));
        Assert.assertEquals(outputMap.get("rg2"), new File("/foo/bar/rg2.bam"));
    }

    @Test
    public void testFilePathsWithMapFile() {
        final Map<String, File> outputMap = RevertSam.createOutputMap(validOutputMap, null, ".bam", Collections.emptyList());
        Assert.assertEquals(outputMap.get("rg1"), new File("/path/to/my_rg_1.ubam"));
        Assert.assertEquals(outputMap.get("rg2"), new File("/path/to/my_rg_2.ubam"));
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
        final String [] args = new String[]{
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
        final String [] args = new String[]{
                "I=testdata/picard/sam/bam2fastq/paired/bad/missing-rg-info.sam",
                "SANITIZE=true",
                "O=" + output.getAbsolutePath()
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);
        verifyPositiveResults(output, new RevertSam(), true, true, false, false, null, 240, null, null);
    }

    @Test
    public void testSanitizeAndDeduplicateRecords() throws Exception {
        final File input  = File.createTempFile("test-input-santize-and-deduplicate-records", ".sam");
        final File output = File.createTempFile("test-output-santize-and-deduplicate-records", ".sam");

        // Create a SAM file that has duplicate records
        final SamReader reader     = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(Paths.get(basicSamToRevert));
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

        final String [] args = new String[]{
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
        final String [] mergeBamAlignmentsArgs = new String[] {
           "UNMAPPED_BAM=" + hardClippedUnmappedSam.getAbsolutePath(),
           "ALIGNED_BAM=" + hardClippedAlignedSam.getAbsolutePath(),
           "OUTPUT=" + outputMBA.getAbsolutePath(),
           "REFERENCE_SEQUENCE=" + hardClipFasta.getAbsolutePath(),
           "HARD_CLIP_OVERLAPPING_READS=true"
        };
        Assert.assertEquals(runPicardCommandLine("MergeBamAlignment", mergeBamAlignmentsArgs), 0);

        final File outputRevert = File.createTempFile("test-output-hard-clipped-round-trip-reverted", ".sam");
        outputRevert.deleteOnExit();
        final String [] revertSamArgs = new String[] {
                "I=" + outputMBA.getAbsolutePath(),
                "RESTORE_HARDCLIPS=true",
                "O=" + outputRevert.getAbsolutePath()
        };
        Assert.assertEquals(runPicardCommandLine("RevertSam", revertSamArgs), 0);

        try(SamReader revertedReader = SamReaderFactory.makeDefault().referenceSequence(referenceFasta).open(outputRevert);
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
}
