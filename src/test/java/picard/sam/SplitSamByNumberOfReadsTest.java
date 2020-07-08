/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.sam.testers.ValidateSamTester;
import com.google.common.collect.Iterators;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

/**
 * Tests for SAMSplitter
 */
public class SplitSamByNumberOfReadsTest extends CommandLineProgramTest {
    private static final String TEST_DATA_DIR = "testdata/picard/sam/SplitSamByNumberOfReads";
    private static final File PAIRED_FILE = new File(TEST_DATA_DIR + "/sorted-pair.sam");
    private static final File THREE_READ_TEMPLATE = new File(TEST_DATA_DIR + "/sorted-triplet.sam");
    private final ValidateSamTester VALIDATE_SAM_TESTER = new ValidateSamTester();
    private final String TMP_DIR_NAME = "temp_dir";

    public String getCommandLineProgramName() {
        return SplitSamByNumberOfReads.class.getSimpleName();
    }

    @Test
    public void testOkFile() throws IOException {
        final File tmpDir = Files.createTempDirectory(TMP_DIR_NAME).toFile();

        final String[] args = {
                "INPUT=" + PAIRED_FILE.getPath(),
                "TOTAL_READS_IN_INPUT=10",
                "SPLIT_TO_N_READS=5",
                "OUTPUT=" + tmpDir
        };

        runPicardCommandLine(args);

        final File out1 = new File(tmpDir, "shard_0001.bam");
        out1.deleteOnExit();
        final File out2 = new File(tmpDir, "shard_0002.bam");
        out2.deleteOnExit();
        final SamReader input = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(PAIRED_FILE);
        final SAMRecordIterator inputIter = input.iterator();
        final SamReader reader1 = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(out1);
        final SamReader reader2 = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(out2);
        final SAMFileHeader inputHeader = input.getFileHeader();
        final SAMFileHeader reader1Header = reader1.getFileHeader();

        // The output BAM versions might not match the input header version if the input is outdated.
        inputHeader.setAttribute(SAMFileHeader.VERSION_TAG, reader1Header.getVersion());

        Assert.assertEquals(inputHeader, reader1Header);
        Assert.assertEquals(reader1.getFileHeader(), reader2.getFileHeader());
        VALIDATE_SAM_TESTER.assertSamValid(out1);
        VALIDATE_SAM_TESTER.assertSamValid(out2);

        compareInputWithOutputs(reader1, reader2, inputIter, 6);
    }

    @Test
    public void testNoTotalReads() throws IOException {
        final File tmpDir = Files.createTempDirectory(TMP_DIR_NAME).toFile();

        final String[] args = {
                "INPUT=" + PAIRED_FILE.getPath(),
                "SPLIT_TO_N_READS=5",
                "OUTPUT=" + tmpDir
        };

        runPicardCommandLine(args);

        final File out1 = new File(tmpDir, "shard_0001.bam");
        out1.deleteOnExit();
        final File out2 = new File(tmpDir, "shard_0002.bam");
        out2.deleteOnExit();
        final SamReader input = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(PAIRED_FILE);
        final SAMRecordIterator inputIter = input.iterator();
        final SamReader reader1 = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(out1);
        final SamReader reader2 = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(out2);
        final SAMFileHeader inputHeader = input.getFileHeader();
        final SAMFileHeader reader1Header = reader1.getFileHeader();

        // The output BAM versions might not match the input header version if the input is outdated.
        inputHeader.setAttribute(SAMFileHeader.VERSION_TAG, reader1Header.getVersion());

        Assert.assertEquals(inputHeader, reader1Header);
        Assert.assertEquals(reader1.getFileHeader(), reader2.getFileHeader());
        VALIDATE_SAM_TESTER.assertSamValid(out1);
        VALIDATE_SAM_TESTER.assertSamValid(out2);

        compareInputWithOutputs(reader1, reader2, inputIter, 6);
    }

    @Test
    public void testNFiles() throws IOException {
        final File tmpDir = Files.createTempDirectory(TMP_DIR_NAME).toFile();

        final String[] args = {
                "INPUT=" + PAIRED_FILE.getPath(),
                "TOTAL_READS_IN_INPUT=10",
                "SPLIT_TO_N_FILES=2",
                "OUTPUT=" + tmpDir
        };

        runPicardCommandLine(args);

        final File out1 = new File(tmpDir, "shard_0001.bam");
        out1.deleteOnExit();
        final File out2 = new File(tmpDir, "shard_0002.bam");
        out2.deleteOnExit();
        final SAMRecordIterator inputIter = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(PAIRED_FILE).iterator();
        final SamReader reader1 = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(out1);
        final SamReader reader2 = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(out2);

        Assert.assertEquals(reader1.getFileHeader(), reader2.getFileHeader());
        VALIDATE_SAM_TESTER.assertSamValid(out1);
        VALIDATE_SAM_TESTER.assertSamValid(out2);

        compareInputWithOutputs(reader1, reader2, inputIter, 6);
    }

    @Test
    public void testOneOutput() throws IOException {
        final File tmpDir = Files.createTempDirectory(TMP_DIR_NAME).toFile();

        final String[] args = {
                "INPUT=" + PAIRED_FILE.getPath(),
                "TOTAL_READS_IN_INPUT=10",
                "SPLIT_TO_N_READS=10",
                "OUTPUT=" + tmpDir
        };

        runPicardCommandLine(args);

        final File out = new File(tmpDir, "shard_0001.bam");
        out.deleteOnExit();
        final SamReader outReader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(out);
        final SamReader inReader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(PAIRED_FILE);
        VALIDATE_SAM_TESTER.assertSamValid(out);
        Assert.assertEquals(Iterators.size(inReader.iterator()), Iterators.size(outReader.iterator()));
    }

    @Test
    public void testWrongTotalReads() throws IOException {
        final File tmpDir = Files.createTempDirectory(TMP_DIR_NAME).toFile();

        final String[] args = {
                "INPUT=" + PAIRED_FILE.getPath(),
                "TOTAL_READS_IN_INPUT=20",
                "SPLIT_TO_N_READS=5",
                "OUTPUT=" + tmpDir
        };

        runPicardCommandLine(args);

        final File out1 = new File(tmpDir, "shard_0001.bam");
        out1.deleteOnExit();
        final File out2 = new File(tmpDir, "shard_0002.bam");
        out2.deleteOnExit();
        final SAMRecordIterator inputIter = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(PAIRED_FILE).iterator();
        final SamReader reader1 = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(out1);
        final SamReader reader2 = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(out2);

        Assert.assertEquals(reader1.getFileHeader(), reader2.getFileHeader());
        VALIDATE_SAM_TESTER.assertSamValid(out1);
        VALIDATE_SAM_TESTER.assertSamValid(out2);

        compareInputWithOutputs(reader1, reader2, inputIter, 6);
    }

    @Test
    public void testStreamWithoutTotalReads() throws IOException {
        final File tmpDir = Files.createTempDirectory(TMP_DIR_NAME).toFile();

        final String[] args = {
                "INPUT=/dev/stdin",
                "SPLIT_TO_N_READS=5",
                "OUTPUT=" + tmpDir
        };

        final int rc = runPicardCommandLine(args);
        Assert.assertEquals(rc, 1);
    }

    @Test
    public void testThreeReadTemplate() throws IOException {
        final File tmpDir = Files.createTempDirectory(TMP_DIR_NAME).toFile();

        final String[] args = {
                "INPUT=" + THREE_READ_TEMPLATE.getPath(),
                "TOTAL_READS_IN_INPUT=11",
                "SPLIT_TO_N_READS=2",
                "OUTPUT=" + tmpDir
        };

        runPicardCommandLine(args);

        final File out1 = new File(tmpDir, "shard_0001.bam");
        out1.deleteOnExit();
        final File out2 = new File(tmpDir, "shard_0002.bam");
        out2.deleteOnExit();
        final File out3 = new File(tmpDir, "shard_0003.bam");
        out3.deleteOnExit();
        final SAMRecordIterator inputIter = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(THREE_READ_TEMPLATE).iterator();
        final SamReader reader1 = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(out1);
        final SamReader reader2 = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(out2);

        Assert.assertEquals(reader1.getFileHeader(), reader2.getFileHeader());
        VALIDATE_SAM_TESTER.assertSamValid(out1);
        VALIDATE_SAM_TESTER.assertSamValid(out2);

        compareInputWithOutputs(reader1, reader2, inputIter, 3);
    }

    private void compareInputWithOutputs(final SamReader reader1, final SamReader reader2, final SAMRecordIterator inputIterator, final int expectedFirstSize) {
        int count = 0;
        for (SAMRecord rec : reader1) {
            SAMRecord inputRec = inputIterator.next();
            Assert.assertEquals(rec.getReadName(), inputRec.getReadName());
            count++;
        }
        Assert.assertEquals(count, expectedFirstSize);

        for (SAMRecord rec : reader2) {
            SAMRecord inputRec = inputIterator.next();
            Assert.assertEquals(rec.getReadName(), inputRec.getReadName());
        }
    }

    @Test
    public void testOutPrefixWithZeros() throws IOException {
        final File tmpDir = Files.createTempDirectory(TMP_DIR_NAME).toFile();

        final String[] args = {
                "INPUT=" + PAIRED_FILE,
                "SPLIT_TO_N_READS=5",
                "OUTPUT=" + tmpDir,
                "OUT_PREFIX=Sample01"
        };

        runPicardCommandLine(args);
        final File out1 = new File(tmpDir, "Sample01_0001.bam");
        out1.deleteOnExit();
        VALIDATE_SAM_TESTER.assertSamValid(out1);
        final File out2 = new File(tmpDir, "Sample01_0002.bam");
        out2.deleteOnExit();
        VALIDATE_SAM_TESTER.assertSamValid(out2);
    }
}
