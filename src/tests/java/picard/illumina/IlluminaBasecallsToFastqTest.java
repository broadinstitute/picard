/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
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
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.illumina.parser.ReadStructure;

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class IlluminaBasecallsToFastqTest extends CommandLineProgramTest {

    private static final File BASECALLS_DIR = new File("testdata/picard/illumina/25T8B25T/Data/Intensities/BaseCalls");
    private static final File DUAL_BASECALLS_DIR = new File("testdata/picard/illumina/25T8B8B25T/Data/Intensities/BaseCalls");
    private static final File TEST_DATA_DIR = new File("testdata/picard/illumina/25T8B25T/fastq");
    private static final File DUAL_TEST_DATA_DIR = new File("testdata/picard/illumina/25T8B8B25T/fastq");

    public String getCommandLineProgramName() {
        return IlluminaBasecallsToFastq.class.getSimpleName();
    }

    @Test
    public void testNonBarcoded() throws Exception {
        final String suffix = ".1.fastq";
        final File outputFastq1 = File.createTempFile("nonBarcoded.", suffix);
        outputFastq1.deleteOnExit();
        final String outputPrefix = outputFastq1.getAbsolutePath().substring(0, outputFastq1.getAbsolutePath().length() - suffix.length());
        final File outputFastq2 = new File(outputPrefix + ".2.fastq");
        outputFastq2.deleteOnExit();
        final int lane = 1;
        runPicardCommandLine(new String[]{
                "BASECALLS_DIR=" + BASECALLS_DIR,
                "LANE=" + lane,
                "READ_STRUCTURE=25T8B25T",
                "OUTPUT_PREFIX=" + outputPrefix,
                "RUN_BARCODE=HiMom",
                "MACHINE_NAME=machine1",
                "FLOWCELL_BARCODE=abcdeACXX"
        });
        IOUtil.assertFilesEqual(outputFastq1, new File(TEST_DATA_DIR, "nonBarcoded.1.fastq"));
        IOUtil.assertFilesEqual(outputFastq2, new File(TEST_DATA_DIR, "nonBarcoded.2.fastq"));
    }

    @Test
    public void testMultiplexWithIlluminaReadNameHeaders() throws Exception {
        final File outputDir = File.createTempFile("testMultiplexRH.", ".dir");
        try {
            outputDir.delete();
            outputDir.mkdir();
            outputDir.deleteOnExit();

            final String filePrefix = "testMultiplexRH";
            final File outputPrefix = new File(outputDir, filePrefix);

            runPicardCommandLine(new String[]{
                    "BASECALLS_DIR=" + BASECALLS_DIR,
                    "LANE=" + 1,
                    "RUN_BARCODE=HiMom",
                    "READ_STRUCTURE=" + "25T8B25T",
                    "OUTPUT_PREFIX=" + outputPrefix.getAbsolutePath(),
                    "MACHINE_NAME=machine1",
                    "FLOWCELL_BARCODE=abcdeACXX",
                    "READ_NAME_FORMAT=" + IlluminaBasecallsToFastq.ReadNameFormat.ILLUMINA
            });

            final String[] filenames = new String[]{
                filePrefix + ".1.fastq",
                filePrefix + ".barcode_1.fastq"
            };
            for (final String filename : filenames) {
                IOUtil.assertFilesEqual(new File(outputDir, filename), new File(TEST_DATA_DIR, filename));
            }

        } finally {
            TestUtil.recursiveDelete(outputDir);
        }
    }

    @Test
    public void testDeMultiplexed() throws Exception {
        runStandardTest(1, "multiplexedBarcode.", "mp_barcode.params", 1, "25T8B25T", BASECALLS_DIR, TEST_DATA_DIR);
    }

    @Test
    public void testDualBarcodes() throws Exception {
        runStandardTest(1, "dualBarcode.", "barcode_double.params", 2, "25T8B8B25T", DUAL_BASECALLS_DIR, DUAL_TEST_DATA_DIR);
    }

    /**
     * This test utility takes a libraryParamsFile and generates output sam files through IlluminaBasecallsToFastq to compare against
     * preloaded test data
     *
     * @param jobName
     * @param libraryParamsFile
     * @param concatNColumnFields
     * @param readStructureString
     * @throws Exception
     */
    private void runStandardTest(final int lane, final String jobName, final String libraryParamsFile,
                                 final int concatNColumnFields, final String readStructureString, final File baseCallsDir,
                                 final File testDataDir) throws Exception {
        final File outputDir = File.createTempFile(jobName, ".dir");
        try {
            outputDir.delete();
            outputDir.mkdir();
            outputDir.deleteOnExit();
            // Create barcode.params with output files in the temp directory
            final File libraryParams = new File(outputDir, libraryParamsFile);
            libraryParams.deleteOnExit();
            final List<File> outputPrefixes = new ArrayList<File>();
            final LineReader reader = new BufferedLineReader(new FileInputStream(new File(testDataDir, libraryParamsFile)));
            final PrintWriter writer = new PrintWriter(libraryParams);
            final String header = reader.readLine();
            writer.println(header + "\tOUTPUT_PREFIX");
            while (true) {
                final String line = reader.readLine();
                if (line == null) {
                    break;
                }
                final String[] fields = line.split("\t");
                final File outputPrefix = new File(outputDir, StringUtil.join("", Arrays.copyOfRange(fields, 0, concatNColumnFields)));
                outputPrefixes.add(outputPrefix);
                writer.println(line + "\t" + outputPrefix);
            }
            writer.close();
            reader.close();

            runPicardCommandLine(new String[]{
                    "BASECALLS_DIR=" + baseCallsDir,
                    "LANE=" + lane,
                    "RUN_BARCODE=HiMom",
                    "READ_STRUCTURE=" + readStructureString,
                    "MULTIPLEX_PARAMS=" + libraryParams,
                    "MACHINE_NAME=machine1",
                    "FLOWCELL_BARCODE=abcdeACXX"
            });

            final ReadStructure readStructure = new ReadStructure(readStructureString);
            for (final File outputSam : outputPrefixes) {
                for (int i = 1; i <= readStructure.templates.length(); ++i) {
                    final String filename = outputSam.getName() + "." + i + ".fastq";
                    IOUtil.assertFilesEqual(new File(outputSam.getParentFile(), filename), new File(testDataDir, filename));
                }
                for (int i = 1; i <= readStructure.barcodes.length(); ++i) {
                    final String filename = outputSam.getName() + ".barcode_" + i + ".fastq";
                    IOUtil.assertFilesEqual(new File(outputSam.getParentFile(), filename), new File(testDataDir, filename));
                }
            }
        } finally {
            TestUtil.recursiveDelete(outputDir);
        }
    }
}
