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
package net.sf.picard.illumina;

import net.sf.picard.illumina.parser.ReadStructure;
import net.sf.picard.io.IoUtil;
import net.sf.samtools.util.BufferedLineReader;
import net.sf.samtools.util.LineReader;
import net.sf.samtools.util.StringUtil;
import net.sf.samtools.util.TestUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class IlluminaBasecallsToFastqTest {

    private static final File BASECALLS_DIR = new File("testdata/net/sf/picard/illumina/IlluminaTests/BasecallsDir");
    private static final File TEST_DATA_DIR = new File("testdata/net/sf/picard/illumina/IlluminaBasecallsToFastqTest");

    @Test
    public void testNonBarcoded() throws Exception {
        final String suffix = ".1.fastq";
        final File outputFastq1 = File.createTempFile("nonBarcoded.", suffix);
        outputFastq1.deleteOnExit();
        final String outputPrefix = outputFastq1.getAbsolutePath().substring(0, outputFastq1.getAbsolutePath().length() - suffix.length());
        final File outputFastq2 = new File(outputPrefix + ".2.fastq");
        outputFastq2.deleteOnExit();
        final int lane = 1;
        new IlluminaBasecallsToFastq().instanceMain(new String[]{
                "BASECALLS_DIR=" + BASECALLS_DIR,
                "LANE=" + lane,
                "READ_STRUCTURE=76T76T",
                "OUTPUT_PREFIX=" + outputPrefix,
                "RUN_BARCODE=HiMom"
        });
        IoUtil.assertFilesEqual(outputFastq1, new File(TEST_DATA_DIR, "nonBarcoded.1.fastq"));
        IoUtil.assertFilesEqual(outputFastq2, new File(TEST_DATA_DIR, "nonBarcoded.2.fastq"));
    }

    @Test
    public void testMultiplex() throws Exception {
        final File outputDir = File.createTempFile("testMultiplex.", ".dir");
        try {
            outputDir.delete();
            outputDir.mkdir();
            outputDir.deleteOnExit();

            final String filePrefix = "testMultiplex";
            final File outputPrefix = new File(outputDir, filePrefix);

            new IlluminaBasecallsToFastq().instanceMain(new String[]{
                    "BASECALLS_DIR=" + BASECALLS_DIR,
                    "LANE=" + 7,
                    "RUN_BARCODE=HiMom",
                    "READ_STRUCTURE=" + "30T8B",
                    "OUTPUT_PREFIX=" + outputPrefix.getAbsolutePath()
            });

            final String[] filenames = new String[]{
                filePrefix + ".1.fastq",
                filePrefix + ".barcode_1.fastq"
            };
            for (String filename : filenames) {
                IoUtil.assertFilesEqual(new File(outputDir, filename), new File(TEST_DATA_DIR, filename));
            }

        } finally {
            TestUtil.recursiveDelete(outputDir);
        }
    }

    @Test
    public void testDeMultiplexed() throws Exception {
        runStandardTest(7, "multiplexedBarcode.", "barcode.params", 1, "30T8B");
    }

    @Test
    public void testDualBarcodes() throws Exception {
        runStandardTest(9, "dualBarcode.", "barcode_double.params", 2, "30T8B8B");
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
                                 final int concatNColumnFields, final String readStructureString) throws Exception {
        final File outputDir = File.createTempFile(jobName, ".dir");
        try {
            outputDir.delete();
            outputDir.mkdir();
            outputDir.deleteOnExit();
            // Create barcode.params with output files in the temp directory
            final File libraryParams = new File(outputDir, libraryParamsFile);
            libraryParams.deleteOnExit();
            final List<File> outputPrefixes = new ArrayList<File>();
            final LineReader reader = new BufferedLineReader(new FileInputStream(new File(TEST_DATA_DIR, libraryParamsFile)));
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

            new IlluminaBasecallsToFastq().instanceMain(new String[]{
                    "BASECALLS_DIR=" + BASECALLS_DIR,
                    "LANE=" + lane,
                    "RUN_BARCODE=HiMom",
                    "READ_STRUCTURE=" + readStructureString,
                    "MULTIPLEX_PARAMS=" + libraryParams
            });

            final ReadStructure readStructure = new ReadStructure(readStructureString);
            for (final File outputSam : outputPrefixes) {
                for (int i = 1; i <= readStructure.templates.length(); ++i) {
                    String filename = outputSam.getName() + "." + i + ".fastq";
                    IoUtil.assertFilesEqual(new File(outputSam.getParentFile(), filename), new File(TEST_DATA_DIR, filename));
                }
                for (int i = 1; i <= readStructure.barcodes.length(); ++i) {
                    String filename = outputSam.getName() + ".barcode_" + i + ".fastq";
                    IoUtil.assertFilesEqual(new File(outputSam.getParentFile(), filename), new File(TEST_DATA_DIR, filename));
                }
            }
        } finally {
            TestUtil.recursiveDelete(outputDir);
        }
    }
}
