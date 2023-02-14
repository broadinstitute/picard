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

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.illumina.parser.ReadStructure;
import picard.util.IlluminaUtil;

import java.io.*;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Stream;

public class IlluminaBasecallsToFastqTest extends CommandLineProgramTest {

    private static final File BASECALLS_DIR = new File("testdata/picard/illumina/25T8B25T/Data/Intensities/BaseCalls");
    private static final File DUAL_BASECALLS_DIR = new File("testdata/picard/illumina/25T8B8B25T/Data/Intensities/BaseCalls");
    private static final File TEST_DATA_DIR = new File("testdata/picard/illumina/25T8B25T/fastq");
    private static final File TEST_DATA_DIR_WITH_4M = new File("testdata/picard/illumina/25T8B25T/fastq_with_4M");
    private static final File TEST_DATA_DIR_WITH_4M4M = new File("testdata/picard/illumina/25T8B25T/fastq_with_4M4M");
    private static final File TEST_DATA_DIR_WITH_CBCLS = new File("testdata/picard/illumina/151T8B8B151T_cbcl/Data/Intensities/BaseCalls");
    private static final File TEST_DATA_HISEQX_SINGLE_LOCS = new File("testdata/picard/illumina/25T8B8B25T_hiseqx/Data/Intensities/BaseCalls");

    private static final File DUAL_TEST_DATA_DIR = new File("testdata/picard/illumina/25T8B8B25T/fastq");
    private static final File DUAL_CBCL_TEST_DATA_DIR = new File("testdata/picard/illumina/151T8B8B151T_cbcl/fastq");
    private static final File HISEQX_TEST_DATA_DIR = new File("testdata/picard/illumina/25T8B8B25T_hiseqx/fastq");

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
                "FLOWCELL_BARCODE=abcdeACXX",
                "MAX_RECORDS_IN_RAM=100" //force spill to disk to test encode/decode
        });
        IOUtil.assertFilesEqual(outputFastq1, new File(TEST_DATA_DIR, "nonBarcoded.1.fastq"));
        IOUtil.assertFilesEqual(outputFastq2, new File(TEST_DATA_DIR, "nonBarcoded.2.fastq"));
    }

    @Test
    public void testAdapterTrimming() throws Exception {
        final String suffix = ".1.fastq";
        final File outputFastq1 = File.createTempFile("adapterTrimmed.", suffix);
        outputFastq1.deleteOnExit();
        final String outputPrefix = outputFastq1.getAbsolutePath().substring(0, outputFastq1.getAbsolutePath().length() - suffix.length());
        final File outputFastq2 = new File(outputPrefix + ".2.fastq");
        outputFastq2.deleteOnExit();
        final int lane = 1;
        List<String> args = new ArrayList<>(CollectionUtil.makeList(
                "BASECALLS_DIR=" + BASECALLS_DIR,
                "LANE=" + lane,
                "READ_STRUCTURE=25T8B25T",
                "OUTPUT_PREFIX=" + outputPrefix,
                "RUN_BARCODE=HiMom",
                "MACHINE_NAME=machine1",
                "FLOWCELL_BARCODE=abcdeACXX",
                "MAX_RECORDS_IN_RAM=100" //force spill to disk to test encode/decode,
        ));
        for (final IlluminaUtil.IlluminaAdapterPair adapterPair : IlluminaUtil.IlluminaAdapterPair.values()) {
            args.add("ADAPTERS_TO_CHECK=" + adapterPair);
        }
        runPicardCommandLine(args);
        IOUtil.assertFilesEqual(outputFastq1, new File(TEST_DATA_DIR, "nonBarcoded.1.fastq"));
        IOUtil.assertFilesEqual(outputFastq2, new File(TEST_DATA_DIR, "nonBarcoded.2.fastq"));
    }

    @Test
    public void testCustomAdapterTrimming() throws Exception {
        final String suffix = ".1.fastq";
        final File outputFastq1 = File.createTempFile("adapterTrimmed.", suffix);
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
                "FLOWCELL_BARCODE=abcdeACXX",
                "MAX_RECORDS_IN_RAM=100", //force spill to disk to test encode/decode,
                "FIVE_PRIME_ADAPTER=CGGCC",
                "THREE_PRIME_ADAPTER=TGGGC"
        });
        IOUtil.assertFilesEqual(outputFastq1, new File(TEST_DATA_DIR, "adapterTrimmed.1.fastq"));
        IOUtil.assertFilesEqual(outputFastq2, new File(TEST_DATA_DIR, "adapterTrimmed.2.fastq"));
    }

    @Test
    public void testQualityTrimming() throws Exception {
        final String suffix = ".1.fastq";
        final File outputFastq1 = File.createTempFile("qualityTrimmed.", suffix);
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
                "FLOWCELL_BARCODE=abcdeACXX",
                "MAX_RECORDS_IN_RAM=100", //force spill to disk to test encode/decode,
                "TRIMMING_QUALITY=10"
        });
        IOUtil.assertFilesEqual(outputFastq1, new File(TEST_DATA_DIR, "qualityTrimmed.1.fastq"));
        IOUtil.assertFilesEqual(outputFastq2, new File(TEST_DATA_DIR, "qualityTrimmed.2.fastq"));
    }

    @Test
    public void testQualityAndAdapterTrimming() throws Exception {
        final String suffix = ".1.fastq";
        final File outputFastq1 = File.createTempFile("qualityAdapterTrimmed.", suffix);
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
                "FLOWCELL_BARCODE=abcdeACXX",
                "MAX_RECORDS_IN_RAM=100", //force spill to disk to test encode/decode,
                "TRIMMING_QUALITY=10",
                "FIVE_PRIME_ADAPTER=CGGCC",
                "THREE_PRIME_ADAPTER=TGGGC"
        });
        IOUtil.assertFilesEqual(outputFastq1, new File(TEST_DATA_DIR, "qualityAdapterTrimmed.1.fastq"));
        IOUtil.assertFilesEqual(outputFastq2, new File(TEST_DATA_DIR, "qualityAdapterTrimmed.2.fastq"));
    }

    @Test
    public void testMultiplexWithIlluminaReadNameHeaders() throws Exception {
        final File outputDir = Files.createTempDirectory("testMultiplexRH.dir").toFile();
        try {
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
                    "READ_NAME_FORMAT=" + IlluminaBasecallsToFastq.ReadNameFormat.ILLUMINA,
                    "MAX_RECORDS_IN_RAM=100" //force spill to disk to test encode/decode
            });

            final String[] filenames = new String[]{
                    filePrefix + ".1.fastq",
                    filePrefix + ".barcode_1.fastq"
            };
            for (final String filename : filenames) {
                IOUtil.assertFilesEqual(new File(outputDir, filename), new File(TEST_DATA_DIR, filename));
            }

        } finally {
            IOUtil.recursiveDelete(outputDir.toPath());
        }
    }

    @Test
    public void testDeMultiplexed() throws Exception {
        runStandardTest(new int[]{1}, "multiplexedBarcode.", "mp_barcode.params", 1, "25T8B25T", BASECALLS_DIR, TEST_DATA_DIR, 7, 0.04);
    }

    @Test
    public void testDeMultiplexedWithIndex() throws Exception {
        runStandardTest(new int[]{1}, "multiplexedBarcodeWithIndex.", "mp_barcode.params", 1, "25T8B4M21T", BASECALLS_DIR, TEST_DATA_DIR_WITH_4M, 7, 0.04);
    }

    @Test
    public void testDeMultiplexedWithtwoIndexes() throws Exception {
        runStandardTest(new int[]{1}, "multiplexedBarcodeWithTwoIndexes.", "mp_barcode.params", 1, "25T8B4M4M17T", BASECALLS_DIR, TEST_DATA_DIR_WITH_4M4M, 7, 0.04);
    }

    @Test
    public void testDualBarcodes() throws Exception {
        runStandardTest(new int[]{1}, "dualBarcode.", "barcode_double.params", 2, "25T8B8B25T", DUAL_BASECALLS_DIR, DUAL_TEST_DATA_DIR, 2, 0.033333);
    }

    @Test
    public void testCbclConvert() throws Exception {
        runStandardTest(new int[]{1}, "dualBarcode.", "barcode_double.params", 2, "151T8B8B151T", TEST_DATA_DIR_WITH_CBCLS, DUAL_CBCL_TEST_DATA_DIR, 1, 0.25);
    }

    @Test
    public void testHiseqxSingleLocs() throws Exception {
        runStandardTest(new int[]{1}, "hiseqxSingleLocs.", "barcode_double.params", 2, "25T8B8B25T",TEST_DATA_HISEQX_SINGLE_LOCS, HISEQX_TEST_DATA_DIR, 4, 0.033333);
    }

    @Test
    public void testMultipleLanes() throws Exception {
        runStandardTest(new int[]{1, 2}, "dualBarcode.", "barcode_double.params", 2, "25T8B8B25T", DUAL_BASECALLS_DIR, DUAL_TEST_DATA_DIR, 4, 0.033333);
    }

    private void compareFastqs(File actual, File expected) {
        List<FastqRecord> actualReads = slurpReads(actual);
        List<FastqRecord> expectedReads = slurpReads(expected);

        actualReads.sort(Comparator.comparing(FastqRecord::getReadName));
        expectedReads.sort(Comparator.comparing(FastqRecord::getReadName));
        Assert.assertEquals(actualReads, expectedReads);
    }

    private List<FastqRecord> slurpReads(final File f) {
        final FastqReader in = new FastqReader(f);
        final List<FastqRecord> recs = new ArrayList<>();
        for (final FastqRecord r : in) {
          recs.add(r);
        }
        in.close();
        return recs;
    }

    /**
     * This test utility takes a libraryParamsFile and generates output sam files through IlluminaBasecallsToFastq to compare against
     * preloaded test data
     *
     * @param lanes               lanes number to use
     * @param jobName             name of job for the temp file
     * @param libraryParamsFile   the params file to use for the de-multiplexing
     * @param concatNColumnFields how many columns to concatenate to get the barcode
     * @param readStructureString what read-structure string to use
     * @param baseCallsDir        what directory can I find the BCLs in
     * @param testDataDir         what directory can I find the expected resulting files
     * @throws Exception
     */

    private void runStandardTest(final int[] lanes, final String jobName, final String libraryParamsFile,
                                 final int concatNColumnFields, final String readStructureString, final File baseCallsDir,
                                 final File testDataDir, final long expectedPfMatches, final Double expectedPctMatches) throws Exception {
        final File outputDir = Files.createTempDirectory(jobName + ".dir").toFile();
        try {

            outputDir.deleteOnExit();
            // Create barcode.params with output files in the temp directory
            final File libraryParams = new File(outputDir, libraryParamsFile);
            libraryParams.deleteOnExit();
            final List<File> outputPrefixes = new ArrayList<>();
            convertParamsFile(libraryParamsFile, concatNColumnFields, testDataDir, outputDir, libraryParams, outputPrefixes, lanes.length > 1);
            String[] laneArgs = new String[lanes.length];
            for(int i = 0; i < lanes.length; i++) {
                laneArgs[i] = "LANE=" + lanes[i];
            }

            for (final boolean sort : new boolean[]{ true, false} ) {
                for (final boolean otfBarcodeExtract : new boolean[]{ false, true }) {
                    final File metricsFile = File.createTempFile("ibtf.", ".metrics");
                    metricsFile.deleteOnExit();
                    String[] args = {
                            "BASECALLS_DIR=" + baseCallsDir,
                            "RUN_BARCODE=HiMom",
                            "READ_STRUCTURE=" + readStructureString,
                            "MULTIPLEX_PARAMS=" + libraryParams,
                            "MACHINE_NAME=machine1",
                            "FLOWCELL_BARCODE=abcdeACXX",
                            "MAX_RECORDS_IN_RAM=1000", //force spill to disk to test encode/decode,
                            "SORT=" + sort,
                            "MATCH_BARCODES_INLINE=" + otfBarcodeExtract,
                            "METRICS_FILE=" + metricsFile
                    };
                    String[] allArgs = Stream.concat(Arrays.stream(args), Arrays.stream(laneArgs))
                            .toArray(String[]::new);
                    runPicardCommandLine(allArgs);

                    final ReadStructure readStructure = new ReadStructure(readStructureString);
                    for (final File prefix : outputPrefixes) {
                        String filePrefix = prefix.getName();
                        for (int i = 1; i <= readStructure.templates.length(); ++i) {
                            final String filename = filePrefix + "." + i + ".fastq";
                            compareResults(testDataDir, sort, prefix, filename);
                        }
                        for (int i = 1; i <= readStructure.sampleBarcodes.length(); ++i) {
                            final String filename = filePrefix + ".barcode_" + i + ".fastq";
                            compareResults(testDataDir, sort, prefix, filename);
                        }
                        for (int i = 1; i <= readStructure.molecularBarcode.length(); ++i) {
                            final String filename = filePrefix + ".index_" + i + ".fastq";
                            compareResults(testDataDir, sort, prefix, filename);
                        }

                        if (otfBarcodeExtract) {
                            final MetricsFile<BarcodeMetric, Integer> metrics = new MetricsFile<>();
                            metrics.read(new FileReader(metricsFile));
                            Assert.assertEquals(metrics.getMetrics().get(3).PERFECT_MATCHES, expectedPfMatches);
                            Assert.assertEquals(metrics.getMetrics().get(3).PCT_MATCHES, expectedPctMatches);
                        }
                    }
                }
            }
        } finally {
            IOUtil.recursiveDelete(outputDir.toPath());
        }
    }

    private void compareResults(File testDataDir, boolean sort, File prefix, String filename) {
        File actual = new File(prefix.getParentFile(), filename);
        File expected = new File(testDataDir, filename);
        if(actual.length() == 0 && expected.length() == 0) return;
        if (sort) {
            IOUtil.assertFilesEqual(actual, expected);
        } else {
            compareFastqs(actual, expected);
        }
    }

    private void convertParamsFile(String libraryParamsFile, int concatNColumnFields, File testDataDir, File outputDir, File libraryParams, List<File> outputPrefixes, boolean multilane) throws FileNotFoundException {
        try (LineReader reader = new BufferedLineReader(new FileInputStream(new File(testDataDir, libraryParamsFile)))) {
            final PrintWriter writer = new PrintWriter(libraryParams);
            final String header = reader.readLine();
            writer.println(header + "\tOUTPUT_PREFIX");
            while (true) {
                final String line = reader.readLine();
                if (line == null) {
                    break;
                }
                final String[] fields = line.split("\t");
                String prefix = StringUtil.join("", Arrays.copyOfRange(fields, 0, concatNColumnFields));
                if (multilane) {
                    prefix += ".multilane";
                }
                final File outputPrefix = new File(outputDir, prefix);
                outputPrefixes.add(outputPrefix);
                writer.println(line + "\t" + outputPrefix);
            }
            writer.close();
        }
    }
}
