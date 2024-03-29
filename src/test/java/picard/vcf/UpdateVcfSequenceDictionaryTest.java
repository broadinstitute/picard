/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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
package picard.vcf;

import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import java.io.*;
import java.lang.reflect.Field;


public class UpdateVcfSequenceDictionaryTest {
    private static final File TEST_DATA_PATH = new File("testdata/picard/vcf/");
    private static final File INPUT_FILE = new File(TEST_DATA_PATH, "vcfFormatTest.vcf");

    // vcfFormatTest.bad_dict.vcf is a vcf with two (2) ##contig lines deleted
    private final File SAM_SEQUENCE_DICTIONARY_VCF = new File(TEST_DATA_PATH, "vcfFormatTest.bad_dict.vcf");

    private static final File OUTPUT_DATA_PATH = IOUtil.createTempDir("UpdateVcfSequenceDictionaryTest.tmp").toFile();
    private static final File STD_OUT_FILE = new File(OUTPUT_DATA_PATH, "stdout.vcf");
    private static final String STD_OUT_NAME = "/dev/stdout";

    @AfterClass
    public void teardown() {
        IOUtil.deleteDirectoryTree(OUTPUT_DATA_PATH);
    }

    @DataProvider(name = "OutputFiles")
    public static Object[][] outputFies() {

        return new Object[][] {
                {OUTPUT_DATA_PATH + "updateVcfSequenceDictionaryTest-delete-me"  + IOUtil.COMPRESSED_VCF_FILE_EXTENSION},
                {OUTPUT_DATA_PATH + "updateVcfSequenceDictionaryTest-delete-me" + IOUtil.VCF_FILE_EXTENSION},
                {STD_OUT_NAME}
        };
    }

    /**
     * Utility for unzipping a zipped file's contents into a human readable, unzipped file
     *
     * @param zippedFile    input zipped file
     * @param unzippedFile  unzipped file
     * @throws IOException
     */
    private void unzipFile(final File zippedFile, final File unzippedFile) throws IOException {
        final InputStream gzInputStream = IOUtil.openFileForReading(zippedFile);
        final FileOutputStream fileOutputStream = new FileOutputStream(unzippedFile.getAbsolutePath());
        final byte[] buffer = new byte[1024];
        int len;
        while ((len = gzInputStream.read(buffer)) > 0) {
            fileOutputStream.write(buffer, 0, len);
        }
        gzInputStream.close();
        fileOutputStream.close();
    }

    @Test(dataProvider = "OutputFiles")
    public void testUpdateVcfSequenceDictionary(final String outputFileName) throws IOException, NoSuchFieldException, IllegalAccessException {
        File outputFile = new File(outputFileName);
        outputFile.deleteOnExit();

        final UpdateVcfSequenceDictionary updateVcfSequenceDictionary = new UpdateVcfSequenceDictionary();

        updateVcfSequenceDictionary.INPUT = INPUT_FILE;
        updateVcfSequenceDictionary.SEQUENCE_DICTIONARY = SAM_SEQUENCE_DICTIONARY_VCF;
        if (outputFileName.equals(STD_OUT_NAME)) {
            final FileOutputStream stream = new FileOutputStream(STD_OUT_FILE);

            // Ugliness required to write to a stream given as a string on the commandline.
            // Since the actual fd number is private inside FileDescriptor, needs reflection
            // in order to pull it out.

            final Field fdField = FileDescriptor.class.getDeclaredField("fd");
            fdField.setAccessible(true);
            updateVcfSequenceDictionary.OUTPUT = new File("/dev/fd/" + fdField.getInt(stream.getFD()));
        } else {
            updateVcfSequenceDictionary.OUTPUT = outputFile;
        }

        Assert.assertEquals(updateVcfSequenceDictionary.instanceMain(new String[0]), 0);

        if (outputFileName.equals(STD_OUT_NAME)) {
            outputFile = STD_OUT_FILE;
        }
        else if (outputFileName.endsWith(IOUtil.COMPRESSED_VCF_FILE_EXTENSION)) {
            // Output is a gzipped vcf, unzip it
            File unzippedOutputFile = new File(OUTPUT_DATA_PATH + ".out.vcf");
            unzippedOutputFile.deleteOnExit();
            unzipFile(outputFile, unzippedOutputFile);
            outputFile = unzippedOutputFile;
        }

        IOUtil.assertFilesEqual(SAM_SEQUENCE_DICTIONARY_VCF, outputFile);

        // A little extra checking.
        Assert.assertEquals(SAMSequenceDictionaryExtractor.extractDictionary(INPUT_FILE.toPath()).size(), 84);
        Assert.assertEquals(SAMSequenceDictionaryExtractor.extractDictionary(SAM_SEQUENCE_DICTIONARY_VCF.toPath()).size(), 82);
        Assert.assertEquals(SAMSequenceDictionaryExtractor.extractDictionary(outputFile.toPath()).size(), 82);
    }


    private String classPath = "\"" + System.getProperty("java.class.path") + "\" ";

    @Test
    public void testCaptureStdout() throws IOException {
        final File outputFile = File.createTempFile("UpdateVcfSequenceDictionaryTest.output.", ".vcf");
        outputFile.deleteOnExit();

        String[] command = {
                "/bin/bash",
                "-c",
                "java -Dpicard.useLegacyParser=false -classpath " +
                        classPath +
                        "picard.cmdline.PicardCommandLine " +
                        "UpdateVcfSequenceDictionary " +
                        "-I " + INPUT_FILE.getAbsolutePath() + " " +
                        "-O /dev/stdout " +
                        "-SD " + SAM_SEQUENCE_DICTIONARY_VCF.getAbsolutePath()
        };

        try {
            ProcessBuilder processBuilder = new ProcessBuilder(command);
            processBuilder.inheritIO();
            processBuilder.redirectOutput(outputFile);

            Process process = processBuilder.start();
            Assert.assertEquals(process.waitFor(), 0);

            IOUtil.assertFilesEqual(SAM_SEQUENCE_DICTIONARY_VCF, outputFile);

            // A little extra checking.
            Assert.assertEquals(SAMSequenceDictionaryExtractor.extractDictionary(INPUT_FILE.toPath()).size(), 84);
            Assert.assertEquals(SAMSequenceDictionaryExtractor.extractDictionary(SAM_SEQUENCE_DICTIONARY_VCF.toPath()).size(), 82);
            Assert.assertEquals(SAMSequenceDictionaryExtractor.extractDictionary(outputFile.toPath()).size(), 82);
        } catch (Exception e) {
            Assert.fail("Failure executing UpdateVcfSequenceDictionary", e);
        }
    }
}