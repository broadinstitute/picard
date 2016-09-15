/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

package picard.analysis;


import htsjdk.samtools.*;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

import java.io.File;
import java.io.IOError;
import java.io.IOException;
import java.util.stream.IntStream;


/**
 * Contains util methods for CollectWgsMetricsTest, CollectWgsMetricsWithNonZeroCoverageTest
 */

public class CollectWgsMetricsTestUtils {
    protected static SAMRecordSetBuilder createTestSAMBuilder(final File reference,
                                                              final String readGroupId,
                                                              final String sample,
                                                              final String platform,
                                                              final String library) {
        final SAMFileHeader header = new SAMFileHeader();

        // check that dictionary file is readable and then set header dictionary
        try {
            header.setSequenceDictionary(SAMSequenceDictionaryExtractor.extractDictionary(reference));
            header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        } catch (final SAMException e) {
            e.printStackTrace();
        }

        // set readGroupRecord
        final SAMReadGroupRecord readGroupRecord = new SAMReadGroupRecord(readGroupId);
        readGroupRecord.setSample(sample);
        readGroupRecord.setPlatform(platform);
        readGroupRecord.setLibrary(library);
        readGroupRecord.setPlatformUnit(readGroupId);
        header.addReadGroup(readGroupRecord);

        final SAMRecordSetBuilder setBuilder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate, true, 100);
        setBuilder.setReadGroup(readGroupRecord);
        setBuilder.setUseNmFlag(true);
        setBuilder.setHeader(header);

        return(setBuilder);
    }


    /**
     * Template code for creating a custom SAM file for testing. Modify to suit your needs.
     */
    private static void createTestSAM(String testSamName) throws IOException {
        final File testDir = new File("testdata/picard/analysis/directed/CollectHsMetrics/");
        final File reference = new File("testdata/picard/quality/chrM.reference.fasta");
        final String readGroupId = "TestReadGroup";
        final String sample = "TestSample";
        final String platform = "Illumina";
        final String library = "TestLibrary";
        final int numReads = 1;
        final int readLength = 10;

        File samFile = File.createTempFile(testSamName, ".bam", testDir);

        final SAMRecordSetBuilder setBuilder = createTestSAMBuilder(reference, readGroupId, sample, platform, library);
        setBuilder.setReadLength(readLength);

        // add reads to the SAMRecordSetBuilder
        IntStream.range(0, numReads).forEach(i -> setBuilder.addPair("MediocreBaseQ" + i, 0, 1, 200, false, false, readLength + "M", readLength + "M", false, true, 40));

        final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(setBuilder.getHeader(), false, samFile);
        setBuilder.forEach(writer::addAlignment);
        writer.close();
    }

    // call main (from IDE for example) to createTestSAM(), which creates a test SAM file
    public static void main(String[] args) {
        try { createTestSAM("TestSam"); } catch(IOException e) { ; }
    }
}
