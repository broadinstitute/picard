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
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryAlignmentFilter;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.AbstractLocusIterator;
import htsjdk.samtools.util.EdgeReadIterator;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import picard.filter.CountingDuplicateFilter;
import picard.filter.CountingFilter;
import picard.filter.CountingMapQFilter;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;


/**
 * Contains util methods for CollectWgsMetricsTest, CollectWgsMetricsWithNonZeroCoverageTest
 */

public class CollectWgsMetricsTestUtils {

    private static final String sqHeaderLN20 = "@HD\tSO:coordinate\tVN:1.0\n@SQ\tSN:chrM\tAS:HG18\tLN:20\n";
    private static final String s1 = "3851612\t16\tchrM\t1\t255\t3M2D10M\t*\t0\t0\tACCTACGTTCAAT\tDDDDDDDDDDDDD\n";
    private static final String sqHeaderLN100 = "@HD\tSO:coordinate\tVN:1.0\n@SQ\tSN:chrM\tAS:HG18\tLN:100\n";
    private static final String s2 = "3851612\t16\tchrM\t2\t255\t10M1D5M\t*\t0\t0\tCCTACGTTCAATATT\tDDDDDDDDDDDDDDD\n";
    static final String exampleSamOneRead = sqHeaderLN20+s1;
    static final String exampleSamTwoReads = sqHeaderLN100+s1+s2;

    protected static SAMRecordSetBuilder createTestSAMBuilder(final File reference,
                                                              final String readGroupId,
                                                              final String sample,
                                                              final String platform,
                                                              final String library) {
        final SAMFileHeader header = new SAMFileHeader();

        // check that dictionary file is readable and then set header dictionary
        try {
            header.setSequenceDictionary(SAMSequenceDictionaryExtractor.extractDictionary(reference.toPath()));
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

    static IntervalList createIntervalList() {
        return new IntervalList(new SAMFileHeader());
    }

    static AbstractLocusIterator createReadEndsIterator(String exampleSam){
        final List<SamRecordFilter> filters = new ArrayList<>();
        final CountingFilter dupeFilter = new CountingDuplicateFilter();
        final CountingFilter mapqFilter = new CountingMapQFilter(0);
        filters.add(new SecondaryAlignmentFilter()); // Not a counting filter because we never want to count reads twice
        filters.add(mapqFilter);
        filters.add(dupeFilter);
        SamReader samReader = createSamReader(exampleSam);
        AbstractLocusIterator iterator =  new EdgeReadIterator(samReader);
        iterator.setSamFilters(filters);
        iterator.setMappingQualityScoreCutoff(0); // Handled separately because we want to count bases
        iterator.setIncludeNonPfReads(false);
        return iterator;
    }

    static SamReader createSamReader(String samExample) {
        ByteArrayInputStream inputStream = new ByteArrayInputStream(samExample.getBytes());
        return SamReaderFactory.makeDefault().open(SamInputResource.of(inputStream));
    }

    static ReferenceSequenceFileWalker getReferenceSequenceFileWalker(String referenceString){
        ReferenceSequenceFile referenceSequenceFile = createReferenceSequenceFile(referenceString);
        return new ReferenceSequenceFileWalker(referenceSequenceFile);
    }

    static ReferenceSequenceFileWalker getReferenceSequenceFileWalker(){
        String referenceString = ">ref\nACCTACGTTCAATATTCTTC";
        return getReferenceSequenceFileWalker(referenceString);
    }

    static ReferenceSequenceFile createReferenceSequenceFile(String referenceString) {
        final SAMSequenceRecord record = new SAMSequenceRecord("ref", referenceString.length());
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(record);
        return new ReferenceSequenceFile() {

            boolean done = false;
            @Override
            public SAMSequenceDictionary getSequenceDictionary() {
                return dictionary;
            }

            @Override
            public ReferenceSequence nextSequence() {
                if (!done) {
                    done = true;
                    return getSequence(record.getSequenceName());
                }
                return null;
            }

            @Override
            public void reset() {
                done=false;
            }

            @Override
            public boolean isIndexed() {
                return false;
            }

            @Override
            public ReferenceSequence getSequence(String contig) {
                if (contig.equals(record.getSequenceName())) {
                    return new ReferenceSequence(record.getSequenceName(), 0, referenceString.getBytes());
                } else {
                    return null;
                }
            }

            @Override
            public ReferenceSequence getSubsequenceAt(String contig, long start, long stop) {
                return null;
            }

            @Override
            public String toString() {
                return null;
            }

            @Override
            public void close() throws IOException {

            }
        };
    }

    // call main (from IDE for example) to createTestSAM(), which creates a test SAM file
    public static void main(String[] args) {
        try { createTestSAM("TestSam"); } catch(IOException e) { ; }
    }
}
