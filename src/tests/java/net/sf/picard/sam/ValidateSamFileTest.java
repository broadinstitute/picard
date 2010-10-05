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

package net.sf.picard.sam;

import net.sf.picard.PicardException;
import net.sf.picard.metrics.MetricBase;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.util.Histogram;
import net.sf.samtools.*;
import net.sf.samtools.util.StringUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.util.Arrays;
import java.util.Iterator;

/**
 * Tests almost all error conditions detected by the sam file validator. The
 * conditions not tested are proactively prevented by sam generation code.
 *
 * @author Doug Voet
 */
public class ValidateSamFileTest {
    private static final File TEST_DATA_DIR = new File("testdata/net/sf/picard/sam/ValidateSamFileTest");

    @Test
    public void testSortOrder() throws IOException {
        Histogram<String> results = executeValidation(new SAMFileReader(new File(TEST_DATA_DIR, "invalid_coord_sort_order.sam")), null);
        Assert.assertEquals(results.get(SAMValidationError.Type.RECORD_OUT_OF_ORDER.getHistogramString()).getValue(), 1.0);
        results = executeValidation(new SAMFileReader(new File(TEST_DATA_DIR, "invalid_queryname_sort_order.sam")), null);
        Assert.assertEquals(results.get(SAMValidationError.Type.RECORD_OUT_OF_ORDER.getHistogramString()).getValue(), 5.0);
    }
    
    @Test
    public void testVerbose() throws IOException {
        final SAMRecordSetBuilder samBuilder = new SAMRecordSetBuilder();
        
        for (int i=0; i<20; i++) {
            samBuilder.addFrag(String.valueOf(i), 1, i, false);
        }
        for (final SAMRecord record : samBuilder) {
            record.setProperPairFlag(true);
        }
        
        final StringWriter results = new StringWriter();
        final SamFileValidator validator = new SamFileValidator(new PrintWriter(results));
        validator.setVerbose(true, 10);
        validator.validateSamFileVerbose(
                samBuilder.getSamReader(), 
                null);
        
        final int lineCount = results.toString().split("\n").length;
        Assert.assertEquals(lineCount, 11);
    }
    
    @Test
    public void testUnpairedRecords() throws IOException {
        final SAMRecordSetBuilder samBuilder = new SAMRecordSetBuilder();
        
        for (int i=0; i<6; i++) {
            samBuilder.addFrag(String.valueOf(i), i, i, false);
        }
        final Iterator<SAMRecord> records = samBuilder.iterator();
        records.next().setProperPairFlag(true);
        records.next().setMateUnmappedFlag(true);
        records.next().setMateNegativeStrandFlag(true);
        records.next().setFirstOfPairFlag(true);
        records.next().setSecondOfPairFlag(true);
        records.next().setMateReferenceIndex(1);
        
        final Histogram<String> results = executeValidation(samBuilder.getSamReader(), null);
        
        Assert.assertEquals(results.get(SAMValidationError.Type.INVALID_FLAG_PROPER_PAIR.getHistogramString()).getValue(), 1.0);
        Assert.assertEquals(results.get(SAMValidationError.Type.INVALID_FLAG_MATE_UNMAPPED.getHistogramString()).getValue(), 1.0);
        Assert.assertEquals(results.get(SAMValidationError.Type.INVALID_FLAG_MATE_NEG_STRAND.getHistogramString()).getValue(), 1.0);
        Assert.assertEquals(results.get(SAMValidationError.Type.INVALID_FLAG_FIRST_OF_PAIR.getHistogramString()).getValue(), 1.0);
        Assert.assertEquals(results.get(SAMValidationError.Type.INVALID_FLAG_SECOND_OF_PAIR.getHistogramString()).getValue(), 1.0);
        Assert.assertEquals(results.get(SAMValidationError.Type.INVALID_MATE_REF_INDEX.getHistogramString()).getValue(), 1.0);
    }

    @Test
    public void testPairedRecords() throws IOException {
        final SAMRecordSetBuilder samBuilder = new SAMRecordSetBuilder();
        
        for (int i=0; i<5; i++) {
            samBuilder.addPair(String.valueOf(i), i, i, i+100);
        }
        final Iterator<SAMRecord> records = samBuilder.iterator();
        records.next().setMateReferenceName("*");
        records.next().setMateAlignmentStart(Integer.MAX_VALUE);
        records.next().setMateAlignmentStart(records.next().getAlignmentStart()+1);
        records.next().setMateNegativeStrandFlag(!records.next().getReadNegativeStrandFlag());
        records.next().setMateReferenceIndex(records.next().getReferenceIndex() + 1);
        records.next().setMateUnmappedFlag(!records.next().getReadUnmappedFlag());

        
        final Histogram<String> results = executeValidation(samBuilder.getSamReader(), null);
        
        Assert.assertEquals(results.get(SAMValidationError.Type.INVALID_ALIGNMENT_START.getHistogramString()).getValue(), 3.0);
        Assert.assertEquals(results.get(SAMValidationError.Type.INVALID_FLAG_MATE_UNMAPPED.getHistogramString()).getValue(), 1.0);
        Assert.assertEquals(results.get(SAMValidationError.Type.MISMATCH_FLAG_MATE_NEG_STRAND.getHistogramString()).getValue(), 1.0);
        Assert.assertEquals(results.get(SAMValidationError.Type.MISMATCH_FLAG_MATE_UNMAPPED.getHistogramString()).getValue(), 1.0);
        Assert.assertEquals(results.get(SAMValidationError.Type.MISMATCH_MATE_ALIGNMENT_START.getHistogramString()).getValue(), 2.0);
        Assert.assertEquals(results.get(SAMValidationError.Type.MISMATCH_MATE_REF_INDEX.getHistogramString()).getValue(), 2.0);
    }

    @Test
    public void testUnmappedRecords() throws IOException {
        final SAMRecordSetBuilder samBuilder = new SAMRecordSetBuilder();
        
        for (int i=0; i<4; i++) {
            samBuilder.addUnmappedFragment(String.valueOf(i));
        }
        final Iterator<SAMRecord> records = samBuilder.iterator();
        records.next().setReadNegativeStrandFlag(true);
        records.next().setNotPrimaryAlignmentFlag(true);
        records.next().setMappingQuality(10);
        records.next().setCigarString("36M");
        
        final Histogram<String> results = executeValidation(samBuilder.getSamReader(), null);
        
        Assert.assertEquals(results.get(SAMValidationError.Type.INVALID_FLAG_NOT_PRIM_ALIGNMENT.getHistogramString()).getValue(), 1.0);
        Assert.assertEquals(results.get(SAMValidationError.Type.INVALID_MAPPING_QUALITY.getHistogramString()).getValue(), 1.0);
        Assert.assertEquals(results.get(SAMValidationError.Type.INVALID_CIGAR.getHistogramString()).getValue(), 1.0);
    }

    @Test
    public void testMappedRecords() throws IOException {
        final SAMRecordSetBuilder samBuilder = new SAMRecordSetBuilder();
        
        for (int i=0; i<2; i++) {
            samBuilder.addFrag(String.valueOf(i), i, i, false);
        }
        final Iterator<SAMRecord> records = samBuilder.iterator();
        records.next().setCigarString("25M3S25M");
        records.next().setReferenceName("*");
        
        final Histogram<String> results = executeValidation(samBuilder.getSamReader(), null);
        
        Assert.assertEquals(results.get(SAMValidationError.Type.INVALID_CIGAR.getHistogramString()).getValue(), 1.0);
        Assert.assertEquals(results.get(SAMValidationError.Type.INVALID_FLAG_READ_UNMAPPED.getHistogramString()).getValue(), 1.0);
        Assert.assertEquals(results.get(SAMValidationError.Type.MISSING_TAG_NM.getHistogramString()).getValue(), 1.0);
    }
    
    @Test
    public void testNmFlagValidation() throws IOException {
        final SAMRecordSetBuilder samBuilder = new SAMRecordSetBuilder();
        
        for (int i=0; i<3; i++) {
            samBuilder.addFrag(String.valueOf(i), i, i+1, false);
        }
        final Iterator<SAMRecord> records = samBuilder.iterator();
        records.next().setAttribute(ReservedTagConstants.NM, 4);

        // PIC-215: Confirm correct NM value when there is an insertion and a deletion.
        final SAMRecord recordWithInsert = records.next();
        final byte[] sequence = recordWithInsert.getReadBases();
        Arrays.fill(sequence, (byte)'A');
        recordWithInsert.setReadBases(sequence);
        recordWithInsert.setCigarString("1D" + Integer.toString(sequence.length-1) + "M1I");
        recordWithInsert.setAttribute(ReservedTagConstants.NM, 2);
        
        final Histogram<String> results = executeValidation(samBuilder.getSamReader(), new ReferenceSequenceFile() {
            private int index=0;
            public SAMSequenceDictionary getSequenceDictionary() {
                return null;
            }

            public ReferenceSequence nextSequence() {
                final byte[] bases = new byte[10000];
                Arrays.fill(bases, (byte) 'A'); 
                return new ReferenceSequence("foo", index++, bases);
            }

            public void reset() {
                this.index = 0;
            }

            public boolean isIndexed() { return false; }

            public ReferenceSequence getSequence(final String contig) {
                throw new UnsupportedOperationException();
            }

            public ReferenceSequence getSubsequenceAt(final String contig, final long start, final long stop) {
                throw new UnsupportedOperationException();
            }
        });
        
        Assert.assertEquals(results.get(SAMValidationError.Type.INVALID_TAG_NM.getHistogramString()).getValue(), 1.0);
        Assert.assertEquals(results.get(SAMValidationError.Type.MISSING_TAG_NM.getHistogramString()).getValue(), 1.0);
    }

    @Test(dataProvider = "testTruncatedScenarios")
    public void testTruncated(final String scenario, final String inputFile, final SAMValidationError.Type expectedError)
            throws Exception {
        final SAMFileReader reader = new SAMFileReader(new File(TEST_DATA_DIR, inputFile));
        final Histogram<String> results = executeValidation(reader, null);
        Assert.assertNotNull(results.get(expectedError.getHistogramString()));
        Assert.assertEquals(results.get(expectedError.getHistogramString()).getValue(), 1.0);
    }

    @DataProvider(name = "testTruncatedScenarios")
    public Object[][] testTruncatedScenarios() {
        return new Object[][] {
                {"truncated bam", "truncated.bam", SAMValidationError.Type.TRUNCATED_FILE},
                {"truncated quals", "truncated_quals.sam", SAMValidationError.Type.MISMATCH_READ_LENGTH_AND_QUALS_LENGTH},
                // TODO: Because validation is turned off when parsing, this error is not detectable currently by validator.
                //{"truncated tag", "truncated_tag.sam", SAMValidationError.Type.TRUNCATED_FILE},
                // TODO: Currently, this is not considered an error.  Should it be?
                //{"hanging tab", "hanging_tab.sam", SAMValidationError.Type.TRUNCATED_FILE},
        };
    }

    @Test(expectedExceptions = PicardException.class, dataProvider = "testFatalParsingErrors")
    public void testFatalParsingErrors(final String scenario, final String inputFile) throws Exception {
        final SAMFileReader reader = new SAMFileReader(new File(TEST_DATA_DIR, inputFile));
        executeValidation(reader, null);
        Assert.fail("Exception should have been thrown.");
    }

    @DataProvider(name = "testFatalParsingErrors")
    public Object[][] testFatalParsingErrorScenarios() {
        return new Object[][] {
                {"missing fields", "missing_fields.sam"},
                {"zero length read", "zero_length_read.sam"}
        };
    }

    @Test
    public void testHeaderVersionValidation() throws Exception {
        final String header = "@HD	VN:Hi,Mom!	SO:queryname";
        final InputStream strm = new ByteArrayInputStream(StringUtil.stringToBytes(header));
        final SAMFileReader samReader = new SAMFileReader(strm);
        final Histogram<String> results = executeValidation(samReader, null);
        Assert.assertEquals(results.get(SAMValidationError.Type.INVALID_VERSION_NUMBER.getHistogramString()).getValue(), 1.0);
    }

    @Test
    public void testCigarOffEndOfReferenceValidation() throws Exception {
        final SAMRecordSetBuilder samBuilder = new SAMRecordSetBuilder();
        samBuilder.addFrag(String.valueOf(0), 0, 1, false);
        final int contigLength = samBuilder.getHeader().getSequence(0).getSequenceLength();
        // Should hang off the end.
        samBuilder.addFrag(String.valueOf(1), 0, contigLength - 1, false);
        final Histogram<String> results = executeValidation(samBuilder.getSamReader(), null);
        Assert.assertNotNull(results.get(SAMValidationError.Type.CIGAR_MAPS_OFF_REFERENCE.getHistogramString()));
        Assert.assertEquals(results.get(SAMValidationError.Type.CIGAR_MAPS_OFF_REFERENCE.getHistogramString()).getValue(), 1.0);
    }

    @Test(expectedExceptions = SAMFormatException.class)
    public void testConflictingTags() throws Exception {
        final String header = "@HD	VN:1.0	SO:queryname	SO:coordinate";
        final InputStream strm = new ByteArrayInputStream(StringUtil.stringToBytes(header));
        final SAMFileReader samReader = new SAMFileReader(strm);
        Assert.fail("Exception should have been thrown.");
    }

    @Test
    public void testRedundantTags() throws Exception {
        final String header = "@HD	VN:1.0	SO:coordinate	SO:coordinate";
        final InputStream strm = new ByteArrayInputStream(StringUtil.stringToBytes(header));
        final SAMFileReader samReader = new SAMFileReader(strm);
        Assert.assertEquals(SAMFileHeader.SortOrder.coordinate, samReader.getFileHeader().getSortOrder());
    }

    @Test
    public void testHeaderValidation() throws Exception {
        SAMFileReader.ValidationStringency saveStringency = SAMFileReader.getDefaultValidationStringency();
        SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        try {
            final SAMFileReader samReader = new SAMFileReader(new File(TEST_DATA_DIR, "buggyHeader.sam"));
            Histogram<String> results = executeValidation(samReader, null);
            Assert.assertEquals(results.get(SAMValidationError.Type.UNRECOGNIZED_HEADER_TYPE.getHistogramString()).getValue(), 3.0);
            Assert.assertEquals(results.get(SAMValidationError.Type.HEADER_TAG_MULTIPLY_DEFINED.getHistogramString()).getValue(), 1.0);
        } finally {
            SAMFileReader.setDefaultValidationStringency(saveStringency);
        }
    }

    @Test
    public void testIndexFileValidation() throws Exception {
        SAMFileReader.ValidationStringency saveStringency = SAMFileReader.getDefaultValidationStringency();
        SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        try {
            final SAMFileReader samReader = new SAMFileReader(new File(TEST_DATA_DIR, "bad_index.bam"));
            samReader.enableIndexCaching(true);
            Histogram<String> results = executeValidation(samReader, null);
            Assert.assertEquals(results.get(SAMValidationError.Type.INVALID_INDEX_FILE_POINTER.getHistogramString()).getValue(), 1.0);
        } finally {
            SAMFileReader.setDefaultValidationStringency(saveStringency);
        }
    }
 
    private Histogram<String> executeValidation(final SAMFileReader samReader, final ReferenceSequenceFile reference) throws IOException {
        final File outFile = File.createTempFile("validation", ".txt");
        final PrintWriter out = new PrintWriter(outFile);
        new SamFileValidator(out).setValidateIndex(true).validateSamFileSummary(samReader, reference);
        LineNumberReader reader = new LineNumberReader(new FileReader(outFile));
        if (reader.readLine().equals("No errors found")) {
            return new Histogram<String>();
        }
        final MetricsFile<MetricBase, String> outputFile = new MetricsFile<MetricBase, String>();
        outputFile.read(new FileReader(outFile));
        Assert.assertNotNull(outputFile.getHistogram());
        return outputFile.getHistogram();
    }

    @Test(dataProvider = "headerVersions")
    public void testHeaderVersion(final String version, boolean expectValid) throws Exception {
        final File samFile = File.createTempFile("validateHeader.", ".sam");
        samFile.deleteOnExit();
        final PrintWriter pw = new PrintWriter(samFile);
        pw.println("@HD\tVN:" + version);
        pw.close();
        final SAMFileReader reader = new SAMFileReader(samFile);
        final Histogram<String> results = executeValidation(reader, null);
        if (expectValid) Assert.assertNull(results.get(SAMValidationError.Type.INVALID_VERSION_NUMBER.getHistogramString()));
        else {
            Assert.assertNotNull(results.get(SAMValidationError.Type.INVALID_VERSION_NUMBER.getHistogramString()));
            Assert.assertEquals(results.get(SAMValidationError.Type.INVALID_VERSION_NUMBER.getHistogramString()).getValue(), 1.0);
        }
    }

    @DataProvider(name = "headerVersions")
    public Object[][] testHeaderVersionScenarios() {
        return new Object[][] {
                {"1.0", true},
                {"1.3", true},
                {"1.4", false},
        };
    }

}
