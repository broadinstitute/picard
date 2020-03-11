package picard.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.util.PeekableIterator;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Unit tests for QuerySortedReadPairIteratorUtil
 */
public class QuerySortedReadPairIteratorUtilTest {
    private static final int READ_LENGTH = 20;

    @Test
    public void testBasicPairedRead() {
        SAMRecordSetBuilder builder = new SAMRecordSetBuilder(false, SAMFileHeader.SortOrder.queryname);
        builder.setReadLength(READ_LENGTH);
        builder.addPair("mapped_paired", 1, 1, 31);
        PeekableIterator<SAMRecord> iterator = new PeekableIterator<SAMRecord>(builder.iterator());

        QuerySortedReadPairIteratorUtil.ReadPair pair = QuerySortedReadPairIteratorUtil.getNextReadPair(iterator);
        Assert.assertNotNull(pair);
        Assert.assertNotNull(pair.read1);
        Assert.assertNotNull(pair.read2);
        Assert.assertEquals("mapped_paired", pair.read1.getReadName());
        Assert.assertEquals("mapped_paired", pair.read2.getReadName());

        pair = QuerySortedReadPairIteratorUtil.getNextReadPair(iterator);
        Assert.assertNull(pair);
    }

    @Test
    public void testBasicUnmappedReadPair() {
        SAMRecordSetBuilder builder = new SAMRecordSetBuilder(false, SAMFileHeader.SortOrder.queryname);
        builder.setReadLength(READ_LENGTH);
        builder.addUnmappedPair("unmapped_paired");
        PeekableIterator<SAMRecord> iterator = new PeekableIterator<SAMRecord>(builder.iterator());

        QuerySortedReadPairIteratorUtil.ReadPair pair = QuerySortedReadPairIteratorUtil.getNextReadPair(iterator);
        Assert.assertNotNull(pair);
        Assert.assertNotNull(pair.read1);
        Assert.assertNotNull(pair.read2);
        Assert.assertEquals("unmapped_paired", pair.read1.getReadName());
        Assert.assertEquals("unmapped_paired", pair.read2.getReadName());

        pair = QuerySortedReadPairIteratorUtil.getNextReadPair(iterator);
        Assert.assertNull(pair);
    }

    @Test
    public void testBasicHalfmappedReadPair() {
        SAMRecordSetBuilder builder = new SAMRecordSetBuilder(false, SAMFileHeader.SortOrder.queryname);
        builder.setReadLength(READ_LENGTH);
        builder.addPair("halfmapped_paired", 1, 1, 31, false, true, "20M", "20M", true, false, 20);
        PeekableIterator<SAMRecord> iterator = new PeekableIterator<SAMRecord>(builder.iterator());

        QuerySortedReadPairIteratorUtil.ReadPair pair = QuerySortedReadPairIteratorUtil.getNextReadPair(iterator);
        Assert.assertNotNull(pair);
        Assert.assertNotNull(pair.read1);
        Assert.assertNotNull(pair.read2);
        Assert.assertEquals("halfmapped_paired", pair.read1.getReadName());
        Assert.assertEquals("halfmapped_paired", pair.read2.getReadName());

        pair = QuerySortedReadPairIteratorUtil.getNextReadPair(iterator);
        Assert.assertNull(pair);
    }

    @Test
    public void testFragmentNoReadPair() {
        SAMRecordSetBuilder builder = new SAMRecordSetBuilder(false, SAMFileHeader.SortOrder.queryname);
        builder.setReadLength(READ_LENGTH);
        builder.addFrag("mapped_frag_a", 1, 1, false);
        builder.addFrag("mapped_frag_b", 1, 1, false);
        PeekableIterator<SAMRecord> iterator = new PeekableIterator<SAMRecord>(builder.iterator());

        QuerySortedReadPairIteratorUtil.ReadPair pair = QuerySortedReadPairIteratorUtil.getNextReadPair(iterator);
        Assert.assertNotNull(pair);
        Assert.assertNotNull(pair.read1);
        Assert.assertNull(pair.read2);
        Assert.assertEquals("mapped_frag_a", pair.read1.getReadName());

        pair = QuerySortedReadPairIteratorUtil.getNextReadPair(iterator);
        Assert.assertNotNull(pair);
        Assert.assertNotNull(pair.read1);
        Assert.assertNull(pair.read2);
        Assert.assertEquals("mapped_frag_b", pair.read1.getReadName());

        pair = QuerySortedReadPairIteratorUtil.getNextReadPair(iterator);
        Assert.assertNull(pair);
    }
}
