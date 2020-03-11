/*
 * The MIT License
 *
 * Copyright (c) 2016 Nils Homer
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
 *
 */

package picard.vcf;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Iterator;

/**
 * Tests for ByIntervalListVariantContextIterator.
 */
public class ByIntervalListVariantContextIteratorTest {

    private static final File TEST_DATA_PATH       = new File("testdata/picard/vcf/");
    private static final File CEU_TRIOS_SNPS_VCF   = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");
    private static final File CEU_TRIOS_INDELS_VCF = new File(TEST_DATA_PATH, "CEUTrio-indels.vcf");
    private static final File EMPTY_VCF            = new File(TEST_DATA_PATH, "empty.vcf");

    private final SAMFileHeader header;
    private final SAMSequenceDictionary dict;

    public ByIntervalListVariantContextIteratorTest() {
        this.header = getSAMFileHeader();
        this.dict   = header.getSequenceDictionary();
    }

    private VCFFileReader getReader(final File file) {
        return new VCFFileReader(file, true);
    }

    private SAMFileHeader getSAMFileHeader() {
        final VCFFileReader reader = getReader(CEU_TRIOS_SNPS_VCF);
        final SAMSequenceDictionary dict = reader.getFileHeader().getSequenceDictionary();
        reader.close();
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(dict);
        return header;
    }

    @Test
    public void testEmptyIntervalList() {
        final IntervalList intervalList         = new IntervalList(header);
        final VCFFileReader reader              = getReader(CEU_TRIOS_SNPS_VCF);
        final Iterator<VariantContext> iterator = new ByIntervalListVariantContextIterator(reader, intervalList);
        Assert.assertFalse(iterator.hasNext());
        reader.close();
    }

    @Test
    public void testNoVariants() {
        final IntervalList intervalList         = new IntervalList(header);
        intervalList.add(new Interval(this.dict.getSequence(0).getSequenceName(), 1, 100));
        final VCFFileReader reader              = getReader(EMPTY_VCF);
        final Iterator<VariantContext> iterator = new ByIntervalListVariantContextIterator(reader, intervalList);
        Assert.assertFalse(iterator.hasNext());
        reader.close();
    }

    @Test
    public void testSimpleOverlap() {
        final IntervalList intervalList         = new IntervalList(header);
        intervalList.add(new Interval("2", 167166899, 167166899));
        final VCFFileReader reader              = getReader(CEU_TRIOS_SNPS_VCF);
        final Iterator<VariantContext> iterator = new ByIntervalListVariantContextIterator(reader, intervalList);
        Assert.assertTrue(iterator.hasNext());
        final VariantContext ctx = iterator.next();
        Assert.assertEquals(ctx.getStart(), 167166899);
        Assert.assertFalse(iterator.hasNext());
        reader.close();
    }

    @Test
    public void testNoOverlapDifferentContig() {
        final IntervalList intervalList         = new IntervalList(header);
        intervalList.add(new Interval("3", 167166899, 167166899));
        final VCFFileReader reader              = getReader(CEU_TRIOS_SNPS_VCF);
        final Iterator<VariantContext> iterator = new ByIntervalListVariantContextIterator(reader, intervalList);
        Assert.assertFalse(iterator.hasNext());
        reader.close();
    }

    @Test
    public void testSimpleEnclosing() {
        final IntervalList intervalList         = new IntervalList(header);
        intervalList.add(new Interval("12", 68921962, 68921962)); // deletion spans this
        final VCFFileReader reader              = getReader(CEU_TRIOS_INDELS_VCF);
        final Iterator<VariantContext> iterator = new ByIntervalListVariantContextIterator(reader, intervalList);
        Assert.assertTrue(iterator.hasNext());
        final VariantContext ctx = iterator.next();
        Assert.assertEquals(ctx.getStart(), 68921960);
        Assert.assertEquals(ctx.getEnd(), 68921966);
        Assert.assertFalse(iterator.hasNext());
        reader.close();
    }

    @Test
    public void testVariantOverlappingMultipleIntervalsIsReturnedOnlyOnce() {
        final IntervalList intervalList         = new IntervalList(header);
        intervalList.add(new Interval("12", 68921962, 68921962)); // deletion spans this
        intervalList.add(new Interval("12", 68921964, 68921964)); // deletion spans this
        final VCFFileReader reader              = getReader(CEU_TRIOS_INDELS_VCF);
        final Iterator<VariantContext> iterator = new ByIntervalListVariantContextIterator(reader, intervalList);
        Assert.assertTrue(iterator.hasNext());
        final VariantContext ctx = iterator.next();
        Assert.assertEquals(ctx.getStart(), 68921960);
        Assert.assertEquals(ctx.getEnd(), 68921966);
        Assert.assertFalse(iterator.hasNext());
        reader.close();
    }
}
