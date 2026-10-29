/*
 * The MIT License
 *
 * Copyright (c) 2026 The Broad Institute
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

package picard.sam.markduplicates.util;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SamRecordTrackingBuffer;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Tests for {@link SamRecordWithOrdinalAndSetDuplicateReadFlag}: the record and ordinal it holds, that its result
 * state is the record's duplicate flag, and that a {@link SamRecordTrackingBuffer} of these emits records in
 * ordinal order with the duplicate flags it was given.
 */
public class SamRecordWithOrdinalAndSetDuplicateReadFlagTest {

    @Test
    public void testConstructors() {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final SAMRecord rec = records.fragment("r", "rgA", 0, 10, false, "10M", 30);
        final SamRecordWithOrdinalAndSetDuplicateReadFlag wrapped = new SamRecordWithOrdinalAndSetDuplicateReadFlag(rec, Long.MAX_VALUE);
        Assert.assertSame(wrapped.getRecord(), rec);
        Assert.assertEquals(wrapped.getRecordOrdinal(), Long.MAX_VALUE);

        final SamRecordWithOrdinalAndSetDuplicateReadFlag empty = new SamRecordWithOrdinalAndSetDuplicateReadFlag();
        Assert.assertNull(empty.getRecord());
        empty.setRecord(rec);
        empty.setRecordOrdinal(3);
        Assert.assertSame(empty.getRecord(), rec);
        Assert.assertEquals(empty.getRecordOrdinal(), 3);
    }

    @Test
    public void testResultStateIsTheDuplicateFlag() {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final SAMRecord rec = records.fragment("r", "rgA", 0, 10, false, "10M", 30);
        final SamRecordWithOrdinalAndSetDuplicateReadFlag wrapped = new SamRecordWithOrdinalAndSetDuplicateReadFlag(rec, 0);
        Assert.assertFalse(rec.getDuplicateReadFlag());
        wrapped.setResultState(true);
        Assert.assertTrue(rec.getDuplicateReadFlag());
        wrapped.setResultState(true);
        Assert.assertTrue(rec.getDuplicateReadFlag());
        wrapped.setResultState(false);
        Assert.assertFalse(rec.getDuplicateReadFlag());

        // a duplicate flag already on the input is cleared when the result is "not a duplicate"
        rec.setDuplicateReadFlag(true);
        wrapped.setResultState(false);
        Assert.assertFalse(rec.getDuplicateReadFlag());
    }

    /** States set out of order; a record is emitted only once it and every earlier record have a state. */
    @Test
    public void testTrackingBufferAppliesStatesInOrder() {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final File tmpDir = IOUtil.createTempDir("trackingBuffer").toFile();
        final SamRecordTrackingBuffer<SamRecordWithOrdinalAndSetDuplicateReadFlag> buffer = new SamRecordTrackingBuffer<>(
                100, 2, Collections.singletonList(tmpDir), records.header, SamRecordWithOrdinalAndSetDuplicateReadFlag.class);
        try {
            final List<SamRecordWithOrdinalAndSetDuplicateReadFlag> wrapped = new ArrayList<>();
            for (int i = 0; i < 5; i++) {
                final SAMRecord rec = records.fragment("r" + i, "rgA", 0, 10 + i, false, "10M", 30);
                // the input flags must not survive
                rec.setDuplicateReadFlag(i % 2 == 0);
                wrapped.add(new SamRecordWithOrdinalAndSetDuplicateReadFlag(rec, i));
                buffer.add(wrapped.get(i));
            }
            Assert.assertFalse(buffer.canEmit());

            buffer.setResultState(wrapped.get(3), true);
            buffer.setResultState(wrapped.get(1), true);
            Assert.assertFalse(buffer.canEmit());
            buffer.setResultState(wrapped.get(0), false);

            final List<String> emitted = new ArrayList<>();
            while (!buffer.isEmpty() && buffer.canEmit()) {
                final SAMRecord rec = buffer.next().getRecord();
                emitted.add(rec.getReadName() + "=" + rec.getDuplicateReadFlag());
            }
            Assert.assertEquals(emitted, List.of("r0=false", "r1=true"));

            buffer.setResultState(wrapped.get(4), false);
            buffer.setResultState(wrapped.get(2), false);
            while (!buffer.isEmpty() && buffer.canEmit()) {
                final SAMRecord rec = buffer.next().getRecord();
                emitted.add(rec.getReadName() + "=" + rec.getDuplicateReadFlag());
            }
            Assert.assertEquals(emitted, List.of("r0=false", "r1=true", "r2=false", "r3=true", "r4=false"));
            Assert.assertTrue(buffer.isEmpty());
        } finally {
            buffer.close();
            IOUtil.deleteDirectoryTree(tmpDir);
        }
    }
}
