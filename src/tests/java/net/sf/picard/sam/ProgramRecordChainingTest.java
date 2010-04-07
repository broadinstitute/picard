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

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Test for SequenceUtil.chainProgramRecord
 */
public class ProgramRecordChainingTest {

    @Test
    public void testChainProgramRecord() {
        SAMFileHeader header = new SAMFileHeader();
        SAMProgramRecord first = header.createProgramRecord();
        SAMUtils.chainSAMProgramRecord(header, first);
        Assert.assertEquals(header.getProgramRecords().size(), 1);
        Assert.assertNull(first.getPreviousProgramGroupId());

        SAMProgramRecord second = header.createProgramRecord();
        SAMUtils.chainSAMProgramRecord(header, second);
        Assert.assertEquals(header.getProgramRecords().size(), 2);
        Assert.assertNull(first.getPreviousProgramGroupId());
        Assert.assertEquals(second.getPreviousProgramGroupId(), first.getProgramGroupId());

        SAMProgramRecord third = new SAMProgramRecord("3");
        SAMUtils.chainSAMProgramRecord(header, third);
        header.addProgramRecord(third);
        Assert.assertEquals(header.getProgramRecords().size(), 3);
        Assert.assertNull(first.getPreviousProgramGroupId());
        Assert.assertEquals(second.getPreviousProgramGroupId(), first.getProgramGroupId());
        Assert.assertEquals(third.getPreviousProgramGroupId(), second.getProgramGroupId());

    }
}
