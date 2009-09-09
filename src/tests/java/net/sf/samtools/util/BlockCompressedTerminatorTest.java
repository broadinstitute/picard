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
package net.sf.samtools.util;

import org.testng.annotations.Test;
import org.testng.Assert;

import java.io.File;

/**
 * @author alecw@broadinstitute.org
 */
public class BlockCompressedTerminatorTest {
    private static final File TEST_DATA_DIR = new File("testdata/net/sf/samtools/util");

    @Test
    public void testFileWithTerminator() throws Exception {
        final File tmpCompressedFile = File.createTempFile("test.", ".bgzf");
        tmpCompressedFile.deleteOnExit();
        final BlockCompressedOutputStream os = new BlockCompressedOutputStream(tmpCompressedFile);
        os.write("Hi, Mom!\n".getBytes());
        os.close();
        Assert.assertEquals(BlockCompressedInputStream.checkTermination(tmpCompressedFile),
                BlockCompressedInputStream.FileTermination.HAS_TERMINATOR_BLOCK);
    }

    @Test
    public void testValidFileWithoutTerminator() throws Exception {
        Assert.assertEquals(BlockCompressedInputStream.checkTermination(new File(TEST_DATA_DIR, "no_bgzf_terminator.bam")),
                BlockCompressedInputStream.FileTermination.HAS_HEALTHY_LAST_BLOCK);
    }

    @Test
    public void testDefectiveFile() throws Exception {
        Assert.assertEquals(BlockCompressedInputStream.checkTermination(new File(TEST_DATA_DIR, "defective_bgzf.bam")),
                BlockCompressedInputStream.FileTermination.DEFECTIVE);
    }
}
