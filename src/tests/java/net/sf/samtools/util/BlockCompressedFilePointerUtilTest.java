/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.util.ArrayList;
import java.util.List;


public class BlockCompressedFilePointerUtilTest
{
    @Test
    public void basicTest() 
    {
        List<Long> pointers = new ArrayList<Long>();
        pointers.add(makeFilePointer(0, 0));
        pointers.add(makeFilePointer(0, BlockCompressedFilePointerUtil.MAX_OFFSET));
        final long BIG_BLOCK_ADDRESS = 1L << 46;
        pointers.add(makeFilePointer(BIG_BLOCK_ADDRESS-1, 0));
        pointers.add(makeFilePointer(BIG_BLOCK_ADDRESS-1, BlockCompressedFilePointerUtil.MAX_OFFSET));
        pointers.add(makeFilePointer(BIG_BLOCK_ADDRESS, 0));
        pointers.add(makeFilePointer(BIG_BLOCK_ADDRESS, BlockCompressedFilePointerUtil.MAX_OFFSET));
        pointers.add(makeFilePointer(BlockCompressedFilePointerUtil.MAX_BLOCK_ADDRESS, 0));
        pointers.add(makeFilePointer(BlockCompressedFilePointerUtil.MAX_BLOCK_ADDRESS, BlockCompressedFilePointerUtil.MAX_OFFSET));
        for (int i = 0; i < pointers.size() - 1; ++i) {
            for (int j = i+1; j < pointers.size(); ++j) {
                Assert.assertTrue(BlockCompressedFilePointerUtil.compare(pointers.get(i), pointers.get(j)) < 0,
                        BlockCompressedFilePointerUtil.asString(pointers.get(i)) + " should be < " +
                                BlockCompressedFilePointerUtil.asString(pointers.get(j)));
                Assert.assertTrue(BlockCompressedFilePointerUtil.compare(pointers.get(j), pointers.get(i)) > 0,
                        BlockCompressedFilePointerUtil.asString(pointers.get(j)) + " should be > " +
                                BlockCompressedFilePointerUtil.asString(pointers.get(i)));
            }
        }
        
    }

    /**
     * Create the virtual file pointer, and also assert that is can be converted back into the input parameters.
     * @param blockAddress
     * @param blockOffset
     * @return block compressed file pointer
     */
    private long makeFilePointer(long blockAddress, int blockOffset)
    {
        final long ret = BlockCompressedFilePointerUtil.makeFilePointer(blockAddress, blockOffset);
        Assert.assertEquals(BlockCompressedFilePointerUtil.getBlockAddress(ret), blockAddress);
        Assert.assertEquals(BlockCompressedFilePointerUtil.getBlockOffset(ret), blockOffset);
        Assert.assertEquals(BlockCompressedFilePointerUtil.compare(ret, ret), 0);
        return ret;
    }

    @Test(dataProvider = "badInputs", expectedExceptions = IllegalArgumentException.class)
    public void negativeTests(long blockAddress, int blockOffset) {
        BlockCompressedFilePointerUtil.makeFilePointer(blockAddress, blockOffset);
        Assert.assertFalse(true, "Should not get here.");
    }

    @DataProvider(name="badInputs")
    public Object[][]  badInputs() {
        return new Object[][]{
                {-1L, 0},
                {0L, -1},
                {BlockCompressedFilePointerUtil.MAX_BLOCK_ADDRESS+1, 0},
                {0L, BlockCompressedFilePointerUtil.MAX_OFFSET+1}
        };
    }
}

/******************************************************************/
/**************************[END OF BlockCompressedFilePointerUtilTest.java]*************************/
/******************************************************************/
