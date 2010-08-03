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

/**
 * Static for manipulating virtual file pointers in BGZF files.
 * TODO: Make these work properly for very big files, where longs become negative. 
 */
public class BlockCompressedFilePointerUtil {
    private static final int SHIFT_AMOUNT = 16;
    private static final int OFFSET_MASK = 0xffff;
    private static final long ADDRESS_MASK = 0xFFFFFFFFFFFFL;

    /**
     * @param vfp1
     * @param vfp2
     * @return negative if vfp1 is earlier in file than vfp2, positive if it is later, 0 if equal.
     */
    public static int compare(final long vfp1, final long vfp2) {
        if (vfp1 < vfp2) return -1;
        if (vfp1 > vfp2) return 1;
        return 0;
    }

    /**
     * @param blockAddress File offset of start of BGZF block.
     * @param blockOffset Offset into uncompressed block.
     * @return Virtual file pointer that embodies the input parameters.
     */
    static long makeFilePointer(final long blockAddress, final int blockOffset) {
        return blockAddress << SHIFT_AMOUNT | blockOffset;
    }

    /**
     * @param virtualFilePointer
     * @return File offset of start of BGZF block for this virtual file pointer.
     */
    static long getBlockAddress(final long virtualFilePointer) {
        return (virtualFilePointer >> SHIFT_AMOUNT) & ADDRESS_MASK;
    }

    /**
     * @param virtualFilePointer
     * @return Offset into uncompressed block for this virtual file pointer.
     */
    static int getBlockOffset(final long virtualFilePointer) {
        return (int) (virtualFilePointer & 0xFFFF);
    }
}
