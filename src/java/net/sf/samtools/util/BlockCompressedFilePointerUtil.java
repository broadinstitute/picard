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
 */
public class BlockCompressedFilePointerUtil {
    private static final int SHIFT_AMOUNT = 16;
    private static final int OFFSET_MASK = 0xffff;
    private static final long ADDRESS_MASK = 0xFFFFFFFFFFFFL;

    public static final long MAX_BLOCK_ADDRESS = ADDRESS_MASK;
    public static final int MAX_OFFSET = OFFSET_MASK;
    
    /**
     * @param vfp1
     * @param vfp2
     * @return negative if vfp1 is earlier in file than vfp2, positive if it is later, 0 if equal.
     */
    public static int compare(final long vfp1, final long vfp2) {
        if (vfp1 == vfp2) return 0;
        // When treating as unsigned, negative number is > positive.
        if (vfp1 < 0 && vfp2 >= 0) return 1;
        if (vfp1 >= 0 && vfp2 < 0) return -1;
        // Either both negative or both non-negative, so regular comparison works.
        if (vfp1 < vfp2) return -1;
        return 1; // vfp1 > vfp2
    }

    /**
     * @return true if vfp2 points to somewhere in the same BGZF block, or the one immediately following vfp1's BGZF block.
     */
    public static boolean areInSameOrAdjacentBlocks(final long vfp1, final long vfp2) {
        final long block1 = getBlockAddress(vfp1);
        final long block2 = getBlockAddress(vfp2);
        return (block1 == block2 || block1 + 1 == block2);        
    }

    /**
     * @param blockAddress File offset of start of BGZF block.
     * @param blockOffset Offset into uncompressed block.
     * @return Virtual file pointer that embodies the input parameters.
     */
    static long makeFilePointer(final long blockAddress, final int blockOffset) {
        if (blockOffset < 0) {
            throw new IllegalArgumentException("Negative blockOffset " + blockOffset + " not allowed.");
        }
        if (blockAddress < 0) {
            throw new IllegalArgumentException("Negative blockAddress " + blockAddress + " not allowed.");
        }
        if (blockOffset > MAX_OFFSET) {
            throw new IllegalArgumentException("blockOffset " + blockOffset + " too large.");
        }
        if (blockAddress > MAX_BLOCK_ADDRESS) {
            throw new IllegalArgumentException("blockAddress " + blockAddress + " too large.");
        }
        return blockAddress << SHIFT_AMOUNT | blockOffset;
    }

    /**
     * @param virtualFilePointer
     * @return File offset of start of BGZF block for this virtual file pointer.
     */
    public static long getBlockAddress(final long virtualFilePointer) {
        return (virtualFilePointer >> SHIFT_AMOUNT) & ADDRESS_MASK;
    }

    /**
     * @param virtualFilePointer
     * @return Offset into uncompressed block for this virtual file pointer.
     */
    public static int getBlockOffset(final long virtualFilePointer) {
        return (int) (virtualFilePointer & OFFSET_MASK);
    }

    public static String asString(final long vfp) {
        return String.format("%d(0x%x): (block address: %d, offset: %d)", vfp, vfp, getBlockAddress(vfp), getBlockOffset(vfp));
    }
}
