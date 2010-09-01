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
package net.sf.samtools;

import org.testng.Assert;
import org.testng.annotations.Test;

public class ChunkTest {
    @Test
    public void testOverlaps() {
        // Test completely disjoint offsets.
        Assert.assertFalse(new Chunk(1,5).overlaps(new Chunk(11,15)),"Test found overlap in non-overlapping offsets");
        Assert.assertFalse(new Chunk(11,15).overlaps(new Chunk(1,5)),"Test found overlap in non-overlapping offsets");

        // Test adjacent offsets
        Assert.assertFalse(new Chunk(1,5).overlaps(new Chunk(6,10)),"Test found overlap in adjacent offsets");
        Assert.assertFalse(new Chunk(6,10).overlaps(new Chunk(1,5)),"Test found overlap in adjacent offsets");

        // Test overlapping offsets
        Assert.assertTrue(new Chunk(1,5).overlaps(new Chunk(2,6)),"Test returned incorrect value for overlapping offsets");
        Assert.assertTrue(new Chunk(2,6).overlaps(new Chunk(1,5)),"Test returned incorrect value for overlapping offsets");
        Assert.assertTrue(new Chunk(1,5).overlaps(new Chunk(2,6)),"Test returned incorrect value for overlapping offsets");

        // Completely disjoint and adjacent blocks
        Assert.assertFalse(new Chunk(1<<16,2<<16).overlaps(new Chunk(3<<16,4<<16)),"Test found overlap in non-overlapping blocks");
        Assert.assertFalse(new Chunk(1<<16,2<<16).overlaps(new Chunk(2<<16,3<<16)),"test found overlap in adjacent blocks");

        // True overlaps in the same block
        Assert.assertTrue(new Chunk(1<<16,2<<16).overlaps(new Chunk(1<<16,2<<16)),"Test failed to find overlap in completely overlapping chunks");
        Assert.assertTrue(new Chunk(1<<16,2<<16).overlaps(new Chunk(1<<16,1<<16|0xFF00)),"Test failed to find overlap in chunks with head overlapping");
        Assert.assertTrue(new Chunk(1<<16,2<<16).overlaps(new Chunk(1<<16|0xFF00,2<<16)),"Test failed to find overlap in overlapping chunks");
        Assert.assertTrue(new Chunk(1<<16,2<<16).overlaps(new Chunk(1<<16|0xFEFF,1<<16|0xFF00)),"Test failed to find overlap in contained chunk");
        Assert.assertTrue(new Chunk(1<<16|0xFEFF,1<<16|0xFF00).overlaps(new Chunk(1<<16,2<<16)),"Test failed to find overlap in enclosing chunk");
        Assert.assertTrue(new Chunk(1<<16,1<<16|0xFF00).overlaps(new Chunk(1<<16|0xFEFF,2<<16)),"Test failed to find tail->head overlap");
        Assert.assertTrue(new Chunk(1<<16|0xFEFF,2<<16).overlaps(new Chunk(1<<16,1<<16|0xFF00)),"Test failed to find head->tail overlap");

        // Test overlaps spanning blocks
        Assert.assertTrue(new Chunk(1<<16,2<<16).overlaps(new Chunk(1<<16|0xFF00,2<<16|0xFF00)),"Test failed to find overlap spanning blocks");
    }

    @Test
    public void testAdjacency() {
        // Test offset adjacency
        Assert.assertTrue(new Chunk(1,5).isAdjacentTo(new Chunk(5,9)),"Offsets which should be adjacent are not");
        Assert.assertTrue(new Chunk(5,9).isAdjacentTo(new Chunk(1,5)),"Offsets which should be adjacent are not");

        // Test block adjacency
        Assert.assertTrue(new Chunk(1<<16,2<<16).isAdjacentTo(new Chunk(2<<16,3<<16)),"Blocks which should be adjacent are not");

        // Test offset-block adjacency
        Assert.assertTrue(new Chunk(1<<16,2<<16|0xFF00).isAdjacentTo(new Chunk(2<<16|0xFF00,3<<16)),"Block-offsets which should be adjacent are not");

        // Discontiguous and overlapping blocks
        Assert.assertFalse(new Chunk(1,5).isAdjacentTo(new Chunk(11,15)),"Disjoint block should not be adjacent");
        Assert.assertFalse(new Chunk(1,5).isAdjacentTo(new Chunk(2,3)),"Contained offset should not be adjacent");
    }
}
