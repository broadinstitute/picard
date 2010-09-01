package net.sf.samtools;

import net.sf.samtools.util.BlockCompressedFilePointerUtil;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * A [start,stop) file pointer pairing into the BAM file, stored
 * as a BAM file index.  A chunk is represented as a single 64-bit
 * value where the high-order 48 bits point to the location of the
 * start of a compressed BGZF block within a BGZF file and the
 * low-order 16 bits point to a position within the decompressed
 * data in the BGZF block.
 *
 * See the SAM/BAM spec for more details.
 */
class Chunk implements Cloneable, Serializable,Comparable<Chunk> {
    private static final long serialVersionUID = 1L;

    /**
     * A pointer to the start of a region in a SAM/BAM file.  The
     * start is inclusive: start reading from this point.
     */
    private long mChunkStart;

    /**
     * A pointer to the end of a region in a SAM/BAM file.  The end
     * is exclusive: this pointer points to one byte past the end
     * of the region of interest inside the file.
     */
    private long mChunkEnd;

    public Chunk(final long start, final long end) {
        mChunkStart = start;
        mChunkEnd = end;
    }

    public Chunk clone() {
        return new Chunk(mChunkStart,mChunkEnd);
    }

    protected long getChunkStart() {
        return mChunkStart;
    }

    protected void setChunkStart(final long value) {
        mChunkStart = value;
    }

    protected long getChunkEnd() {
        return mChunkEnd;
    }

    protected void setChunkEnd(final long value) {
        mChunkEnd = value;
    }

    public int compareTo(final Chunk chunk) {
        int result = Long.signum(mChunkStart - chunk.mChunkStart);
        if (result == 0) {
            result = Long.signum(mChunkEnd - chunk.mChunkEnd);
        }
        return result;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final Chunk chunk = (Chunk) o;

        if (mChunkEnd != chunk.mChunkEnd) return false;
        if (mChunkStart != chunk.mChunkStart) return false;

        return true;
    }

    /**
     * Returns whether two chunks overlap.
     * @param other Chunk to which this should be compared.
     * @return True if the chunks overlap.  Returns false if the two chunks abut or are disjoint.
     */
    public boolean overlaps(final Chunk other) {
        int comparison = this.compareTo(other);
        if(comparison == 0)
            return true;

        // "sort" the two chunks using the comparator.
        Chunk leftMost = comparison==-1 ? this : other;
        Chunk rightMost = comparison==1 ? this : other;

        long leftMostBlockAddress = BlockCompressedFilePointerUtil.getBlockAddress(leftMost.getChunkEnd());
        long rightMostBlockAddress = BlockCompressedFilePointerUtil.getBlockAddress(rightMost.getChunkStart());

        // If the left block's address is after the right block's address, compare the two blocks.
        // If the two blocks are identical, compare the block offsets.
        // If the right block is after the left block, no overlap is possible.
        if(leftMostBlockAddress > rightMostBlockAddress)
            return true;
        else if(leftMostBlockAddress == rightMostBlockAddress) {
            int leftMostOffset = BlockCompressedFilePointerUtil.getBlockOffset(leftMost.getChunkEnd());
            int rightMostOffset = BlockCompressedFilePointerUtil.getBlockOffset(rightMost.getChunkStart());
            return leftMostOffset > rightMostOffset;
        }
        else
            return false;
    }

    /**
     * Returns whether two chunks overlap.
     * @param other Chunk to which this should be compared.
     * @return True if the two chunks are adjacent.  Returns false if the chunks overlap or are discontinuous.
     */
    public boolean isAdjacentTo(final Chunk other) {
        // Simpler implementation would be to == the chunk end of one to the chunk start of the other.  Chose this implementation to ensure that all chunk
        // comparisons point directly to the 
        return (BlockCompressedFilePointerUtil.getBlockAddress(this.getChunkEnd()) == BlockCompressedFilePointerUtil.getBlockAddress(other.getChunkStart()) &&
                BlockCompressedFilePointerUtil.getBlockOffset(this.getChunkEnd()) == BlockCompressedFilePointerUtil.getBlockOffset(other.getChunkStart())) ||
               (BlockCompressedFilePointerUtil.getBlockAddress(this.getChunkStart()) == BlockCompressedFilePointerUtil.getBlockAddress(other.getChunkEnd()) &&
                BlockCompressedFilePointerUtil.getBlockOffset(this.getChunkStart()) == BlockCompressedFilePointerUtil.getBlockOffset(other.getChunkEnd()));
    }

    @Override
    public int hashCode() {
        int result = (int) (mChunkStart ^ (mChunkStart >>> 32));
        result = 31 * result + (int) (mChunkEnd ^ (mChunkEnd >>> 32));
        return result;
    }

    @Override
    public String toString() {
        return String.format("%d:%d-%d:%d",mChunkStart >> 16,mChunkStart & 0xFFFF,mChunkEnd >> 16,mChunkEnd & 0xFFFF);
    }
}
