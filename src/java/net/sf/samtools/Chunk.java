package net.sf.samtools;

import java.io.Serializable;

/**
 * A [start,stop) file pointer pairing into the BAM file, stored
 * as a BAM file index.
 *
 * @author mhanna
 * @version 0.1
 */
class Chunk implements Cloneable, Serializable,Comparable<Chunk> {
    private static final long serialVersionUID = 1L;
    
    private long mChunkStart;
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
