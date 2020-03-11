package picard.util;

import picard.PicardException;

/**
 * Implementation of a circular byte buffer that uses a large byte[] internally and supports basic
 * read/write operations from/to other byte[]s passed as arguments. Uses wait/nofity() to manage
 * cross-thread coordination when the buffer is either full or empty.
 */
public class CircularByteBuffer {
    private final byte[] bytes;
    private final int capacity;

    private int nextWritePos = 0;
    private int bytesAvailableToWrite;
    private int nextReadPos = 0;
    private int bytesAvailableToRead = 0;

    private boolean closed = false;

    /** Constructs a buffer capable of holding the given number of bytes. */
    public CircularByteBuffer(final int size) {
        this.bytes = new byte[size];
        this.capacity = this.bytes.length;
        this.bytesAvailableToWrite = this.capacity;
    }

    /**
     * Write bytes into the buffer from the supplied array.  Will attempt to read 'size' bytes beginning
     * at 'start' in the supplied array and copy them into the buffer.  If the buffer is near full or
     * cannot write 'size' bytes contiguously it may write fewer than 'size' bytes.
     *
     * @return the number of bytes read from the input array and copied into the buffer
     */
    synchronized public int write(final byte[] bytes, final int start, final int size) {
        if (closed) throw new IllegalStateException("Cannot write to closed buffer.");

        try { if (this.bytesAvailableToWrite == 0) wait(); }
        catch (final InterruptedException ie) {throw new PicardException("Interrupted while waiting to write to fifo.", ie); }

        final int writePos      = this.nextWritePos;
        final int distanceToEnd = this.capacity - writePos;
        final int available     = distanceToEnd < this.bytesAvailableToWrite ? distanceToEnd : this.bytesAvailableToWrite;
        final int length        = available < size ? available : size;

        System.arraycopy(bytes, start, this.bytes, writePos, length);
        this.bytesAvailableToWrite -= length;
        this.bytesAvailableToRead  += length;
        this.nextWritePos = (writePos + length) % this.capacity;
        notify();
        return length;
    }

    /**
     * Read bytes from the buffer into the supplied array.  Will attempt to read 'size' bytes and write
     * them into the supplied array beginning at index 'start' in the supplied array.  If the buffer is
     * near empty or cannot read 'size' bytes contiguously it may write fewer than 'size' bytes.
     *
     * @return the number of bytes read from the buffer and copied into the input array
     */
    synchronized public int read(final byte[] bytes, final int start, final int size) {
        try { if (this.bytesAvailableToRead == 0 && !closed) wait(); }
        catch (final InterruptedException ie) {throw new PicardException("Interrupted while waiting to read from fifo.", ie); }

        final int readPos       = this.nextReadPos;
        final int distanceToEnd = this.capacity - readPos;
        final int available     = distanceToEnd < this.bytesAvailableToRead ? distanceToEnd : this.bytesAvailableToRead;
        final int length        = available < size ? available : size;

        System.arraycopy(this.bytes, readPos, bytes, start, length);
        this.bytesAvailableToRead  -= length;
        this.bytesAvailableToWrite += length;
        this.nextReadPos = (readPos + length) % this.capacity;
        notify();
        return length;
    }

    /** Signals that the buffer is closed and no further writes will occur. */
    synchronized public void close() {
        this.closed = true;
        notify();
    }

    /** Returns true if the buffer is closed, false otherwise. */
    synchronized public boolean isClosed() {
        return this.closed;
    }

    /** Returns the total capacity of the buffer (empty+filled). */
    public int getCapacity() { return this.capacity; }

    /** Returns the number of bytes that are in the buffer at the time of the method invocation. */
    synchronized public int getBytesAvailableToRead() { return this.bytesAvailableToRead; }
}