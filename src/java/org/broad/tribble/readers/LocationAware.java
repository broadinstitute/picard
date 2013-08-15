package org.broad.tribble.readers;

/**
 * Describes API for getting current position in a stream, writer, or underlying file.
 * 
 * The expected functionality is simple: if you are a output stream / writer, and you've written 50 bytes to the stream, then 
 * {@link #getPosition()} should return 50; if you are an input stream or file reader, and you've read 25 bytes from the object, then it 
 * should return 25.
 * 
 * In the context of an iterator or any producer-like object that doesn't map directly to a byte stream, {@link #getPosition()} should
 * return the position (in the underlying stream being read/written to) of the most-recently read/written element.  For example, if you
 * are reading lines from a file with a {@link AsciiLineReaderIterator}, calling {@link #getPosition()} should return the byte position 
 * of the start of the most recent line returned by {@link org.broad.tribble.readers.AsciiLineReaderIterator#next()}.
 * 
 * @author mccowan
 */
public interface LocationAware {
    /**
     * The current offset, in bytes, of this stream/writer/file.  Or, if this is an iterator/producer, the offset (in bytes) of the
     * END of the most recently returned record (since a produced record corresponds to something that has been read already). See class
     * javadoc for more.
     */
    public long getPosition();
}
