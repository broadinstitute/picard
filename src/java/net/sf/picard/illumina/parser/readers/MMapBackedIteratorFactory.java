/*
 * The MIT License
 *
 * Copyright (c) 2012 The Broad Institute
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
package net.sf.picard.illumina.parser.readers;


import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;
import net.sf.samtools.util.CloserUtil;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.*;
import java.nio.channels.FileChannel;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * MMapBackedIteratorFactory a file reader that takes a header size and a binary file, maps the file to
 * a read-only byte buffer and provides methods to retrieve the header as it's own bytebuffer and create
 * iterators of different data types over the values of file (starting after the end of the header).
 * Values provided by the MMappedBinaryFileReader are read as if they are little endian.
 *
 * Note (read to end):
 * This class IS thread-safe and immutable though the iterator and ByteBuffers it produces are NOT.
 * The values read are assumed to be signed, NO promoting/sign conversion happens in this class.
 */
public class MMapBackedIteratorFactory {
    private static int BYTE_SIZE  = 1;
    private static int INT_SIZE   = 4;
    private static int FLOAT_SIZE = 4;

    public static BinaryFileIterator<Integer> getIntegerIterator(final int headerSize, final File binaryFile) {
        checkFactoryVars(headerSize, binaryFile);
        final ByteBuffer buf = getBuffer(binaryFile);
        final byte [] header = getHeader(buf, headerSize);

        return new IntegerMMapIterator(header, binaryFile, buf);
    }

    public static BinaryFileIterator<Byte> getByteIterator(final int headerSize, final File binaryFile) {
        checkFactoryVars(headerSize, binaryFile);
        final ByteBuffer buf = getBuffer(binaryFile);
        final byte [] header = getHeader(buf, headerSize);

        return new ByteMMapIterator(header, binaryFile, buf);
    }

    public static BinaryFileIterator<Float> getFloatIterator(final int headerSize, final File binaryFile) {
        checkFactoryVars(headerSize, binaryFile);
        final ByteBuffer buf = getBuffer(binaryFile);
        final byte [] header = getHeader(buf, headerSize);

        return new FloatMMapIterator(header, binaryFile, buf);
    }

    private static void checkFactoryVars(final int headerSize, final File binaryFile) {
        IoUtil.assertFileIsReadable(binaryFile);

        if(headerSize < 0) {
            throw new PicardException("Header size cannot be negative.  HeaderSize(" + headerSize + ")");
        }

        if(headerSize > binaryFile.length()) {
            throw new PicardException("Header size(" + headerSize + ") is greater than file size(" + binaryFile.length() + ")");
        }
    }

    private static ByteBuffer getBuffer(final File binaryFile) {
        final ByteBuffer buf;
        try {
            final FileInputStream is = new FileInputStream(binaryFile);
            final FileChannel channel = is.getChannel();
            final long fileSize = channel.size();
            buf = channel.map(FileChannel.MapMode.READ_ONLY, 0, fileSize);
            buf.order(ByteOrder.LITTLE_ENDIAN);
            CloserUtil.close(channel);
            CloserUtil.close(is);
        } catch (IOException e) {
            throw new PicardException("IOException opening cluster binary file " + binaryFile, e);
        }

        return buf;
    }

    private static byte [] getHeader(final ByteBuffer buf, final int headerSize) {
        final byte [] headerBytes = new byte[headerSize];
        if(headerSize > 0) {
            buf.get(headerBytes);
        }
        return headerBytes;
    }

    /** A simple iterator that uses a reference to the enclosing ByteBuffer and a member position
     * value to iterate over values in the buffer, starting after headerSize bytes */
    static abstract class MMapBackedIterator<TYPE> extends BinaryFileIterator<TYPE>{
        protected final ByteBuffer buffer;

        protected MMapBackedIterator(final byte[] header, final File file, final int elementSize, final ByteBuffer buffer) {
            super(header, file, elementSize);
            this.buffer = buffer;
        }

        public boolean hasNext() {
            return buffer.limit() - buffer.position() >= elementSize;
        }

        /** The method that actually retrieves the data from the enclosing buffer */
        protected abstract TYPE getElement();

        public Iterator<TYPE> iterator() {
            return this;
        }
    }

    private static class IntegerMMapIterator extends MMapBackedIterator<Integer> {
        public IntegerMMapIterator(final byte[] header, final File file, final ByteBuffer buf) {
            super(header, file, INT_SIZE, buf);
        }

        @Override
        protected Integer getElement() {
            return buffer.getInt();
        }
    }

    private static class ByteMMapIterator extends MMapBackedIterator<Byte> {
        public ByteMMapIterator(final byte[] header, final File file, final ByteBuffer buf) {
            super(header, file, BYTE_SIZE, buf);
        }

        @Override
        protected Byte getElement() {
            return buffer.get();
        }
    }

    private static class FloatMMapIterator extends MMapBackedIterator<Float> {
        public FloatMMapIterator(final byte[] header, final File file, final ByteBuffer buf) {
            super(header, file, FLOAT_SIZE, buf);
        }

        @Override
        protected Float getElement() {
            return buffer.getFloat();
        }
    }
}


abstract class BinaryFileIterator<TYPE> implements Iterator<TYPE>, Iterable<TYPE> {
    protected final File file;
    protected final long fileSize;
    protected final int elementSize;
    private final byte [] header;

    public BinaryFileIterator(final byte[] header, final File file, final int elementSize) {
        this.header = header;
        this.file   = file;
        this.fileSize = file.length();
        this.elementSize = elementSize;
    }
    /** Return the bytes found in the first headerSize bytes of the file, wrapped as a
     * ByteBuffer */
    public ByteBuffer getHeaderBytes() {
        final ByteBuffer bb = ByteBuffer.allocate(header.length);
        bb.order(ByteOrder.LITTLE_ENDIAN);
        bb.put(header);
        bb.position(0);
        return bb;
    }

    public void assertTotalElementsEqual(final long numElements) {
        if(getElementsInFile() != numElements) {
            throw new PicardException("Expected " + numElements + " elements in file but found " + getElementsInFile() + " elements! File(" + file.getAbsolutePath() +  ")");
        }

        if(getExtraBytes() != 0) {
            throw new PicardException("Malformed file, expected " + (header.length + numElements * elementSize) + " bytes in file, found " + fileSize + " bytes for file("
                    + file.getAbsolutePath() + ")");
        }
    }

    public int getElementSize() {
        return elementSize;
    }

    public long getExtraBytes() {
        return fileSize - header.length - (getElementsInFile() * elementSize);
    }

    public long getElementsInFile() {
        return (fileSize - header.length) / elementSize;
    }

    public File getFile() {
        return file;
    }


    public TYPE next() {
        if(!hasNext()) {
            throw new NoSuchElementException();
        }
        return getElement();
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }


    public Iterator<TYPE> iterator() {
        return this;
    }

    /** The method that actually retrieves the data from the enclosing buffer */
    protected abstract TYPE getElement();
    public abstract boolean hasNext();
}



