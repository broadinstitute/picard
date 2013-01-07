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


import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.net.URL;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Arrays;

import net.sf.samtools.FileTruncatedException;
import net.sf.samtools.seekablestream.SeekableBufferedStream;
import net.sf.samtools.seekablestream.SeekableFileStream;
import net.sf.samtools.seekablestream.SeekableHTTPStream;
import net.sf.samtools.seekablestream.SeekableStream;

/*
 * Utility class for reading BGZF block compressed files.  The caller can treat this file like any other InputStream.
 * It probably is not necessary to wrap this stream in a buffering stream, because there is internal buffering.
 * The advantage of BGZF over conventional GZip format is that BGZF allows for seeking without having to read the
 * entire file up to the location being sought.  Note that seeking is only possible if the ctor(File) is used.
 *
 * c.f. http://samtools.sourceforge.net/SAM1.pdf for details of BGZF format
 */
public class BlockCompressedInputStream extends InputStream {
    private InputStream mStream = null;
    private SeekableStream mFile = null;
    private byte[] mFileBuffer = null;
    private byte[] mCurrentBlock = null;
    private int mCurrentOffset = 0;
    private long mBlockAddress = 0;
    private int mLastBlockLength = 0;
    private final BlockGunzipper blockGunzipper = new BlockGunzipper();


    /**
     * Note that seek() is not supported if this ctor is used.
     */
    public BlockCompressedInputStream(final InputStream stream) {
        mStream = IOUtil.toBufferedStream(stream);
        mFile = null;
    }

    /**
     * Use this ctor if you wish to call seek()
     */
    public BlockCompressedInputStream(final File file)
        throws IOException {
        mFile = new SeekableFileStream(file);
        mStream = null;

    }

    public BlockCompressedInputStream(final URL url) {
        mFile = new SeekableBufferedStream(new SeekableHTTPStream(url));
        mStream = null;
    }

    /**
     * For providing some arbitrary data source.  No additional buffering is
     * provided, so if the underlying source is not buffered, wrap it in a
     * SeekableBufferedStream before passing to this ctor.
     */
    public BlockCompressedInputStream(final SeekableStream strm) {
        mFile = strm;
        mStream = null;
    }

    /**
     * Determines whether or not the inflater will re-calculated the CRC on the decompressed data
     * and check it against the value stored in the GZIP header.  CRC checking is an expensive
     * operation and should be used accordingly.
     */
    public void setCheckCrcs(final boolean check) {
        this.blockGunzipper.setCheckCrcs(check);
    }

    /**
     * @return the number of bytes that can be read (or skipped over) from this input stream without blocking by the
     * next caller of a method for this input stream. The next caller might be the same thread or another thread.
     * Note that although the next caller can read this many bytes without blocking, the available() method call itself
     * may block in order to fill an internal buffer if it has been exhausted.
     */
    public int available()
        throws IOException {
        if (mCurrentBlock == null || mCurrentOffset == mCurrentBlock.length) {
            readBlock();
        }
        if (mCurrentBlock == null) {
            return 0;
        }
        return mCurrentBlock.length - mCurrentOffset;
    }

    /**
     * Closes the underlying InputStream or RandomAccessFile
     */
    public void close()
        throws IOException {
        if (mFile != null) {
            mFile.close();
            mFile = null;
        } else if (mStream != null) {
            mStream.close();
            mStream = null;
        }
        // Encourage garbage collection
        mFileBuffer = null;
        mCurrentBlock = null;
    }

    /**
     * Reads the next byte of data from the input stream. The value byte is returned as an int in the range 0 to 255.
     * If no byte is available because the end of the stream has been reached, the value -1 is returned.
     * This method blocks until input data is available, the end of the stream is detected, or an exception is thrown.

     * @return the next byte of data, or -1 if the end of the stream is reached.
     */
    public int read()
        throws IOException {
        return (available() > 0) ? (mCurrentBlock[mCurrentOffset++] & 0xFF) : -1;
    }

    /**
     * Reads some number of bytes from the input stream and stores them into the buffer array b. The number of bytes
     * actually read is returned as an integer. This method blocks until input data is available, end of file is detected,
     * or an exception is thrown.
     *
     * read(buf) has the same effect as read(buf, 0, buf.length).
     *
     * @param buffer the buffer into which the data is read.
     * @return the total number of bytes read into the buffer, or -1 is there is no more data because the end of
     * the stream has been reached.
     */
    public int read(final byte[] buffer)
        throws IOException {
        return read(buffer, 0, buffer.length);
    }

    private volatile ByteArrayOutputStream buf = null;
    private static final byte eol = '\n';
    private static final byte eolCr = '\r';
    
    /**
     * Reads a whole line. A line is considered to be terminated by either a line feed ('\n'), 
     * carriage return ('\r') or carriage return followed by a line feed ("\r\n").
     *
     * @return  A String containing the contents of the line, excluding the line terminating
     *          character, or null if the end of the stream has been reached
     *
     * @exception  IOException  If an I/O error occurs
     * @
     */
    public String readLine() throws IOException {
    	int available = available();
        if (available == 0) {
            return null;
        }
        if(null == buf){ // lazy initialisation 
        	buf = new ByteArrayOutputStream(8192);
        }
        buf.reset();
    	boolean done = false;
    	boolean foundCr = false; // \r found flag
        while (!done) {
        	int linetmpPos = mCurrentOffset;
        	int bCnt = 0;
        	while((available-- > 0)){
        		final byte c = mCurrentBlock[linetmpPos++];
        		if(c == eol){ // found \n
        			done = true;
        			break;
        		} else if(foundCr){  // previous char was \r
        			--linetmpPos; // current char is not \n so put it back
        			done = true;
        			break;
        		} else if(c == eolCr){ // found \r
					foundCr = true;
        			continue; // no ++bCnt
        		}
				++bCnt;
        	}
        	if(mCurrentOffset < linetmpPos){
				buf.write(mCurrentBlock, mCurrentOffset, bCnt);
	        	mCurrentOffset = linetmpPos;
        	}
        	available = available();    
        	if(available == 0){
        		// EOF
        		done = true;
        	}
        }
    	return buf.toString();
    }

    /**
     * Reads up to len bytes of data from the input stream into an array of bytes. An attempt is made to read
     * as many as len bytes, but a smaller number may be read. The number of bytes actually read is returned as an integer.
     *
     * This method blocks until input data is available, end of file is detected, or an exception is thrown.
     *
     * @param buffer buffer into which data is read.
     * @param offset the start offset in array b  at which the data is written.
     * @param length the maximum number of bytes to read.
     * @return the total number of bytes read into the buffer, or -1 if there is no more data because the end of
     * the stream has been reached.
     */
    public int read(final byte[] buffer, int offset, int length)
        throws IOException {
        final int originalLength = length;
        while (length > 0) {
            final int available = available();
            if (available == 0) {
                // Signal EOF to caller
                if (originalLength == length) {
                    return -1;
                }
                break;
            }
            final int copyLength = Math.min(length, available);
            System.arraycopy(mCurrentBlock, mCurrentOffset, buffer, offset, copyLength);
            mCurrentOffset += copyLength;
            offset += copyLength;
            length -= copyLength;
        }
        return originalLength - length;
    }

    /**
     * Seek to the given position in the file.  Note that pos is a special virtual file pointer,
     * not an actual byte offset.
     *
     * @param pos virtual file pointer
     */
    public void seek(final long pos)
        throws IOException {
        if (mFile == null) {
            throw new IOException("Cannot seek on stream based file");
        }
        // Decode virtual file pointer
        // Upper 48 bits is the byte offset into the compressed stream of a block.
        // Lower 16 bits is the byte offset into the uncompressed stream inside the block.
        final long compressedOffset = BlockCompressedFilePointerUtil.getBlockAddress(pos);
        final int uncompressedOffset = BlockCompressedFilePointerUtil.getBlockOffset(pos);
        final int available;
        if (mBlockAddress == compressedOffset && mCurrentBlock != null) {
            available = mCurrentBlock.length;
        } else {
            mFile.seek(compressedOffset);
            mBlockAddress = compressedOffset;
            mLastBlockLength = 0;
            readBlock();
            available = available();
        }
        if (uncompressedOffset > available ||
                (uncompressedOffset == available && !eof())) {
            throw new IOException("Invalid file pointer: " + pos);
        }
        mCurrentOffset = uncompressedOffset;
    }

    private boolean eof() throws IOException {
        if (mFile.eof()) {
            return true;
        }
        // If the last remaining block is the size of the EMPTY_GZIP_BLOCK, this is the same as being at EOF.
        return (mFile.length() - (mBlockAddress + mLastBlockLength) == BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK.length);
    }

    /**
     * @return virtual file pointer that can be passed to seek() to return to the current position.  This is
     * not an actual byte offset, so arithmetic on file pointers cannot be done to determine the distance between
     * the two.
     */
    public long getFilePointer() {
        if (mCurrentOffset == mCurrentBlock.length) {
            // If current offset is at the end of the current block, file pointer should point
            // to the beginning of the next block.
            return BlockCompressedFilePointerUtil.makeFilePointer(mBlockAddress + mLastBlockLength, 0);
        }
        return BlockCompressedFilePointerUtil.makeFilePointer(mBlockAddress, mCurrentOffset);
    }

    public static long getFileBlock(final long bgzfOffset) {
        return BlockCompressedFilePointerUtil.getBlockAddress(bgzfOffset);
    }
    
    /**
     * @param stream Must be at start of file.  Throws RuntimeException if !stream.markSupported().
     * @return true if the given file looks like a valid BGZF file.
     */
    public static boolean isValidFile(final InputStream stream)
        throws IOException {
        if (!stream.markSupported()) {
            throw new RuntimeException("Cannot test non-buffered stream");
        }
        stream.mark(BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH);
        final byte[] buffer = new byte[BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH];
        final int count = readBytes(stream, buffer, 0, BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH);
        stream.reset();
        return count == BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH && isValidBlockHeader(buffer);
    }

    private static boolean isValidBlockHeader(final byte[] buffer) {
        return (buffer[0] == BlockCompressedStreamConstants.GZIP_ID1 &&
                (buffer[1] & 0xFF) == BlockCompressedStreamConstants.GZIP_ID2 &&
                (buffer[3] & BlockCompressedStreamConstants.GZIP_FLG) != 0 &&
                buffer[10] == BlockCompressedStreamConstants.GZIP_XLEN &&
                buffer[12] == BlockCompressedStreamConstants.BGZF_ID1 &&
                buffer[13] == BlockCompressedStreamConstants.BGZF_ID2);
    }

    private void readBlock()
        throws IOException {

        if (mFileBuffer == null) {
            mFileBuffer = new byte[BlockCompressedStreamConstants.MAX_COMPRESSED_BLOCK_SIZE];
        }
        int count = readBytes(mFileBuffer, 0, BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH);
        if (count == 0) {
            // Handle case where there is no empty gzip block at end.
            mCurrentOffset = 0;
            mBlockAddress += mLastBlockLength;
            mCurrentBlock = new byte[0];
            return;
        }
        if (count != BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH) {
            throw new IOException("Premature end of file");
        }
        final int blockLength = unpackInt16(mFileBuffer, BlockCompressedStreamConstants.BLOCK_LENGTH_OFFSET) + 1;
        if (blockLength < BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH || blockLength > mFileBuffer.length) {
            throw new IOException("Unexpected compressed block length: " + blockLength);
        }
        final int remaining = blockLength - BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH;
        count = readBytes(mFileBuffer, BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH, remaining);
        if (count != remaining) {
            throw new FileTruncatedException("Premature end of file");
        }
        inflateBlock(mFileBuffer, blockLength);
        mCurrentOffset = 0;
        mBlockAddress += mLastBlockLength;
        mLastBlockLength = blockLength;
    }

    private void inflateBlock(final byte[] compressedBlock, final int compressedLength)
        throws IOException {
        final int uncompressedLength = unpackInt32(compressedBlock, compressedLength-4);
        byte[] buffer = mCurrentBlock;
        mCurrentBlock = null;
        if (buffer == null || buffer.length != uncompressedLength) {
            try {
                buffer = new byte[uncompressedLength];
            } catch (NegativeArraySizeException e) {
                throw new RuntimeException("BGZF file has invalid uncompressedLength: " + uncompressedLength, e);
            }
        }
        blockGunzipper.unzipBlock(buffer, compressedBlock, compressedLength);
        mCurrentBlock = buffer;
    }

    private int readBytes(final byte[] buffer, final int offset, final int length)
        throws IOException {
        if (mFile != null) {
            return readBytes(mFile, buffer, offset, length);
        } else if (mStream != null) {
            return readBytes(mStream, buffer, offset, length);
        } else {
            return 0;
        }
    }

    private static int readBytes(final SeekableStream file, final byte[] buffer, final int offset, final int length)
        throws IOException {
        int bytesRead = 0;
        while (bytesRead < length) {
            final int count = file.read(buffer, offset + bytesRead, length - bytesRead);
            if (count <= 0) {
                break;
            }
            bytesRead += count;
        }
        return bytesRead;
    }

    private static int readBytes(final InputStream stream, final byte[] buffer, final int offset, final int length)
        throws IOException {
        int bytesRead = 0;
        while (bytesRead < length) {
            final int count = stream.read(buffer, offset + bytesRead, length - bytesRead);
            if (count <= 0) {
                break;
            }
            bytesRead += count;
        }
        return bytesRead;
    }

    private int unpackInt16(final byte[] buffer, final int offset) {
        return ((buffer[offset] & 0xFF) |
                ((buffer[offset+1] & 0xFF) << 8));
    }

    private int unpackInt32(final byte[] buffer, final int offset) {
        return ((buffer[offset] & 0xFF) |
                ((buffer[offset+1] & 0xFF) << 8) |
                ((buffer[offset+2] & 0xFF) << 16) |
                ((buffer[offset+3] & 0xFF) << 24));
    }

    public enum FileTermination {HAS_TERMINATOR_BLOCK, HAS_HEALTHY_LAST_BLOCK, DEFECTIVE}

    public static FileTermination checkTermination(final File file)
        throws IOException {
        final long fileSize = file.length();
        if (fileSize < BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK.length) {
            return FileTermination.DEFECTIVE;
        }
        final RandomAccessFile raFile = new RandomAccessFile(file, "r");
        try {
            raFile.seek(fileSize - BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK.length);
            byte[] buf = new byte[BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK.length];
            raFile.readFully(buf);
            if (Arrays.equals(buf, BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK)) {
                return FileTermination.HAS_TERMINATOR_BLOCK;
            }
            final int bufsize = (int)Math.min(fileSize, BlockCompressedStreamConstants.MAX_COMPRESSED_BLOCK_SIZE);
            buf = new byte[bufsize];
            raFile.seek(fileSize - bufsize);
            raFile.read(buf);
            for (int i = buf.length - BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK.length;
                    i >= 0; --i) {
                if (!preambleEqual(BlockCompressedStreamConstants.GZIP_BLOCK_PREAMBLE,
                        buf, i, BlockCompressedStreamConstants.GZIP_BLOCK_PREAMBLE.length)) {
                    continue;
                }
                final ByteBuffer byteBuffer = ByteBuffer.wrap(buf, i + BlockCompressedStreamConstants.GZIP_BLOCK_PREAMBLE.length, 4);
                byteBuffer.order(ByteOrder.LITTLE_ENDIAN);
                final int totalBlockSizeMinusOne =  byteBuffer.getShort() & 0xFFFF;
                if (buf.length - i == totalBlockSizeMinusOne + 1) {
                    return FileTermination.HAS_HEALTHY_LAST_BLOCK;
                } else {
                    return FileTermination.DEFECTIVE;
                }
            }
            return FileTermination.DEFECTIVE;
        } finally {
            raFile.close();
        }
    }

    private static boolean preambleEqual(final byte[] preamble, final byte[] buf, final int startOffset, final int length) {
        for (int i = 0; i < length; ++i) {
            if (preamble[i] != buf[i + startOffset]) {
                return false;
            }
        }
        return true;
    }
}


