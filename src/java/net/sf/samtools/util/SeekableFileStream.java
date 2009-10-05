/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.sf.samtools.util;

import java.io.*;

/**
 *
 * @author jrobinso
 */
public class SeekableFileStream extends SeekableStream {

    File file;
    RandomAccessFile fis;

    public SeekableFileStream(final File file) throws FileNotFoundException {
        this.file = file;
        fis = new RandomAccessFile(file, "r");
    }

    public long length() {
        return file.length();
    }

    public boolean eof() throws IOException {
        return fis.length() == fis.getFilePointer();
    }

    public void seek(final long position) throws IOException {
        fis.seek(position);
    }
    
    public int read(final byte[] buffer, final int offset, final int length) throws IOException {
        if (length < 0) {
            throw new IndexOutOfBoundsException();
        }
        int n = 0;
        while (n < length) {
            final int count = fis.read(buffer, offset + n, length - n);
            if (count < 0) {
                return n;
            }
            n += count;
        }
        return n;

    }


    public void close() throws IOException {
        fis.close();

    }

    public byte[] readBytes(final long position, final int nBytes) throws IOException {
        seek(position);

        final byte[] buffer = new byte[nBytes];
        final int bytesRead = read(buffer, 0, nBytes);
        if (bytesRead != nBytes) {
            throw new EOFException("Trying to read " + nBytes + " from " + file + " at position " + position +
            ", but only read " + bytesRead + " bytes.");
        }
        return buffer;
    }

    public int read() throws IOException {
        return fis.read();  
    }
}