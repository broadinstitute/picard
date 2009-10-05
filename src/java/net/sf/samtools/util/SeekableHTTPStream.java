/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.sf.samtools.util;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.URL;

/**
 * @author jrobinso
 */
public class SeekableHTTPStream extends SeekableStream {

    private long position = 0;
    private long contentLength = -1;
    private final URL url;

    public SeekableHTTPStream(final URL url) {
        this.url = url;

        // Try to get the file length
        final String contentLengthString = HttpUtils.getHeaderField(url, "Content-Length");
        if(contentLengthString != null) {
            try {
                contentLength = Long.parseLong(contentLengthString);
            }
            catch(NumberFormatException ignored) {

            }
        }

    }

    public long length() {
        return contentLength;
    }

    public boolean eof() throws IOException {
        return position >= contentLength;
    }

    public void seek(final long position) {
        this.position = position;
    }

    public int read(final byte[] buffer, final int offset, final int length) throws IOException {
        if (offset < 0 || length < 0 || (offset + length) > buffer.length) {
            throw new IndexOutOfBoundsException();
        }

        HttpURLConnection connection = null;
        InputStream is = null;

        int n = 0;
        try {

            connection = (HttpURLConnection) url.openConnection();
            final String byteRange = "bytes=" + position + "-" + (position + length - 1);
            connection.setRequestProperty("Range", byteRange);
            is = connection.getInputStream();

            while (n < length) {
                final int count = is.read(buffer, offset + n, length - n);
                if (count < 0) {
                    throw new EOFException();
                }
                n += count;
            }

            position += n;
            
            return n;

        } catch (IOException e) {
            // THis is a bit of a hack, but its not clear how else to handle this.  If a byte range is specified
            // that goes past the end of the file the response code will be 416.  The MAC os translates this to
            // an IOException with the 416 code in the message.  Windows translates the error to an EOFException.
            //
            //  The BAM file iterator  uses the return value to detect end of file (specifically looks for n == 0).
            if (e.getMessage().contains("416") || (e instanceof EOFException)) {
                return n;
            } else {
                throw e;
            }
        }

        finally {
            if (is != null) {
                is.close();
            }
            if (connection != null) {
                connection.disconnect();
            }
        }
    }


    public void close() throws IOException {
        // Nothing to do
    }


    public byte[] readBytes(final long position, final int nBytes) throws IOException {
        this.position = position;
        final byte[] buffer = new byte[nBytes];
        final int bytesRead = read(buffer, 0, nBytes);
        if (bytesRead != nBytes) {
            throw new EOFException("Trying to read " + nBytes + " from " + url + " at position " + position +
            ", but only read " + bytesRead + " bytes.");
        }
        return buffer;
    }

    public int read() throws IOException {
        throw new UnsupportedOperationException("read() not support for SeekableHTTPStreams");
    }
}