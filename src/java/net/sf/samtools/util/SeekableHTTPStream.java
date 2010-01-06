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

    public static int readCount = 0;

    public int read(byte[] buffer, int offset, int len) throws IOException {

        if (offset < 0 || len < 0 || (offset + len) > buffer.length) {
            throw new IndexOutOfBoundsException();
        }
        if (len == 0) {
            return 0;
        }

        readCount++;

        HttpURLConnection connection = null;
        InputStream is = null;
        String byteRange = "";
        int n = 0;
        try {
            connection = (HttpURLConnection) url.openConnection();
            byteRange = "bytes=" + position + "-" + (position + len - 1);
            connection.setRequestProperty("Range", byteRange);
            is = connection.getInputStream();

            while (n < len) {
                int count = is.read(buffer, offset + n, len - n);
                if (count < 0) {
                    return (n == 0 ? -1 : n);
                }
                n += count;
            }

            position += n;

            return n;

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


    public int read() throws IOException {
        throw new UnsupportedOperationException("read() not support for SeekableHTTPStreams");
    }
}