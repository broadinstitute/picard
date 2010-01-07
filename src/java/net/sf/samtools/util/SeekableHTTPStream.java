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
        if (contentLengthString != null) {
            try {
                contentLength = Long.parseLong(contentLengthString);
            }
            catch (NumberFormatException ignored) {
                System.out.println("WARNING: Invalid content length (" + contentLengthString + "  for: " + url);
                contentLength = -1;
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

    public int read(byte[] buffer, int offset, int len) throws IOException {

        if (offset < 0 || len < 0 || (offset + len) > buffer.length) {
            throw new IndexOutOfBoundsException();
        }
        if (len == 0) {
            return 0;
        }

        HttpURLConnection connection = null;
        InputStream is = null;
        String byteRange = "";
        int n = 0;
        try {
            connection = (HttpURLConnection) url.openConnection();

            long endRange = position + len - 1;
            // IF we know the total content length, limit the end range to that.
            if(contentLength > 0) {
                endRange = Math.min(endRange, contentLength);
            }
            byteRange = "bytes=" + position + "-" + endRange;
            connection.setRequestProperty("Range", byteRange);
            is = connection.getInputStream();

            while (n < len) {
                int count = is.read(buffer, offset + n, len - n);
                if (count < 0) {
                    if (n == 0) {
                        return -1;
                    } else {
                        break;
                    }
                }
                n += count;
            }

            position += n;

            return n;

        }

        catch (IOException e) {
            // THis is a bit of a hack, but its not clear how else to handle this.  If a byte range is specified
            // that goes past the end of the file the response code will be 416.  The MAC os translates this to
            // an IOException with the 416 code in the message.  Windows translates the error to an EOFException.
            //
            //  The BAM file iterator  uses the return value to detect end of file (specifically looks for n == 0).
            if (e.getMessage().contains("416") || (e instanceof EOFException)) {
                if (n < 0) {
                    return -1;
                } else {
                    position += n;
                    // As we are at EOF, the contentLength and position are by definition =
                    contentLength = position;
                    return n;
                }
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


    public int read() throws IOException {
        throw new UnsupportedOperationException("read() not support for SeekableHTTPStreams");
    }
}