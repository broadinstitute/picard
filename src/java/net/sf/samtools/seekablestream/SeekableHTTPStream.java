/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
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
package net.sf.samtools.seekablestream;

import net.sf.samtools.util.HttpUtils;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.Proxy;
import java.net.URL;

/**
 * @author jrobinso
 */
public class SeekableHTTPStream extends SeekableStream {

    private long position = 0;
    private long contentLength = -1;
    private final URL url;
    private final Proxy proxy;

    public SeekableHTTPStream(final URL url) {
        this(url, null);

    }

    public SeekableHTTPStream(final URL url, Proxy proxy) {

        this.proxy = proxy;
        this.url = url;

        // Try to get the file length
        final String contentLengthString = HttpUtils.getHeaderField(url, "Content-Length");
        if (contentLengthString != null) {
            try {
                contentLength = Long.parseLong(contentLengthString);
            }
            catch (NumberFormatException ignored) {
                System.err.println("WARNING: Invalid content length (" + contentLengthString + "  for: " + url);
                contentLength = -1;
            }
        }

    }

    public long position() {
        return position;
    }

    public long length() {
        return contentLength;
    }

    @Override
    public long skip(long n) throws IOException {
        long bytesToSkip = Math.min(n, contentLength - position);
        position += bytesToSkip;
        return bytesToSkip;
    }

    public boolean eof() throws IOException {
        return contentLength > 0 && position >= contentLength;
    }

    public void seek(final long position) {
        this.position = position;
    }

    public int read(byte[] buffer, int offset, int len) throws IOException {

        if (offset < 0 || len < 0 || (offset + len) > buffer.length) {
            throw new IndexOutOfBoundsException("Offset="+offset+",len="+len+",buflen="+buffer.length);
        }
        if (len == 0) {
            return 0;
        }

        HttpURLConnection connection = null;
        InputStream is = null;
        String byteRange = "";
        int n = 0;
        try {
            connection = proxy == null ?
                    (HttpURLConnection) url.openConnection() :
                    (HttpURLConnection) url.openConnection(proxy);

            long endRange = position + len - 1;
            // IF we know the total content length, limit the end range to that.
            if (contentLength > 0) {
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
                if (n == 0) {
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
    	byte []tmp=new byte[1];
    	read(tmp,0,1);
    	return (int) tmp[0] & 0xFF; 
    }

    @Override
    public String getSource() {
        return url.toString();
    }
}
