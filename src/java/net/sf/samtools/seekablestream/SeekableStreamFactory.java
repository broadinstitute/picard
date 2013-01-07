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

import java.io.File;
import java.io.IOException;
import java.net.URL;

/**
 * @author jrobinso
 * @date Nov 30, 2009
 */
public class SeekableStreamFactory {

    public static SeekableStream getStreamFor(URL url) throws IOException{
        return getStreamFor(url.toExternalForm());
    }

    public static SeekableStream getStreamFor(String path) throws IOException {
        // todo -- add support for SeekableBlockInputStream

        if (path.startsWith("http:") || path.startsWith("https:")) {
            final URL url = new URL(path);
            return new SeekableHTTPStream(url);
        } else if (path.startsWith("ftp:")) {
            return new SeekableFTPStream(new URL(path));
        } else {
            return new SeekableFileStream(new File(path));
        }
    }

}
