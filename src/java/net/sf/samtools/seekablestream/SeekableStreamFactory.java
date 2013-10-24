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
 * Singleton class for getting {@link SeekableStream}s from URL/paths
 * Applications using this library can set their own factory
 * @author jrobinso
 * @date Nov 30, 2009
 */
public class SeekableStreamFactory{

    private static final ISeekableStreamFactory DEFAULT_FACTORY;
    private static ISeekableStreamFactory currentFactory;

    static{
        DEFAULT_FACTORY = new DefaultSeekableStreamFactory();
        currentFactory = DEFAULT_FACTORY;
    }

    private SeekableStreamFactory(){}

    public static void setInstance(final ISeekableStreamFactory factory){
        currentFactory = factory;
    }

    public static ISeekableStreamFactory getInstance(){
        return currentFactory;
    }

    /**
     * Does this path point to a regular file on disk and not something like a URL?
     * @param path the path to test
     * @return true if the path is to a file on disk
     */
    public static boolean isFilePath(final String path) {
        return ! ( path.startsWith("http:") || path.startsWith("https:") || path.startsWith("ftp:") );
    }

    private static class DefaultSeekableStreamFactory implements ISeekableStreamFactory {

        public SeekableStream getStreamFor(final URL url) throws IOException {
            return getStreamFor(url.toExternalForm());
        }

        public SeekableStream getStreamFor(final String path) throws IOException {
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

}
