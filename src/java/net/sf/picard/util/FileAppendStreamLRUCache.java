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
package net.sf.picard.util;

import net.sf.picard.PicardException;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.FileNotFoundException;

/**
 * LRU cache of FileOutputStreams to handle situation in which it is necessary to have more FileOuputStreams
 * than resource limits will allow.  Least-recently-used FileOutputStream is closed when it is pushed out of
 * the cache.  When adding a new element to the cache, the file is opened in append mode. 
 *
 * @author alecw@broadinstitute.org
 */
public class FileAppendStreamLRUCache extends ResourceLimitedMap<File, FileOutputStream> {
    public FileAppendStreamLRUCache(final int cacheSize) {
        super(cacheSize, new Functor());
    }

    private static class Functor implements ResourceLimitedMapFunctor<File, FileOutputStream> {

        // Explicitly GC after this many calls to close() in order to force file handles to truly be released.

        private static final int GC_FREQUENCY = 10000;
        private int numCloses = 0;

        public FileOutputStream makeValue(final File file) {
            try {
                return new FileOutputStream(file, true);
            } catch (FileNotFoundException e) {
                // In case the file could not be opened because of too many file handles, try to force
                // file handles to be closed.
                System.gc();
                System.runFinalization();
                try {
                    return new FileOutputStream(file, true);
                } catch (FileNotFoundException e2) {
                    throw new PicardException(file + "not found", e2);
                }
            }
        }

        public void finalizeValue(final File file, final FileOutputStream fileOutputStream) {
            try {
                fileOutputStream.close();
            } catch (IOException e) {
                throw new PicardException("Exception closing FileOutputStream for " + file, e);
            }
        }
    }
}
