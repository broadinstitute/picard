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


import net.sf.samtools.Defaults;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.File;

/**
 * Miscellaneous stateless static IO-oriented methods.
 */
public class IOUtil {
    /**
     * @deprecated Use Defaults.BUFFER_SIZE instead.
     */
    @Deprecated public static final int STANDARD_BUFFER_SIZE = Defaults.BUFFER_SIZE;

    public static final long ONE_GB   = 1024 * 1024 * 1024;
    public static final long TWO_GBS  = 2 * ONE_GB;
    public static final long FIVE_GBS = 5 * ONE_GB;

    /**
     * Wrap the given stream in a BufferedInputStream, if it isn't already wrapper
     * @param stream stream to be wrapped
     * @return A BufferedInputStream wrapping stream, or stream itself if stream instanceof BufferedInputStream.
     */
    public static BufferedInputStream toBufferedStream(final InputStream stream) {
        if (stream instanceof BufferedInputStream) {
            return (BufferedInputStream) stream;
        } else {
            return new BufferedInputStream(stream, STANDARD_BUFFER_SIZE);
        }
    }

    /**
     * Delete a list of files, and write a warning message if one could not be deleted.
     * @param files Files to be deleted.
     */
    public static void deleteFiles(final File... files) {
        for (final File f : files) {
            if (!f.delete()) {
                System.err.println("Could not delete file " + f);
            }
        }
    }

    public static void deleteFiles(final Iterable<File> files) {
        for (final File f : files) {
            if (!f.delete()) {
                System.err.println("Could not delete file " + f);
            }
        }
    }


    /**
     * @return true if the path is not a device (e.g. /dev/null or /dev/stdin), and is not
     * an existing directory.  I.e. is is a regular path that may correspond to an existing
     * file, or a path that could be a regular output file.
     */
    public static boolean isRegularPath(final File file) {
        return !file.exists() || file.isFile();
    }

    /**
     * Creates a new tmp file on one of the available temp filesystems, registers it for deletion
     * on JVM exit and then returns it.
     */
    public static File newTempFile(final String prefix, final String suffix,
                                   final File[] tmpDirs, final long minBytesFree) throws IOException {
        File f = null;

        for (int i=0; i<tmpDirs.length; ++i) {
            if (tmpDirs[i].getUsableSpace() > minBytesFree || i == tmpDirs.length-1) {
                f = File.createTempFile(prefix, suffix, tmpDirs[i]);
                f.deleteOnExit();
                break;
            }
        }

        return f;
    }

    /** Creates a new tmp file on one of the potential filesystems that has at least 5GB free. */
    public static File newTempFile(final String prefix, final String suffix,
                                   final File[] tmpDirs) throws IOException {
        return newTempFile(prefix, suffix, tmpDirs, FIVE_GBS);
    }


    /** Returns a default tmp directory. */
    public static File getDefaultTmpDir() {
        final String user = System.getProperty("user.name");
        final String tmp = System.getProperty("java.io.tmpdir");

        if (tmp.endsWith("/" + user)) return new File(tmp);
        else return new File(tmp, user);
    }
}
