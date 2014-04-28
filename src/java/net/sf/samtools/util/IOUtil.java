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
import net.sf.samtools.seekablestream.SeekableBufferedStream;
import net.sf.samtools.seekablestream.SeekableFileStream;
import net.sf.samtools.seekablestream.SeekableStream;

import java.io.*;

/**
 * Miscellaneous stateless static IO-oriented methods.
 */
public class IOUtil {
    /**
     * @deprecated Use Defaults.NON_ZERO_BUFFER_SIZE instead.
     */
    @Deprecated public static final int STANDARD_BUFFER_SIZE = Defaults.NON_ZERO_BUFFER_SIZE;

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
            return new BufferedInputStream(stream, Defaults.NON_ZERO_BUFFER_SIZE);
        }
    }

    /**
     * Transfers from the input stream to the output stream using stream operations and a buffer.
     */
    public static void transferByStream(final InputStream in, final OutputStream out, final long bytes) {
        final byte[] buffer = new byte[Defaults.NON_ZERO_BUFFER_SIZE];
        long remaining = bytes;

        try {
            while (remaining > 0) {
                final int read = in.read(buffer, 0, (int) Math.min(buffer.length, remaining));
                out.write(buffer, 0, read);
                remaining -= read;
            }
        }
        catch (final IOException ioe) {
            throw new RuntimeIOException(ioe);
        }
    }

    /**
     * @return If Defaults.BUFFER_SIZE > 0, wrap os in BufferedOutputStream, else return os itself.
     */
    public static OutputStream maybeBufferOutputStream(final OutputStream os) {
        return maybeBufferOutputStream(os, Defaults.BUFFER_SIZE);
    }

    /**
     * @return If bufferSize > 0, wrap os in BufferedOutputStream, else return os itself.
     */
    public static OutputStream maybeBufferOutputStream(final OutputStream os, final int bufferSize) {
        if (bufferSize > 0) return new BufferedOutputStream(os, bufferSize);
        else return os;
    }

    public static SeekableStream maybeBufferedSeekableStream(final SeekableStream stream, final int bufferSize) {
        return bufferSize > 0 ? new SeekableBufferedStream(stream, bufferSize) : stream; 
    }
    
    public static SeekableStream maybeBufferedSeekableStream(final SeekableStream stream) {
        return maybeBufferedSeekableStream(stream, Defaults.BUFFER_SIZE);
    }
    
    public static SeekableStream maybeBufferedSeekableStream(final File file) {
        try {
            return maybeBufferedSeekableStream(new SeekableFileStream(file));
        } catch (final FileNotFoundException e) {
            throw new RuntimeIOException(e);
        }
    }
    
    /**
     * @return If Defaults.BUFFER_SIZE > 0, wrap is in BufferedInputStream, else return is itself.
     */
    public static InputStream maybeBufferInputStream(final InputStream is) {
        return maybeBufferInputStream(is, Defaults.BUFFER_SIZE);
    }
    
    /**
     * @return If bufferSize > 0, wrap is in BufferedInputStream, else return is itself.
     */
    public static InputStream maybeBufferInputStream(final InputStream is, final int bufferSize) {
        if (bufferSize > 0) return new BufferedInputStream(is, bufferSize);
        else return is;
    }

    public static Reader maybeBufferReader(Reader reader, final int bufferSize) {
        if (bufferSize > 0) reader = new BufferedReader(reader, bufferSize);
        return reader;
    }

    public static Reader maybeBufferReader(final Reader reader) {
        return maybeBufferReader(reader, Defaults.BUFFER_SIZE);
    }

    public static Writer maybeBufferWriter(Writer writer, final int bufferSize) {
        if (bufferSize > 0) writer = new BufferedWriter(writer, bufferSize);
        return writer;
    }

    public static Writer maybeBufferWriter(final Writer writer) {
        return maybeBufferWriter(writer, Defaults.BUFFER_SIZE);
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
            if ( i == tmpDirs.length-1 || tmpDirs[i].getUsableSpace() > minBytesFree) {
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

        if (tmp.endsWith(File.separatorChar + user)) return new File(tmp);
        else return new File(tmp, user);
    }

    /** Returns the name of the file minus the extension (i.e. text after the last "." in the filename). */
    public static String basename(final File f) {
        final String full = f.getName();
        final int index = full.lastIndexOf(".");
        if (index > 0  && index > full.lastIndexOf(File.separator)) {
            return full.substring(0, index);
        }
        else {
            return full;
        }
    }
}
