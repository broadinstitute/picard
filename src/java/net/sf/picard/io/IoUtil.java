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
package net.sf.picard.io;

import net.sf.picard.PicardException;

import java.io.*;
import java.util.Arrays;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * A class for utility methods that wrap or aggregate functionality in Java IO.
 *
 * @author Tim Fennell
 */
public class IoUtil {
    /**
     * Checks that a file is non-null, exists, is not a directory and is readable.  If any
     * condition is false then a runtime exception is thrown.
     *
     * @param file the file to check for readability
     */
    public static void assertFileIsReadable(final File file) {
        if (file == null) {
			throw new IllegalArgumentException("Cannot check readability of null file.");
		} else if (!file.exists()) {
            throw new PicardException("Cannot read non-existent file: " + file.getAbsolutePath());
        }
        else if (file.isDirectory()) {
            throw new PicardException("Cannot read file because it is a directory: " + file.getAbsolutePath());
        }
        else if (!file.canRead()) {
            throw new PicardException("File exists but is not readable: " + file.getAbsolutePath());
        }
    }

    /**
     * Checks that a file is non-null, and is either extent and writable, or non-existent but
     * that the parent directory exists and is writable. If any
     * condition is false then a runtime exception is thrown.
     *
     * @param file the file to check for writability
     */
    public static void assertFileIsWritable(final File file) {
        if (file == null) {
			throw new IllegalArgumentException("Cannot check readability of null file.");
		} else if (!file.exists()) {
            // If the file doesn't exist, check that it's parent directory does and is writable
            final File parent = file.getAbsoluteFile().getParentFile();
            if (!parent.exists()) {
                throw new PicardException("Cannot write file: " + file.getAbsolutePath() + ". " +
                        "Neither file nor parent directory exist.");
            }
            else if (!parent.isDirectory()) {
                throw new PicardException("Cannot write file: " + file.getAbsolutePath() + ". " +
                        "File does not exist and parent is not a directory.");
            }
            else if (!parent.canWrite()) {
                throw new PicardException("Cannot write file: " + file.getAbsolutePath() + ". " +
                        "File does not exist and parent directory is not writable..");
            }
        }
        else if (file.isDirectory()) {
            throw new PicardException("Cannot write file because it is a directory: " + file.getAbsolutePath());
        }
        else if (!file.canWrite()) {
            throw new PicardException("File exists but is not writable: " + file.getAbsolutePath());
        }
    }

    /**
     * Checks that a directory is non-null, extent, writable and a directory 
     * otherwise a runtime exception is thrown.
     *
     * @param dir the dir to check for writability
     */
    public static void assertDirectoryIsWritable(final File dir) {
        if (dir == null) {
            throw new IllegalArgumentException("Cannot check readability of null file.");
        } 
        else if (!dir.exists()) {
            throw new PicardException("Directory does not exist: " + dir.getAbsolutePath());
        }
        else if (!dir.isDirectory()) {
            throw new PicardException("Cannot write to directory because it is not a directory: " + dir.getAbsolutePath());
        }
        else if (!dir.canWrite()) {
            throw new PicardException("Directory exists but is not writable: " + dir.getAbsolutePath());
        }
    }

    /**
     * Checks that the two files are the same length, and have the same content, otherwise throws a runtime exception.
     */
    public static void assertFilesEqual(final File f1, final File f2) {
        try {
            if (f1.length() != f2.length()) {
                throw new PicardException("Files " + f1 + " and " + f2 + " are different lengths.");
            }
            final FileInputStream s1 = new FileInputStream(f1);
            final FileInputStream s2 = new FileInputStream(f2);
            final byte[] buf1 = new byte[1024 * 1024];
            final byte[] buf2 = new byte[1024 * 1024];
            int len1;
            while ((len1 = s1.read(buf1)) != -1) {
                final int len2 = s2.read(buf2);
                if (len1 != len2) {
                    throw new PicardException("Unexpected EOF comparing files that are supposed to be the same length.");
                }
                if (!Arrays.equals(buf1, buf2)) {
                    throw new PicardException("Files " + f1 + " and " + f2 + " differ.");
                }
            }
            s1.close();
            s2.close();
        } catch (IOException e) {
            throw new PicardException("Exception comparing files " + f1 + " and " + f2, e);
        }

    }

    /**
     * Checks that a file is of non-zero length
     */
    public static void assertFileSizeNonZero(final File file) {
        if (file.length() == 0) {
            throw new PicardException(file.getAbsolutePath() + " has length 0");
        }
    }

    /**
     * Opens a file for reading, decompressing it if necessary
     *
     * @param file  The file to open
     * @return the input stream to read from
     */
    public static InputStream openFileForReading(final File file) {

        try {
            if (file.getName().endsWith(".gz") ||
                file.getName().endsWith(".bfq"))  {
                return openGzipFileForReading(file);
            }
            //TODO: Other compression formats
            else {
                return new FileInputStream(file);
            }
        }
        catch (IOException ioe) {
            throw new PicardException("Error opening file: " + file.getName(), ioe);
        }

    }

    /**
     * Opens a GZIP-encoded file for reading, decompressing it if necessary
     *
     * @param file  The file to open
     * @return the input stream to read from
     */
    public static InputStream openGzipFileForReading(final File file) {

        try {
            return new GZIPInputStream(new FileInputStream(file));
        }
        catch (IOException ioe) {
            throw new PicardException("Error opening file: " + file.getName(), ioe);
        }
   }

    /**
     * Opens a file for writing, overwriting the file if it already exists
     *
     * @param file  the file to write to
     * @return the output stream to write to
     */
    public static OutputStream openFileForWriting(final File file) {
        return openFileForWriting(file, false);
    }

    /**
     * Opens a file for writing
     *
     * @param file  the file to write to
     * @param append    whether to append to the file if it already exists (we overwrite it if false)
     * @return the output stream to write to
     */
    public static OutputStream openFileForWriting(final File file, final boolean append) {

        try {
            if (file.getName().endsWith(".gz") ||
                file.getName().endsWith(".bfq")) {
                return openGzipFileForWriting(file, append);
            }
            //TODO: Other compression formats
            else {
                return new FileOutputStream(file, append);
            }
        }
        catch (IOException ioe) {
            throw new PicardException("Error opening file for writing: " + file.getName(), ioe);
        }
    }

    /**
     * Preferred over PrintStream and PrintWriter because an exception is thrown on I/O error
     */
    public static BufferedWriter openFileForBufferedWriting(final File file, final boolean append) {
        try {
            return new BufferedWriter(new FileWriter(file, append));
        } catch (IOException ioe) {
            throw new PicardException("Error opening file for writing: " + file.getName(), ioe);
        }
    }

    /**
     * Preferred over PrintStream and PrintWriter because an exception is thrown on I/O error
     */
    public static BufferedWriter openFileForBufferedWriting(final File file) {
        try {
            return new BufferedWriter(new FileWriter(file, false));
        } catch (IOException ioe) {
            throw new PicardException("Error opening file for writing: " + file.getName(), ioe);
        }
    }

    /**
     * Opens a GZIP encoded file for writing
     *
     * @param file  the file to write to
     * @param append    whether to append to the file if it already exists (we overwrite it if false)
     * @return the output stream to write to
     */
    public static OutputStream openGzipFileForWriting(final File file, final boolean append) {

        try {
            return new GZIPOutputStream(new FileOutputStream(file, append));
        }
        catch (IOException ioe) {
            throw new PicardException("Error opening file for writing: " + file.getName(), ioe);
        }
    }

    /**
     * Utility method to copy the contents of input to output. The caller is responsible for
     * opening and closing both streams.
     * 
     * @param input contents to be copied
     * @param output destination
     */
    public static void copyStream(final InputStream input, final OutputStream output) {
        try {
            final byte[] buffer = new byte[1024];
            int bytesRead = 0;
            while((bytesRead = input.read(buffer)) > 0) {
                output.write(buffer, 0, bytesRead);
            }
        } catch (IOException e) {
            throw new PicardException("Exception copying stream", e);
        }
    }

    /**
     * Copy input to output, overwriting output if it already exists.
     */
    public static void copyFile(final File input, final File output) {
        try {
            final InputStream is = new FileInputStream(input);
            final OutputStream os = new FileOutputStream(output);
            copyStream(is, os);
            os.close();
            is.close();
        } catch (IOException e) {
            throw new PicardException("Error copying " + input + " to " + output, e);
        }
    }

    /**
     * 
     * @param directory
     * @param regexp
     * @return
     */
    public static File[] getFilesMatchingRegexp(final File directory, final String regexp) {
        final Pattern pattern = Pattern.compile(regexp);
        return getFilesMatchingRegexp(directory, pattern);
    }

    public static File[] getFilesMatchingRegexp(final File directory, final Pattern regexp) {
        return directory.listFiles( new FilenameFilter() {
            public boolean accept(final File dir, final String name) {
                return regexp.matcher(name).matches();
            }
        });
    }

    /**
     * Delete the given file or directory.  If a directory, all enclosing files and subdirs are also deleted.
     */
    public static void deleteDirectoryTree(final File fileOrDirectory) {
        if (fileOrDirectory.isDirectory()) {
            for (final File child : fileOrDirectory.listFiles()) {
                deleteDirectoryTree(child);
            }
        }
        net.sf.samtools.util.IOUtil.deleteFiles(fileOrDirectory);
    }

    /**
     * Create a temporary subdirectory in the default temporary-file directory, using the given prefix and suffix to generate the name.
     * Note that this method is not completely safe, because it create a temporary file, deletes it, and then creates
     * a directory with the same name as the file.  Should be good enough.
     *
     * @param prefix The prefix string to be used in generating the file's name; must be at least three characters long
     * @param suffix The suffix string to be used in generating the file's name; may be null, in which case the suffix ".tmp" will be used
     * @return File object for new directory
     */
    public static File createTempDir(final String prefix, final String suffix) {
        try {
            final File tmp = File.createTempFile(prefix, suffix);
            if (!tmp.delete()) {
                throw new PicardException("Could not delete temporary file " + tmp);
            }
            if (!tmp.mkdir()) {
                throw new PicardException("Could not create temporary directory " + tmp);
            }
            return tmp;
        } catch (IOException e) {
            throw new PicardException("Exception creating temporary directory.", e);
        }
    }

    /** Checks that a file exists and is readable, and then returns a buffered reader for it. */
    public static BufferedReader openFileForBufferedReading(final File file) throws IOException {
        return new BufferedReader(new InputStreamReader(openFileForReading(file)));
	}

    /** Takes a string and replaces any characters that are not safe for filenames with an underscore */
    public static String makeFileNameSafe(String str) {
        return str.trim().replaceAll("[\\s!\"#$%&'()*/:;<=>?@\\[\\]\\\\^`{|}~]", "_");
    }

}
