/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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
import net.sf.picard.io.IoUtil;

import java.io.*;
import java.util.Arrays;
import java.util.ArrayList;

import net.sf.samtools.util.AsciiLineReader;
import net.sf.samtools.util.RuntimeIOException;

/**
 * TextFileParser which reads a single text file.
 *
 * @author Kathleen Tibbetts
 */
public class BasicTextFileParser extends AbstractTextFileParser
{
    private AsciiLineReader reader;
    private final ArrayList<File> files = new ArrayList<File>();
    String currentFileName = null;
    private String currentLine;
    private int currentLineNumber;

    /**
     * Constructor.  Opens up a buffered reader and reads the first line.
     *
     * @param files  the file(s) to parse, in order
     */
    public BasicTextFileParser(final boolean treatGroupedDelimitersAsOne, final File... files) {
        if (files.length == 0) {
            throw new IllegalArgumentException("At least one file must be specified.");
        }
        this.files.addAll(Arrays.asList(files));
        final File f = this.files.remove(0);
        currentFileName = f.getAbsolutePath();
        reader = new AsciiLineReader(IoUtil.openFileForReading(f));
        this.setTreatGroupedDelimitersAsOne(treatGroupedDelimitersAsOne);
    }

    /**
     * Constructor.  In addition to opening and priming the files, it sets the number of
     * whitespace-separated "words" per line.
     *
     * @param files      the file(s) to parse
     * @param wordCount number of whitespace-separated "words" per line
     */
    public BasicTextFileParser(final boolean treatGroupedDelimitersAsOne, final int wordCount, final File... files) {
        this(treatGroupedDelimitersAsOne, files);
        setWordCount(wordCount);
    }
    /**
     * Workhorse method that reads the next line from the underlying reader
     *
     * @return  String or null if there is no next line
     */
    protected byte[] readNextLine()
    {
        try {
            final String line = reader.readLine();
            if (line != null) {
                currentLineNumber++;
                currentLine = line;
                return line.getBytes();
            }
            if (files.size() > 0) {
                currentFileName = files.get(0).getAbsolutePath();
                currentLineNumber = 0;
                reader = new AsciiLineReader(IoUtil.openFileForReading(files.remove(0)));
                return readNextLine();
            }
            return null;
        }
        catch(RuntimeIOException ioe) {
            throw new PicardException("Error reading from file " + currentFileName, ioe);
        }
    }

    /**
     * Closes the underlying stream
     */
    public void close() {
        if (reader != null)  {
            reader.close();
        }
    }

    /**
     * Gets the name of the file being parsed
     *
     * @return  the name of the file being parsed
     */
    public String getFileName() {
        return this.currentFileName;
    }

    /**
     * Provides access to the current (just parsed) line in pre-parsed format.
     * NOTE: Because AbstractTextFileParser pre-fetches the next line, this method actually returns the
     * next line, not the most recent line returned by next().
     */
    public String getCurrentLine() {
        return this.currentLine;
    }

    /**
     * NOTE: Because AbstractTextFileParser pre-fetches the next line, this method actually returns the
     * next line, not the most recent line returned by next().
     */
    public int getCurrentLineNumber() {
        return currentLineNumber;
    }
}
