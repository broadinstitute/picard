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
package picard.util;

import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import picard.PicardException;

import java.io.File;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * TextFileParser which reads a single text file.
 *
 * @author Kathleen Tibbetts
 */
public class BasicInputParser extends AbstractInputParser
{
    private BufferedLineReader reader;
    private final ArrayList<InputStream> inputs = new ArrayList<InputStream>();
    private final ArrayList<String> fileNames = new ArrayList<String>();
    String currentFileName = null;
    private String currentLine = null;
    private String nextLine = null;
    private int currentLineNumber = 0;
    private int nextLineNumber = 0;

    /**
     * Constructor.  Opens up a buffered reader and reads the first line.
     *
     * @param inputStreams  the file(s) to parse, in order
     */
    public BasicInputParser(final boolean treatGroupedDelimitersAsOne, final InputStream... inputStreams) {
        if (inputStreams.length == 0) {
            throw new IllegalArgumentException("At least one input must be specified.");
        }
        this.inputs.addAll(Arrays.asList(inputStreams));
        reader = new BufferedLineReader(this.inputs.remove(0));
        this.setTreatGroupedDelimitersAsOne(treatGroupedDelimitersAsOne);
    }

    public BasicInputParser(final boolean treatGroupedDelimitersAsOne, final int wordCount, final InputStream... inputStreams) {
        this(treatGroupedDelimitersAsOne, inputStreams);
        setWordCount(wordCount);
    }

    /**
     * Constructor.  Opens up a buffered reader and reads the first line.
     *
     * @param files  the file(s) to parse, in order
     */
    public BasicInputParser(final boolean treatGroupedDelimitersAsOne, final File... files) {
        this(treatGroupedDelimitersAsOne, filesToInputStreams(files));
        for (File f : files) fileNames.add(f.getAbsolutePath());
        this.currentFileName = fileNames.remove(0);
    }

    /**
     * Constructor.  In addition to opening and priming the files, it sets the number of
     * whitespace-separated "words" per line.
     *
     * @param files      the file(s) to parse
     * @param wordCount number of whitespace-separated "words" per line
     */
    public BasicInputParser(final boolean treatGroupedDelimitersAsOne, final int wordCount, final File... files) {
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
            if (nextLine != null && !isComment(nextLine.getBytes())) {
                currentLineNumber = nextLineNumber;
                currentLine = nextLine;
            }
            if (line != null) {
                nextLineNumber++;
                nextLine = line;
                return line.getBytes();
            }
            if (!inputs.isEmpty()) {
                advanceFile();
                return readNextLine();
            }
            return null;
        }
        catch(RuntimeIOException ioe) {
            throw new PicardException("Error reading from file " + currentFileName, ioe);
        }
    }

    protected void advanceFile() {
        currentFileName = !fileNames.isEmpty() ? fileNames.remove(0) : null;
        nextLineNumber = 0;
        nextLine = null;
        reader = new BufferedLineReader(inputs.remove(0));
    }

    /**
     * Closes the underlying stream
     */
    public void close() {
        if (reader != null)  {
            reader.close();
        }
        for(final InputStream stream : inputs){
            CloserUtil.close(stream);
        }
    }

    /**
     * Gets the name of the file being parsed
     *
     * @return  the name of the file being parsed
     */
    public String getFileName() {
        return this.currentFileName != null ? this.currentFileName : "(file name unavailable)";
    }

    /**
     * Provides access to the current (just parsed) line in pre-parsed format.
     * NOTE: Because AbstractInputParser pre-fetches the next line, this method actually returns the
     * next line, not the most recent line returned by next().
     */
    public String getCurrentLine() {
        return this.currentLine;
    }

    /**
     * NOTE: Because AbstractInputParser pre-fetches the next line, this method actually returns the
     * next line, not the most recent line returned by next().
     */
    public int getCurrentLineNumber() {
        return currentLineNumber;
    }

    private static InputStream[] filesToInputStreams(final File[] files) {
        final InputStream[] result = new InputStream[files.length];
        for (int i = 0; i < files.length; i++) {
            result[i] = IOUtil.openFileForReading(files[i]);
        }
        return result;
    }
}
