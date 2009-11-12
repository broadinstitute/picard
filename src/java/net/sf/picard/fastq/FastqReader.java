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
package net.sf.picard.fastq;

import net.sf.picard.PicardException;
import net.sf.samtools.util.StringUtil;
import net.sf.picard.io.IoUtil;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.io.*;

/**
 * Reads a fastq file.
 */
public class FastqReader implements Iterator<FastqRecord>, Iterable<FastqRecord>, Closeable {
    final private File fastqFile;
    final private BufferedReader reader;
    private FastqRecord nextRecord;
    private int line=1; 

    public FastqReader(final File file) {
        fastqFile = file;
        reader = new BufferedReader(new InputStreamReader(IoUtil.openFileForReading(fastqFile)));
        nextRecord = readNextRecord();
    }

    private FastqRecord readNextRecord() {
        try {

            // Read sequence header
            final String seqHeader = reader.readLine();
            if (seqHeader == null) return null ;
            if (StringUtil.isBlank(seqHeader)) {
                throw new PicardException(error("Missing sequence header"));
            }
            if (!seqHeader.startsWith(FastqConstants.SEQUENCE_HEADER)) {
                throw new PicardException(error("Sequence header must start with "+ FastqConstants.SEQUENCE_HEADER+": "+seqHeader));
            }

            // Read sequence line
            final String seqLine = reader.readLine();
            checkLine(seqLine,"sequence line");

            // Read quality header
            final String qualHeader = reader.readLine();
            checkLine(qualHeader,"quality header");
            if (!qualHeader.startsWith(FastqConstants.QUALITY_HEADER)) {
                throw new PicardException(error("Quality header must start with "+ FastqConstants.QUALITY_HEADER+": "+qualHeader));
            }

            // Read quality line
            final String qualLine = reader.readLine();
            checkLine(qualLine,"quality line");

            // Check sequence and quality lines are same length
            if (seqLine.length() != qualLine.length()) {
                throw new PicardException(error("Sequence and quality line must be the same length"));
            }

            final FastqRecord frec = new FastqRecord(seqHeader.substring(1, seqHeader.length()), seqLine,
                    qualHeader.substring(1, qualHeader.length()), qualLine);
            line += 4 ;
            return frec ;

        } catch (IOException e) {
            throw new PicardException(String.format("Error reading '%s'", fastqFile.getAbsolutePath()),e);
        }
    }

    public boolean hasNext() { return nextRecord != null; }

    public FastqRecord next() {
        if (!hasNext()) {
            throw new NoSuchElementException("next() called when !hasNext()");
        }
        final FastqRecord rec = nextRecord;
        nextRecord = readNextRecord();
        return rec;
    }

    public void remove() { throw new UnsupportedOperationException("Unsupported operation"); }

    public Iterator<FastqRecord> iterator() { return this; }

    public int getLineNumber() { return line ; }
    public File getFile() { return fastqFile ; }

    public void close() {
        try {
            reader.close();
        } catch (IOException e) {
            throw new PicardException("IO problem in file "+fastqFile.getAbsolutePath(),e);
        }
    }

    private void checkLine(final String line, final String kind) {
        if (line == null) {
            throw new PicardException(error("File is too short - missing "+kind+" line"));
        }
        if (StringUtil.isBlank(line)) {
            throw new PicardException(error("Missing "+kind));
        }
    }

    private String error(final String msg) {
        return msg + " at line "+line+" in "+fastqFile.getAbsolutePath();
    }
}
