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
import net.sf.picard.io.IoUtil;

import java.io.File;
import java.io.PrintStream;

/**
 * @author alecw@broadinstitute.org
 */
public class FastqWriter {
    private final File file;
    private final PrintStream writer;

    public FastqWriter(final File file) {
        this.file = file;
        this.writer = new PrintStream(IoUtil.openFileForWriting(file));
    }

    public void write(final FastqRecord rec) {
        writer.print(FastqConstants.SEQUENCE_HEADER);
        writer.println(rec.getReadHeader());
        writer.println(rec.getReadString());
        writer.print(FastqConstants.QUALITY_HEADER);
        writer.println(rec.getBaseQualityHeader() == null ? "" : rec.getBaseQualityHeader());
        writer.println(rec.getBaseQualityString());
        if (writer.checkError()) {
            throw new PicardException("Error in writing file "+file);
        }
    }

    public void close() {
        writer.close();
    }
}
