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
import java.io.OutputStream;
import java.io.PrintStream;

/**
 * In general FastqWriterFactory should be used so that AsyncFastqWriter can be enabled, but there are some
 * cases in which that behavior is explicitly not wanted.
 */
public class BasicFastqWriter implements FastqWriter {
    private final String path;
    private final PrintStream writer;

    public BasicFastqWriter(final File file) {
        this(file, false);
    }

    public BasicFastqWriter(final File file, final boolean createMd5) {
        this(file, new PrintStream(maybeMd5Wrap(file, createMd5)));
    }

    private BasicFastqWriter(final File file, final PrintStream writer) {
        this.path = (file != null? file.getAbsolutePath(): "");
        this.writer = writer;
    }

    public BasicFastqWriter(final PrintStream writer) {
        this(null, writer);
    }

    public void write(final FastqRecord rec) {
        writer.print(FastqConstants.SEQUENCE_HEADER);
        writer.println(rec.getReadHeader());
        writer.println(rec.getReadString());
        writer.print(FastqConstants.QUALITY_HEADER);
        writer.println(rec.getBaseQualityHeader() == null ? "" : rec.getBaseQualityHeader());
        writer.println(rec.getBaseQualityString());
        if (writer.checkError()) {
            throw new PicardException("Error in writing fastq file " + path);
        }
    }

    public void flush() {
        writer.flush();
    }

    public void close() {
        writer.close();
    }

    private static OutputStream maybeMd5Wrap(final File file, final boolean createMd5) {
        if (createMd5) {
            return IoUtil.openFileForMd5CalculatingWriting(file);
        } else {
            return IoUtil.openFileForWriting(file);
        }
    }
}
