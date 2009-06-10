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
package net.sf.samtools;

import net.sf.samtools.util.AsciiWriter;
import net.sf.samtools.util.RuntimeIOException;

import java.io.*;

/**
 * Writer for text-format SAM files.
 */
class SAMTextWriter extends SAMFileWriterImpl {
    private static final String FIELD_SEPARATOR = "\t";

    private final Writer out;
    // For error reporting only.
    private final File file;
    private final TextTagCodec tagCodec = new TextTagCodec();
    private final SAMTagUtil tagUtil = new SAMTagUtil();

    /**
     * Prepare to write SAM text file.
     * @param file Where to write the output.
     */
    SAMTextWriter(final File file) {
        try {
            this.file = file;
            this.out = new AsciiWriter(new FileOutputStream(file));
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    /**
     * Constructs a SAMTextWriter for outputting to a stream instead of to a file.
     * @param stream Need not be buffered because this class provides buffering. 
     */
    SAMTextWriter(final OutputStream stream) {
        this.file = null;
        this.out = new AsciiWriter(stream);
    }

    /**
     * Writes the record to disk.  Sort order has been taken care of by the time
     * this method is called.
     *
     * @param alignment
     */
    protected void writeAlignment(final SAMRecord alignment) {
        try {
            out.write(alignment.getReadName());
            out.write(FIELD_SEPARATOR);
            out.write(Integer.toString(alignment.getFlags()));
            out.write(FIELD_SEPARATOR);
            out.write(alignment.getReferenceName());
            out.write(FIELD_SEPARATOR);
            out.write(Integer.toString(alignment.getAlignmentStart()));
            out.write(FIELD_SEPARATOR);
            out.write(Integer.toString(alignment.getMappingQuality()));
            out.write(FIELD_SEPARATOR);
            out.write(alignment.getCigarString());
            out.write(FIELD_SEPARATOR);

            //  == is OK here because these strings are interned
            if (alignment.getReferenceName() == alignment.getMateReferenceName() &&
                    SAMRecord.NO_ALIGNMENT_REFERENCE_NAME != alignment.getReferenceName()) {
                out.write("=");
            } else {
                out.write(alignment.getMateReferenceName());
            }
            out.write(FIELD_SEPARATOR);
            out.write(Integer.toString(alignment.getMateAlignmentStart()));
            out.write(FIELD_SEPARATOR);
            out.write(Integer.toString(alignment.getInferredInsertSize()));
            out.write(FIELD_SEPARATOR);
            out.write(alignment.getReadString());
            out.write(FIELD_SEPARATOR);
            out.write(alignment.getBaseQualityString());
            if (alignment.getBinaryAttributes() != null) {
                for (final SAMBinaryTagAndValue attribute : alignment.getBinaryAttributes()) {
                    out.write(FIELD_SEPARATOR);
                    out.write(tagCodec.encode(tagUtil.makeStringTag(attribute.tag), attribute.value));
                }
            }
            out.write("\n");

        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    /**
     * Write the header to disk.  Header object is available via getHeader().
     *
     * @param textHeader for convenience if the implementation needs it.  Must be newline-terminated.
     */
    protected void writeHeader(final String textHeader) {
        try {
            out.write(textHeader);
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    /**
     * Do any required flushing here.
     */
    protected void finish() {
        try {
            out.close();
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    /**
     * For producing error messages.
     *
     * @return Output filename, or null if there isn't one.
     */
    protected String getFilename() {
        if (file == null) {
            return null;
        }
        return file.getAbsolutePath();
    }
}
