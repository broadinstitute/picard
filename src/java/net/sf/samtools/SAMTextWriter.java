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
public class SAMTextWriter extends SAMFileWriterImpl {
    private static final String FIELD_SEPARATOR = "\t";

    private final Writer out;
    // For error reporting only.
    private final File file;
    private final TextTagCodec tagCodec = new TextTagCodec();
    private final SAMTagUtil tagUtil = new SAMTagUtil();

    /**
     * Constructs a SAMTextWriter that outputs to a Writer.
     * @param out Writer.
     */
    public SAMTextWriter(Writer out) {
	this.out = out;
	this.file = null;
    }

    /**
     * Constructs a SAMTextWriter that writes to a File.
     * @param file Where to write the output.
     */
    public SAMTextWriter(final File file) {
        try {
            this.file = file;
            this.out = new AsciiWriter(new FileOutputStream(file));
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    /**
     * Returns the Writer used by this instance.  Useful for flushing the output.
     */
    public Writer getWriter() {
	return out;
    }

    /**
     * Constructs a SAMTextWriter that writes to an OutputStream.  The OutputStream
     * is wrapped in an AsciiWriter, which can be retrieved with getWriter().
     * @param stream Need not be buffered because this class provides buffering. 
     */
    public SAMTextWriter(final OutputStream stream) {
        this.file = null;
        this.out = new AsciiWriter(stream);
    }

    /**
     * Write the record.
     *
     * @param alignment SAMRecord.
     */
    public void writeAlignment(final SAMRecord alignment) {
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
            SAMBinaryTagAndValue attribute = alignment.getBinaryAttributes();
            while (attribute != null) {
                out.write(FIELD_SEPARATOR);
                final String encodedTag;
                if (attribute.isUnsignedArray()) {
                    encodedTag = tagCodec.encodeUnsignedArray(tagUtil.makeStringTag(attribute.tag), attribute.value);
                } else {
                    encodedTag = tagCodec.encode(tagUtil.makeStringTag(attribute.tag), attribute.value);
                }
                out.write(encodedTag);
                attribute = attribute.getNext();
            }
            out.write("\n");

        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    /* This method is called by SAMRecord.getSAMString(). */
    private static SAMTextWriter textWriter = null;
    private static StringWriter stringWriter = null;
    static synchronized String getSAMString(final SAMRecord alignment) {
	if (stringWriter == null) stringWriter = new StringWriter();
	if (textWriter == null) textWriter = new SAMTextWriter(stringWriter);
	stringWriter.getBuffer().setLength(0);
	textWriter.writeAlignment(alignment);
	return stringWriter.toString();
    }

    /**
     * Write the header text.  This method can also be used to write
     * an arbitrary String, not necessarily the header.
     *
     * @param textHeader String containing the text to write.
     */
    public void writeHeader(final String textHeader) {
        try {
            out.write(textHeader);
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    /**
     * Do any required flushing here.
     */
    public void finish() {
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
    public String getFilename() {
        if (file == null) {
            return null;
        }
        return file.getAbsolutePath();
    }
}
