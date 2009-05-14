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

import java.io.File;
import java.io.OutputStream;

/**
 * Create a SAMFileWriter for writing SAM or BAM.
 */
public class SAMFileWriterFactory {

    /**
     * Create a BAMFileWriter that is ready to receive SAMRecords.  Uses default compression level.
     * @param header entire header. Sort order is determined by the sortOrder property of this arg.
     * @param presorted if true, SAMRecords must be added to the SAMFileWriter in order that agrees with header.sortOrder.
     * @param outputFile where to write the output.
     */
    public SAMFileWriter makeBAMWriter(final SAMFileHeader header, final boolean presorted, final File outputFile) {
        final BAMFileWriter ret = new BAMFileWriter(outputFile);
        ret.setSortOrder(header.getSortOrder(), presorted);
        ret.setHeader(header);
        return ret;
    }

    /**
     *
     * Create a BAMFileWriter that is ready to receive SAMRecords.
     * @param header entire header. Sort order is determined by the sortOrder property of this arg.
     * @param presorted if true, SAMRecords must be added to the SAMFileWriter in order that agrees with header.sortOrder.
     * @param outputFile where to write the output.
     * @param compressionLevel Override default compression level with the given value, between 0 (fastest) and 9 (smallest).
     */
    public SAMFileWriter makeBAMWriter(final SAMFileHeader header, final boolean presorted, final File outputFile,
                                       final int compressionLevel) {
        final BAMFileWriter ret = new BAMFileWriter(outputFile, compressionLevel);
        ret.setSortOrder(header.getSortOrder(), presorted);
        ret.setHeader(header);
        return ret;
    }

    /**
     * Create a SAMTextWriter that is ready to receive SAMRecords.
     * @param header entire header. Sort order is determined by the sortOrder property of this arg.
     * @param presorted if true, SAMRecords must be added to the SAMFileWriter in order that agrees with header.sortOrder.
     * @param outputFile where to write the output.
     */
    public SAMFileWriter makeSAMWriter(final SAMFileHeader header, final boolean presorted, final File outputFile) {
        final SAMTextWriter ret = new SAMTextWriter(outputFile);
        ret.setSortOrder(header.getSortOrder(), presorted);
        ret.setHeader(header);
        return ret;
    }

    /**
     * Create a SAMTextWriter for writing to a stream that is ready to receive SAMRecords.
     * @param header entire header. Sort order is determined by the sortOrder property of this arg.
     * @param presorted if true, SAMRecords must be added to the SAMFileWriter in order that agrees with header.sortOrder.
     * @param stream the stream to write records to.
     */
    public SAMFileWriter makeSAMWriter(final SAMFileHeader header, final boolean presorted, final OutputStream stream) {
        final SAMTextWriter ret = new SAMTextWriter(stream);
        ret.setSortOrder(header.getSortOrder(), presorted);
        ret.setHeader(header);
        return ret;
    }


    /**
     * Create either a SAM or a BAM writer based on examination of the outputFile extension.
     * @param header entire header. Sort order is determined by the sortOrder property of this arg.
     * @param presorted presorted if true, SAMRecords must be added to the SAMFileWriter in order that agrees with header.sortOrder.
     * @param outputFile where to write the output.  Must end with .sam or .bam.
     * @return SAM or BAM writer based on file extension of outputFile.
     */
    public SAMFileWriter makeSAMOrBAMWriter(final SAMFileHeader header, final boolean presorted, final File outputFile) {
        final String filename = outputFile.getName();
        if (filename.endsWith(".bam")) {
            return makeBAMWriter(header, presorted, outputFile);
        }
        if (filename.endsWith(".sam")) {
            return makeSAMWriter(header, presorted, outputFile);
        }
        throw new IllegalArgumentException("SAM/BAM file should end with .sam or .bam: " + outputFile);
    }
}
