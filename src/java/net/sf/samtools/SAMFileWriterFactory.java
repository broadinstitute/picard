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

    private static boolean DefaultCreateIndexWhileWriting = false;
    private boolean createIndex = false;

    private Integer maxRecordsInRam;

    /**
     * Construct a factory that uses the default to decide whether or not to create bam index files while creating the bam file.
     */
    public SAMFileWriterFactory() {
        this.createIndex = DefaultCreateIndexWhileWriting;
    }

    /**
     * Sets the default for subsequent SAMFileWriterFactories
     * that do not specify whether to create an index
     *
     * @param setting whether to create the bam index on the fly while creating the bam file
     */
    public static void setDefaultCreateIndexWhileWriting(final boolean setting) {
        DefaultCreateIndexWhileWriting = setting;
    }

    /**
     * Convenience method allowing newSAMFileWriterFactory().setCreateIndex(true);
     * Equivalent to SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(true); newSAMFileWriterFactory();
     * @param setting whether to attempt to create a BAM index while creating the BAM file
     * @return this factory object
     */
    public SAMFileWriterFactory setCreateIndex(final boolean setting){
        this.createIndex = setting;
        return this;
    }

    /**
     * Create a BAMFileWriter that is ready to receive SAMRecords.  Uses default compression level.
     * @param header entire header. Sort order is determined by the sortOrder property of this arg.
     * @param presorted if true, SAMRecords must be added to the SAMFileWriter in order that agrees with header.sortOrder.
     * @param outputFile where to write the output.
     */
    public SAMFileWriter makeBAMWriter(final SAMFileHeader header, final boolean presorted, final File outputFile) {
        final BAMFileWriter ret = new BAMFileWriter(outputFile);
        initializeBAMWriter(ret, header, presorted);
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
        initializeBAMWriter(ret, header, presorted);
        return ret;
    }

    private void initializeBAMWriter(BAMFileWriter writer, SAMFileHeader header, boolean presorted) {
        writer.setSortOrder(header.getSortOrder(), presorted);
        if (createIndex && writer.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)){
            writer.enableBamIndexConstruction();
        }
        if (maxRecordsInRam != null) {
            writer.setMaxRecordsInRam(maxRecordsInRam);
        }
        writer.setHeader(header);
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
        if (maxRecordsInRam != null) {
            ret.setMaxRecordsInRam(maxRecordsInRam);
        }
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
        if (maxRecordsInRam != null) {
            ret.setMaxRecordsInRam(maxRecordsInRam);
        }
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
        return makeBAMWriter(header, presorted, outputFile);
    }

    /**
     * Before creating a writer that is not presorted, this method may be called in order to override
     * the default number of SAMRecords stored in RAM before spilling to disk
     * (c.f. SAMFileWriterImpl.MAX_RECORDS_IN_RAM).  When writing very large sorted SAM files, you may need
     * call this method in order to avoid running out of file handles.  The RAM available to the JVM may need
     * to be increased in order to hold the specified number of records in RAM.  This value affects the number
     * of records stored in subsequent calls to one of the make...() methods.
     *
     * @param maxRecordsInRam Number of records to store in RAM before spilling to temporary file when
     * creating a sorted SAM or BAM file.
     */
    public void setMaxRecordsInRam(int maxRecordsInRam) {
        this.maxRecordsInRam = maxRecordsInRam;
    }
}
