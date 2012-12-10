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

import net.sf.samtools.util.BlockCompressedOutputStream;
import net.sf.samtools.util.IOUtil;
import net.sf.samtools.util.Md5CalculatingOutputStream;
import net.sf.samtools.util.RuntimeIOException;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

/**
 * Create a SAMFileWriter for writing SAM or BAM.
 */
public class SAMFileWriterFactory {
    private static boolean defaultCreateIndexWhileWriting = Defaults.CREATE_INDEX;
    private boolean createIndex = defaultCreateIndexWhileWriting ;
    private static boolean defaultCreateMd5File = Defaults.CREATE_MD5;
    private boolean createMd5File = defaultCreateMd5File;
    private boolean useAsyncIo = Defaults.USE_ASYNC_IO;
    private int asyncOutputBufferSize = AsyncSAMFileWriter.DEFAULT_QUEUE_SIZE;
    private File tmpDir;


    private Integer maxRecordsInRam;

    /** Sets the default for whether to create md5Files for BAM files this factory. */
    public static void setDefaultCreateMd5File(final boolean createMd5File) {
        defaultCreateMd5File = createMd5File;
    }

    /** Sets whether to create md5Files for BAMs from this factory. */
    public SAMFileWriterFactory setCreateMd5File(final boolean createMd5File) {
        this.createMd5File = createMd5File;
        return this;
    }

    /**
     * Sets the default for subsequent SAMFileWriterFactories
     * that do not specify whether to create an index.
     * If a BAM (not SAM) file is created, the setting is true, and the file header specifies coordinate order,
     * then a BAM index file will be written along with the BAM file.
     *
     * @param setting whether to attempt to create a BAM index while creating the BAM file
     */
    public static void setDefaultCreateIndexWhileWriting(final boolean setting) {
        defaultCreateIndexWhileWriting = setting;
    }

    /**
     * Convenience method allowing newSAMFileWriterFactory().setCreateIndex(true);
     * Equivalent to SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(true); newSAMFileWriterFactory();
     * If a BAM (not SAM) file is created, the setting is true, and the file header specifies coordinate order,
     * then a BAM index file will be written along with the BAM file.
     *
     * @param setting whether to attempt to create a BAM index while creating the BAM file.
     * @return this factory object
     */
    public SAMFileWriterFactory setCreateIndex(final boolean setting){
        this.createIndex = setting;
        return this;
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
    public SAMFileWriterFactory setMaxRecordsInRam(int maxRecordsInRam) {
        this.maxRecordsInRam = maxRecordsInRam;
        return this;
    }

    /**
     * Turn on or off the use of asynchronous IO for writing output SAM and BAM files.  If true then
     * each SAMFileWriter creates a dedicated thread which is used for compression and IO activities.
     */
    public void setUseAsyncIo(final boolean useAsyncIo) {
        this.useAsyncIo = useAsyncIo;
    }

    /**
     * If and only if using asynchronous IO then sets the maximum number of records that can be buffered per
     * SAMFileWriter before producers will block when trying to write another SAMRecord.
     */
    public void setAsyncOutputBufferSize(final int asyncOutputBufferSize) {
        this.asyncOutputBufferSize = asyncOutputBufferSize;
    }
    
    /**
     * Set the temporary directory to use when sort data.
     * @param tmpDir Path to the temporary directory
     */
    public SAMFileWriterFactory setTempDirectory(final File tmpDir) {
        this.tmpDir = tmpDir;
        return this;
    }

    /**
     * Create a BAMFileWriter that is ready to receive SAMRecords.  Uses default compression level.
     * @param header entire header. Sort order is determined by the sortOrder property of this arg.
     * @param presorted if true, SAMRecords must be added to the SAMFileWriter in order that agrees with header.sortOrder.
     * @param outputFile where to write the output.
     */
    public SAMFileWriter makeBAMWriter(final SAMFileHeader header, final boolean presorted, final File outputFile) {
        return makeBAMWriter(header, presorted, outputFile, BlockCompressedOutputStream.getDefaultCompressionLevel());
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
        try {
            boolean createMd5File = this.createMd5File && IOUtil.isRegularPath(outputFile);
            if (this.createMd5File && !createMd5File) {
                System.err.println("Cannot create MD5 file for BAM because output file is not a regular file: " + outputFile.getAbsolutePath());
            }
            final BAMFileWriter ret = createMd5File
                    ? new BAMFileWriter(new Md5CalculatingOutputStream(new FileOutputStream(outputFile, false),
                        new File(outputFile.getAbsolutePath() + ".md5")), outputFile, compressionLevel)
                    : new BAMFileWriter(outputFile, compressionLevel);
            boolean createIndex = this.createIndex && IOUtil.isRegularPath(outputFile);
            if (this.createIndex && !createIndex) {
                System.err.println("Cannot create index for BAM because output file is not a regular file: " + outputFile.getAbsolutePath());
            }
            if (this.tmpDir!=null) ret.setTempDirectory(this.tmpDir);
            initializeBAMWriter(ret, header, presorted, createIndex);

            if (this.useAsyncIo) return new AsyncSAMFileWriter(ret, this.asyncOutputBufferSize);
            else return ret;
        }
        catch (IOException ioe) {
            throw new RuntimeIOException("Error opening file: " + outputFile.getAbsolutePath());
        }
    }

    private void initializeBAMWriter(final BAMFileWriter writer, final SAMFileHeader header, final boolean presorted, final boolean createIndex) {
        writer.setSortOrder(header.getSortOrder(), presorted);
        if (maxRecordsInRam != null) {
            writer.setMaxRecordsInRam(maxRecordsInRam);
        }
        writer.setHeader(header);
        if (createIndex && writer.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)){
            writer.enableBamIndexConstruction();
        }
    }

    /**
     * Create a SAMTextWriter that is ready to receive SAMRecords.
     * @param header entire header. Sort order is determined by the sortOrder property of this arg.
     * @param presorted if true, SAMRecords must be added to the SAMFileWriter in order that agrees with header.sortOrder.
     * @param outputFile where to write the output.
     */
    public SAMFileWriter makeSAMWriter(final SAMFileHeader header, final boolean presorted, final File outputFile) {
        try {
            final SAMTextWriter ret = this.createMd5File
                ? new SAMTextWriter(new Md5CalculatingOutputStream(new FileOutputStream(outputFile, false),
                    new File(outputFile.getAbsolutePath() + ".md5")))
                : new SAMTextWriter(outputFile);
            ret.setSortOrder(header.getSortOrder(), presorted);
            if (maxRecordsInRam != null) {
                ret.setMaxRecordsInRam(maxRecordsInRam);
            }
            ret.setHeader(header);

            if (this.useAsyncIo) return new AsyncSAMFileWriter(ret, this.asyncOutputBufferSize);
            else return ret;
        }
        catch (IOException ioe) {
            throw new RuntimeIOException("Error opening file: " + outputFile.getAbsolutePath());
        }
    }

    /**
     * Create a SAMTextWriter for writing to a stream that is ready to receive SAMRecords.
     * This method does not support the creation of an MD5 file
     * 
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

        if (this.useAsyncIo) return new AsyncSAMFileWriter(ret, this.asyncOutputBufferSize);
        else return ret;
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

}
