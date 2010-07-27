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

import net.sf.samtools.util.BinaryCodec;
import net.sf.samtools.util.BlockCompressedOutputStream;

import java.io.DataOutputStream;
import java.io.File;

/**
 * Concrete implementation of SAMFileWriter for writing gzipped BAM files.
 */
public class BAMFileWriter extends SAMFileWriterImpl {

    private final BinaryCodec outputBinaryCodec;
    private BAMRecordCodec bamRecordCodec = null;
    private final BlockCompressedOutputStream blockCompressedOutputStream;
    private BAMIndexer bamIndexer = null;

    private boolean makeBamIndex = false;

    public BAMFileWriter(final File path) {
        blockCompressedOutputStream = new BlockCompressedOutputStream(path);
        outputBinaryCodec = new BinaryCodec(new DataOutputStream(blockCompressedOutputStream));
        outputBinaryCodec.setOutputFileName(path.toString());
    }

    public BAMFileWriter(final File path, final int compressionLevel) {
        blockCompressedOutputStream = new BlockCompressedOutputStream(path, compressionLevel);
        outputBinaryCodec = new BinaryCodec(new DataOutputStream(blockCompressedOutputStream));
        outputBinaryCodec.setOutputFileName(path.toString());
    }

    // Allow enabling or disabling the bam index construction
    public void enableBamIndexConstruction (boolean setting) {
        if (setting && !getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)){
           throw new SAMException("Not creating BAM index since not sorted by coordinates: " + getSortOrder());
        }
        makeBamIndex = setting;
    }

    private void prepareToWriteAlignments() {
        if (bamRecordCodec == null) {
            bamRecordCodec = new BAMRecordCodec(getFileHeader());
            bamRecordCodec.setOutputStream(outputBinaryCodec.getOutputStream());
            if (makeBamIndex && bamIndexer == null) {
                bamIndexer = createBamIndex(outputBinaryCodec.getOutputFileName());
            }
        }
    }

    private BAMIndexer createBamIndex(String path) {
        try {
            String indexFileName = path + BAMIndex.BAMIndexSuffix;
            File indexFile = new File(indexFileName);
            if (indexFile.exists()) {
                if (!indexFile.canWrite()) {
                    throw new SAMException("Not creating BAM index since unable to write index file " + indexFileName);
                }
                indexFile.delete();
            }
            return new BAMIndexer(indexFile, getFileHeader().getSequenceDictionary().size(), false);
        } catch (Exception e) {
            throw new SAMException("Not creating BAM index", e);
        }
    }

    protected void writeAlignment(final SAMRecord alignment) {
        prepareToWriteAlignments();

        if (bamIndexer != null && blockCompressedOutputStream != null) {
            try {
                final long startCoordinate = blockCompressedOutputStream.getFilePointer();
                bamRecordCodec.encode(alignment);
                final long stopCoordinate = blockCompressedOutputStream.getFilePointer();
                // set the alignment's SourceInfo and then prepare its index information
                alignment.setFileSource(new SAMFileSource(null, new BAMFileSpan(new Chunk(startCoordinate, stopCoordinate))));
                bamIndexer.processAlignment(alignment);
            } catch (Exception e) {
                bamIndexer.deleteIndex();
                bamIndexer = null;
                throw new SAMException("Exception when processing alignment for BAM index " + alignment, e);
            }
        } else {
            bamRecordCodec.encode(alignment);
        }
    }

    protected void writeHeader(final String textHeader) {
        outputBinaryCodec.writeBytes(BAMFileConstants.BAM_MAGIC);

        // calculate and write the length of the SAM file header text and the header text
        outputBinaryCodec.writeString(textHeader, true, false);

        // write the sequences binarily.  This is redundant with the text header
        outputBinaryCodec.writeInt(getFileHeader().getSequenceDictionary().size());
        for (final SAMSequenceRecord sequenceRecord: getFileHeader().getSequenceDictionary().getSequences()) {
            outputBinaryCodec.writeString(sequenceRecord.getSequenceName(), true, true);
            outputBinaryCodec.writeInt(sequenceRecord.getSequenceLength());
        }
    }

    protected void finish() {
        outputBinaryCodec.close();
        if (bamIndexer != null){
            try {
                bamIndexer.finish();
            } catch (Exception e) {
                bamIndexer.deleteIndex(); // don't want a partial bai file around
                throw new SAMException("Exception writing BAM index file", e);
            }
        }
    }

    protected String getFilename() {
        return outputBinaryCodec.getOutputFileName();
    }
}
