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
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;

/**
 * Class for writing BAI (BAM file indexes) as text or binary
 */
public class BAMIndexTextWriter extends CachingBAMFileIndex {

    /**
     * The number of references (chromosomes) in the BAI file
     * The name n_ref corresponds to the name in the SAM Format Specification Document
     */
    public final int n_ref;

    // private final Log log = Log.getInstance(getClass());

    private final File OUTPUT;  // todo could use mFile in super, though currently private

    /**
     *
     * @param INPUT     A BAM Index File, .bai
     * @param OUTPUT    A Textual BAM Index File, .bai.txt, or a binary bai file .generated.bai
     */
    public BAMIndexTextWriter(final File INPUT, final File OUTPUT) {
        super(INPUT);
        this.OUTPUT = OUTPUT;
        if (OUTPUT.exists()) {
            OUTPUT.delete();
        }
        n_ref = getNumberOfReferences();
    }

    public void writeText(final boolean sortBins) throws Exception {

        final PrintWriter pw = new PrintWriter(OUTPUT);
        pw.println("n_ref=" + n_ref);
        for (int i = 0 ; i < n_ref; i++){
            getQueryResults(i).writeText(pw, sortBins);
        }
        pw.close();
    }

     public void writeBinary(final boolean sortBins, final long bamFileSize) throws Exception {

        final int bufferSize; //  = 1000000; // 1M  works, but doesn't need to be this big
        final int defaultBufferSize = 1000000;  // 1M
        if (bamFileSize < defaultBufferSize  && bamFileSize != 0) {
            bufferSize = (int) bamFileSize;
        } else {
            bufferSize = defaultBufferSize;
        }
        // log.info("ByteBuffer size is " + bufferSize);

        final FileOutputStream stream = new FileOutputStream(OUTPUT, true);
        final FileChannel fileChannel = stream.getChannel();
        final ByteBuffer bb = ByteBuffer.allocateDirect(bufferSize);
        bb.order(ByteOrder.LITTLE_ENDIAN);

        // magic string
        final byte[] magic = BAMFileConstants.BAM_INDEX_MAGIC;
        bb.put(magic);
        // n_ref
        bb.putInt(n_ref);
        for (int i = 0; i < n_ref; i++) {
            getQueryResults(i).writeBinary(bb, sortBins);
            //  write out data and reset the buffer for each reference
            bb.flip();
            fileChannel.write(bb, 0);
            // stream.flush();    // todo will flushing the stream at every reference help memory?
            bb.position(0);
            bb.limit(bufferSize);
        }
        bb.flip();
        fileChannel.write(bb, 0);
        fileChannel.close();
        stream.close();
    }
}