package net.sf.picard.illumina.parser.fakers;

import java.nio.ByteBuffer;

/**
 * Created by jcarey on 3/13/14.
 */
public class QSeqFileFaker extends FileFaker {
    private final String qseqString = "Can not make qseq file";

    @Override
    protected void fakeFile(final ByteBuffer buffer) {
        buffer.put(qseqString.getBytes());
    }

    @Override
    protected boolean addLeadingZeros() {
        return false;
    }

    @Override
    protected int bufferSize() {
        return qseqString.getBytes().length;
    }
}