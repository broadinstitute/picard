package net.sf.picard.illumina.parser.fakers;

import java.nio.ByteBuffer;

public class ClocsFileFaker extends FileFaker {

    @Override
    protected void fakeFile(final ByteBuffer buffer) {
        buffer.put((byte) 1);
        buffer.putInt(1);
        buffer.put((byte) (0xff & 1));
        buffer.putFloat((byte) (0xff & 5));
        buffer.putFloat((byte) (0xff & 5));
    }

    @Override
    protected boolean addLeadingZeros() {
        return false;
    }

    @Override
    protected int bufferSize() {
        return (Integer.SIZE * 2) + (Float.SIZE * 3);
    }
}