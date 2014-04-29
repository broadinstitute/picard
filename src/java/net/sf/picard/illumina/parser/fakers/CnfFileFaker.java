package net.sf.picard.illumina.parser.fakers;

import java.nio.ByteBuffer;

/**
 * Created by jcarey on 3/13/14.
 */
public class CnfFileFaker extends FileFaker {

    @Override
    protected void fakeFile(final ByteBuffer buffer) {
        buffer.putChar('C');
        buffer.putChar('I');
        buffer.putChar('F');
        buffer.put((byte) 1);
        buffer.put((byte) 1);
        buffer.putShort((short) 1);
        buffer.putShort((short) 1);
        buffer.putInt(1);
    }

    @Override
    protected boolean addLeadingZeros() {
        return false;
    }

    @Override
    protected int bufferSize() {
        return Integer.SIZE + (Character.SIZE * 3) + (Short.SIZE * 2) + 2;
    }
}