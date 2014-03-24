package net.sf.picard.illumina.parser.fakers;

import java.nio.ByteBuffer;

/**
 * Created by jcarey on 3/13/14.
 */
public class BclFileFaker extends FileFaker {

    @Override
    public void fakeFile(final ByteBuffer buffer) {
        buffer.putInt(1);
        size -= 4;
        while (size > 0) {
            //fill the file with no calls
            buffer.put((byte) 0);
            size--;
        }
    }

    @Override
    protected boolean addLeadingZeros() {
        return false;
    }

    protected int bufferSize() {
        return size;
    }
}
