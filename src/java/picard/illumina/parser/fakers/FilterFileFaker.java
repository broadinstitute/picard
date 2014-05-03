package picard.illumina.parser.fakers;

import java.nio.ByteBuffer;

/**
 * Created by jcarey on 3/13/14.
 */
public class FilterFileFaker extends FileFaker {

    @Override
    protected void fakeFile(final ByteBuffer buffer) {
        buffer.putInt(0);
        buffer.putInt(3);
        buffer.putInt(1);
    }

    @Override
    protected boolean addLeadingZeros() {
        return true;
    }

    @Override
    protected int bufferSize() {
        return Integer.SIZE * 3;
    }
}
