package picard.illumina.parser.fakers;

import java.nio.ByteBuffer;

/**
 * Created by jcarey on 3/13/14.
 */
public class PosFileFaker extends FileFaker {
    private static final String POS_FILE_STRING = "102.0\t303.3\n";

    @Override
    protected void fakeFile(final ByteBuffer buffer) {
        buffer.put(POS_FILE_STRING.getBytes());
    }

    @Override
    protected boolean addLeadingZeros() {
        return false;
    }

    @Override
    protected int bufferSize() {
        return POS_FILE_STRING.getBytes().length;
    }
}
