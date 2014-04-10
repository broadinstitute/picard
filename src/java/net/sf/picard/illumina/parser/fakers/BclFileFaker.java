package net.sf.picard.illumina.parser.fakers;

import java.nio.ByteBuffer;

/**
 * Created by jcarey on 3/13/14.
 */
public class BclFileFaker extends FileFaker {

    @Override
    public void fakeFile(final ByteBuffer buffer) {

        // Write the number of elements to the header. There's one element per byte in the file,
        // minus the size of the header. Thus we remove four bytes from the size and THEN write
        // its value.
        size -= 4;
        buffer.putInt(size);

        while (size > 0) {
            // Fill the file with no calls
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
