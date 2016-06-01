package picard.illumina.parser.fakers;

import java.nio.ByteBuffer;

/**
 * Created by jcarey on 3/13/14.
 */
public class BarcodeFileFaker extends FileFaker {
    private static final String BARCODE_STRING = "1\tn\t \n";

    @Override
    protected void fakeFile(final ByteBuffer buffer) {
        buffer.put(BARCODE_STRING.getBytes());
    }

    @Override
    protected boolean addLeadingZeros() {
        return false;
    }

    @Override
    protected int bufferSize() {
        return BARCODE_STRING.getBytes().length;
    }
}
