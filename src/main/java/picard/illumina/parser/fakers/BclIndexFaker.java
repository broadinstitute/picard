package picard.illumina.parser.fakers;

import htsjdk.samtools.util.CloserUtil;
import picard.illumina.parser.TileIndex;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;

public class BclIndexFaker extends FileFaker {

    @Override
    protected void fakeFile(final ByteBuffer buffer) {
        buffer.putInt(0);
        buffer.putInt(tiles.size());

        for(int i = 0; i < tiles.size(); i++) {
            buffer.putLong(1);
        }
    }

    @Override
    protected boolean addLeadingZeros() {
        return false;
    }

    @Override
    protected int bufferSize() {
        return 8 * tiles.size();
    }

    public void fakeBciFile(final File bci, final TileIndex tileIndex) throws IOException {
        tiles = tileIndex.getTiles();
        final FileOutputStream fileOutputStream = new FileOutputStream(bci);
        final FileChannel channel = fileOutputStream.getChannel();
        final ByteBuffer buffer = ByteBuffer.allocate(8 + (8 * tiles.size()));
        buffer.order(ByteOrder.LITTLE_ENDIAN);

        fakeFile(buffer);

        buffer.flip();

        channel.write(buffer);
        channel.force(true);

        CloserUtil.close(channel);
        CloserUtil.close(fileOutputStream);
    }
}
