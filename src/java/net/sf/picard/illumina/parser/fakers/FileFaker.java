package net.sf.picard.illumina.parser.fakers;

import net.sf.samtools.util.CloserUtil;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.Collections;
import java.util.List;

/**
 * Created by jcarey on 3/13/14.
 */
public abstract class FileFaker {

    int size;
    List<Integer> tiles;

    protected abstract void fakeFile(ByteBuffer buffer);

    protected abstract boolean addLeadingZeros();

    protected abstract int bufferSize();

    public void fakeFile(final File base, final int tile, final int lane, final String extension) throws IOException {
        fakeFile(base, Collections.singletonList(tile), lane, extension);
    }

    public void fakeFile(final File base, final List<Integer> expectedTiles, final int lane, final String extension)
            throws IOException {
        if (base.exists() || base.mkdirs()) {
            this.tiles = expectedTiles;
            final File fakeFile;
            if (expectedTiles.size() == 1) {
                String longTileName = String.valueOf(tiles.get(0));
                if (addLeadingZeros()) {
                    while (longTileName.length() < 4) {
                        longTileName = "0" + longTileName;
                    }
                }
                fakeFile = new File(base, String.format("s_%d_%s%s", lane, longTileName, extension));
            } else {
                fakeFile = new File(base, String.format("s_%s%s", lane, extension));
            }

            fakeFile(fakeFile, bufferSize());
        }

    }

    public void fakeFile(final File cycleFile, Integer size) throws IOException {
        if (size == null) {
            size = 1;
        }
        this.size = size;
        final FileOutputStream fileOutputStream = new FileOutputStream(cycleFile);
        final FileChannel channel = fileOutputStream.getChannel();
        final ByteBuffer buffer = ByteBuffer.allocate(size);
        buffer.order(ByteOrder.LITTLE_ENDIAN);

        fakeFile(buffer);

        buffer.flip();

        channel.write(buffer);
        channel.force(true);

        CloserUtil.close(channel);
        CloserUtil.close(fileOutputStream);
    }
}
