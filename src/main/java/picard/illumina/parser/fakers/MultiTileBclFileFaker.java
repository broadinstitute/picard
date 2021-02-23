package picard.illumina.parser.fakers;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.CloserUtil;
import picard.illumina.parser.TileIndex;
import picard.illumina.parser.readers.BclReader;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.zip.GZIPOutputStream;

/**
 * Created by jcarey on 3/13/14.
 */
public class MultiTileBclFileFaker extends FileFaker {

    @Override
    protected void fakeFile(final ByteBuffer buffer) {
        buffer.putInt(1);
        long perTileSize = size;
        while (perTileSize > 0) {
            //fill the file with no calls
            buffer.put((byte) 0);
            perTileSize--;
        }
    }

    @Override
    protected boolean addLeadingZeros() {
        return false;
    }

    @Override
    protected int bufferSize() {
        return ((size - Integer.SIZE) * tiles.size()) + Integer.SIZE;
    }

    public void fakeMultiTileBclFile(final File bcl, final TileIndex tileIndex) throws IOException {
        this.tiles = tileIndex.getTiles();
        final OutputStream outputStream;
        if (BclReader.isGzipped(bcl)) {
            outputStream = new GZIPOutputStream(new FileOutputStream(bcl));
        } else if (BclReader.isBlockGzipped(bcl)) {
            outputStream = new BlockCompressedOutputStream(bcl);
        } else {
            outputStream = new FileOutputStream(bcl);
        }

        for (TileIndex.TileIndexRecord tileRecord : tileIndex) {
            final ByteBuffer buffer = ByteBuffer.allocate(tileRecord.getNumClustersInTile() + 4);
            buffer.order(ByteOrder.LITTLE_ENDIAN);

            this.size = tileRecord.getNumClustersInTile();
            fakeFile(buffer);
            buffer.flip();

            outputStream.write(buffer.array());
        }

        CloserUtil.close(outputStream);
    }
}
