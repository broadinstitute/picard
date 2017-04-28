package picard.illumina.parser.readers;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.RuntimeIOException;
import picard.illumina.parser.BclData;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/**
 * ------------------------------------- CBCL Header -----------------------------------
 * Bytes 0 - 1     Version number, current version is 1    unsigned 16 bits little endian integer
 * Bytes 2 - 5     Header size                             unsigned 32 bits little endian integer
 * Byte 6          Number of bits per basecall             unsigned
 * Byte 7          Number of bits per q-score              unsigned
 * <p>
 * q-val mapping info
 * Bytes 0-3   Number of bins (B), zero indicates no mapping
 * B pairs of 4 byte values (if B > 0) {from, to}, {from, to}, {from, to} …from: quality score bin to: quality score
 * <p>
 * Number of tile records                                  unsigned 32bits little endian integer
 * <p>
 * gzip virtual file offsets, one record per tile
 * Bytes 0-3:  tile number
 * Bytes 4-7   Number of clusters that were written into the current block (required due to bit-packed q-scores)
 * unsigned 32 bit integer
 * <p>
 * Bytes 8-11      Uncompressed block size of the tile data (useful for sanity check whenexcluding non-PF clusters)
 * unsigned 32 bit integer
 * <p>
 * Bytes 12-15     Compressed block size of the tile data  unsigned 32 bit integer
 * <p>
 * non-PF clusters excluded flag 1: non-PF clusters are excluded 0: non-PF clusters are included
 * <p>
 * ------------------------------------- CBCL File Content -----------------------------------
 * <p>
 * N blocks of gzip files, where N is the number of tiles.
 * <p>
 * Each block consists of C number of basecall, quality score pairs where C is the number of clusters for the given tile.
 * <p>
 * Each basecall, quality score pair has the following format (assuming 2 bits are used for the basecalls):
 * Bits 0-1: Basecalls (respectively [A, C, G, T] for [00, 01, 10, 11])
 * Bits 2 and up: Quality score (unsigned Q bit little endian integer where Q is the number of bits per q-score).
 * For a two bit quality score, this is two clusters per byte where the bottom 4 bits are the first cluster and the
 * higher 4 bits are the second cluster.
 **/

public class CbclReader extends BaseBclReader implements CloseableIterator<BclData> {

    private static Log log = Log.getInstance(CbclReader.class);

    private ByteBuffer[] cachedTile;

    private int[] currentTile;

    private List<BclData> queue = new ArrayList<>();

    private final CycleData[] cycleData;
    private Map<Integer, File> filters;
    private Map<Integer, boolean[]> cachedFilter = new HashMap<>();

    void addFilters(Map<Integer, File> filters) {
        this.filters = filters;
    }

    private static final int INITIAL_HEADER_SIZE = 6;

    public CbclReader(final List<File> cbcls, final int[] outputLengths) {
        super(outputLengths);
        cycleData = new CycleData[cycles];
        cachedTile = new ByteBuffer[cycles];
        currentTile = new int[cycles];

        try {
            final ByteBuffer byteBuffer = ByteBuffer.allocate(INITIAL_HEADER_SIZE);
            byteBuffer.order(ByteOrder.LITTLE_ENDIAN);
            for (int i = 0; i < cycles; i++) {
                currentTile[i] = 0;
                final File bclFile = cbcls.get(i);

                final InputStream stream = open(bclFile, false, false, false);
                int read = stream.read(byteBuffer.array());

                //we need to read the first 6 bytes to determine the header size
                if (read != INITIAL_HEADER_SIZE) {
                    close();
                    throw new RuntimeIOException(String.format("BCL %s has invalid header structure.", bclFile.getAbsoluteFile()));
                }

                short version = byteBuffer.getShort();
                int headerSize = byteBuffer.getInt();

                final ByteBuffer headerBuffer = ByteBuffer.allocate(headerSize - INITIAL_HEADER_SIZE);
                headerBuffer.order(ByteOrder.LITTLE_ENDIAN);

                read = stream.read(headerBuffer.array());
                if (read != headerSize - INITIAL_HEADER_SIZE) {
                    close();
                    throw new RuntimeIOException(String.format("BCL %s has invalid header structure.", bclFile.getAbsoluteFile()));
                }

                byte bitsPerBasecall = headerBuffer.get();
                byte bitsPerQualityScore = headerBuffer.get();

                if (bitsPerBasecall != 2 && bitsPerBasecall != bitsPerQualityScore) {
                    close();
                    throw new RuntimeIOException("CBCL data not encoded in nibbles. (not currently supported)");
                }

                int numberOfBins = headerBuffer.getInt();

                byte[] qualityBins = new byte[numberOfBins];
                //each bin has a pair of 4 byte mappings
                for (int j = 0; j < numberOfBins; j++) {
                    headerBuffer.getInt(); // first int is "from" value, which we don't need
                    int to = headerBuffer.getInt();
                    qualityBins[j] = (byte) to;
                }

                int numTiles = headerBuffer.getInt();
                TileData[] tileInfo = new TileData[numTiles];
                for (int j = 0; j < numTiles; j++) {
                    int tileNum = headerBuffer.getInt();
                    int numClustersInTile = headerBuffer.getInt();
                    int uncompressedBlockSize = headerBuffer.getInt();
                    int compressedBlockSize = headerBuffer.getInt();
                    tileInfo[j] = new TileData(tileNum, numClustersInTile, uncompressedBlockSize, compressedBlockSize);
                }

                boolean pfExcluded = headerBuffer.get() == 1;
                cycleData[i] = new CycleData(version, headerSize, bitsPerBasecall, bitsPerQualityScore, numberOfBins, qualityBins, numTiles, tileInfo, pfExcluded);
                this.streams[i] = stream;
                this.streamFiles[i] = bclFile;
                byteBuffer.clear();
                headerBuffer.clear();
            }


        } catch (final IOException ioe) {
            throw new RuntimeIOException(ioe);
        }
    }

    @Override
    public boolean hasNext() {
        if (queue.isEmpty()) {
            advance();
        }
        return !queue.isEmpty();
    }

    public BclData next() {
        if (queue.isEmpty()) {
            advance();
        }

        final BclData data = queue.get(0);
        queue.remove(0);
        return data;
    }

    @Override
    public void close() {
        for (final InputStream stream : this.streams) {
            CloserUtil.close(stream);
        }
    }

    private void advance() {
        int totalCycleCount = 0;
        BclData data = new BclData(outputLengths);

        for (int read = 0; read < outputLengths.length; read++) {
            for (int cycle = 0; cycle < outputLengths[read]; cycle++) {
                try {
                    CycleData currentCycleData = cycleData[totalCycleCount];
                    TileData currentTileData = currentCycleData.tileInfo[currentTile[totalCycleCount]];
                    try {
                        if (cachedTile[totalCycleCount] == null) {
                            if (!cachedFilter.containsKey(currentTileData.tileNum)) {
                                cacheFilter(currentTileData);
                            }
                            cacheTile(totalCycleCount, currentTileData, currentCycleData);
                        }
                    } catch (IOException e) {
                        // when logging the error, increment cycle by 1, since totalCycleCount is zero-indexed but Illumina directories are 1-indexed.
                        throw new IOException(String.format("Error while reading from BCL file for cycle %d. Offending file on disk is %s",
                                (totalCycleCount + 1), this.streamFiles[totalCycleCount].getAbsolutePath()), e);
                    }

                    if (!cachedTile[totalCycleCount].hasRemaining()) {
                        //on to the next tile
                        currentTile[totalCycleCount]++;
                        //if past the last tile then return
                        if (currentTile[totalCycleCount] > currentCycleData.tileInfo.length - 1) {
                            return;
                        }
                        currentTileData = currentCycleData.tileInfo[currentTile[totalCycleCount]];
                        if (!cachedFilter.containsKey(currentTileData.tileNum)) {
                            cacheFilter(currentTileData);
                        }
                        cacheTile(totalCycleCount, currentTileData, currentCycleData);
                    }

                    int singleByte = cachedTile[totalCycleCount].get();

                    decodeQualityBinnedBasecall(data, read, cycle, singleByte, currentCycleData);

                    totalCycleCount++;
                } catch (final IOException ioe) {
                    throw new RuntimeIOException(ioe);
                }

            }
        }
        this.queue.add(data);
    }

    private void cacheFilter(TileData currentTileData) {
        boolean[] filterValues = new boolean[currentTileData.numClustersInTile];
        FilterFileReader reader = new FilterFileReader(filters.get(currentTileData.tileNum));
        int count = 0;
        while (reader.hasNext()) {
            filterValues[count] = reader.next();
            count++;
        }
        cachedFilter.put(currentTileData.tileNum, filterValues);
    }

    private void cacheTile(int totalCycleCount, TileData tileData, CycleData currentCycleData) throws IOException {
        ByteBuffer tileByteBuffer = ByteBuffer.allocate(tileData.compressedBlockSize);
        //we are going to explode the nibbles in to bytes to make PF filtering easier
        ByteBuffer tempUncompressedByteBuffer = ByteBuffer.allocate(tileData.uncompressedBlockSize);
        ByteBuffer uncompressedByteBuffer = ByteBuffer.allocate(tileData.uncompressedBlockSize);
        ByteBuffer unNibbledByteBuffer = ByteBuffer.allocate(tileData.uncompressedBlockSize * 2);

        tileByteBuffer.order(ByteOrder.LITTLE_ENDIAN);
        uncompressedByteBuffer.order(ByteOrder.LITTLE_ENDIAN);

        // Read the whole compressed block into a buffer, then sanity check the length
        int readBytes = this.streams[totalCycleCount].read(tileByteBuffer.array());
        if (readBytes != tileData.compressedBlockSize) {
            throw new IOException(String.format("Error while reading from BCL file for cycle %d. Offending file on disk is %s",
                    (totalCycleCount + 1), this.streamFiles[totalCycleCount].getAbsolutePath()));
        }

        // Uncompress the data from the buffer we just wrote - use gzip input stream to write to uncompressed buffer
        ByteArrayInputStream byteInputStream = new ByteArrayInputStream(tileByteBuffer.array());
        if (cachedTile[totalCycleCount] != null) {
            //clear the old cache
            cachedTile[totalCycleCount] = null;
        }
        GZIPInputStream gzipInputStream = new GZIPInputStream(byteInputStream, uncompressedByteBuffer.capacity());
        int read = 0;
        int totalRead = 0;
        while ((read = gzipInputStream.read(tempUncompressedByteBuffer.array(), 0, tempUncompressedByteBuffer.capacity())) != -1) {
            uncompressedByteBuffer.put(tempUncompressedByteBuffer.array(), 0, read);
            totalRead += read;
        }
        if (totalRead != tileData.uncompressedBlockSize) {
            throw new IOException(String.format("Error while decompressing from BCL file for cycle %d. Offending file on disk is %s",
                    (totalCycleCount + 1), this.streamFiles[totalCycleCount].getAbsolutePath()));
        }

        // Read uncompressed data from the buffer and expand each nibble into a full byte for ease of use
        uncompressedByteBuffer.flip();
        while (uncompressedByteBuffer.hasRemaining()) {
            byte singleByte = uncompressedByteBuffer.get();
            unNibbledByteBuffer.put((byte) (singleByte & 0x0f));
            unNibbledByteBuffer.put((byte) ((singleByte >> 4) & 0x0f));
        }
        gzipInputStream.close();

        // Write buffer contents to cached tile array
        unNibbledByteBuffer.flip();
        //if nonPF reads are included we need to strip them out
        if (!currentCycleData.pfExcluded) {
            ByteBuffer uncompressedFilteredByteBuffer = ByteBuffer.allocate(tileData.uncompressedBlockSize * 2);
            FilterFileReader reader = new FilterFileReader(filters.get(tileData.tileNum));
            while (reader.hasNext()) {
                byte readByte = unNibbledByteBuffer.get();
                if (reader.next()) {
                    uncompressedFilteredByteBuffer.put(readByte);
                }
            }
            uncompressedFilteredByteBuffer.flip();
            cachedTile[totalCycleCount] = uncompressedFilteredByteBuffer.slice();
        } else {
            cachedTile[totalCycleCount] = unNibbledByteBuffer;
        }
    }

}
