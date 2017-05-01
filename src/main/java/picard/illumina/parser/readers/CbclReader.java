package picard.illumina.parser.readers;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;
import picard.PicardException;
import picard.illumina.parser.CbclData;

import java.io.ByteArrayInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
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
 * B pairs of 4 byte values (if B > 0) {from, to}, {from, to}, {from, to} from: quality score bin to: quality score
 * <p>
 * Number of tile records                                  unsigned 32bits little endian integer
 * <p>
 * gzip virtual file offsets, one record per tile
 * Bytes 0-3:  tile number
 * Bytes 4-7   Number of clusters that were written into the current block (required due to bit-packed q-scores)
 * unsigned 32 bit integer
 * <p>
 * Bytes 8-11      Uncompressed block size of the tile data (useful for sanity check when excluding non-PF clusters)
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

public class CbclReader extends BaseBclReader implements CloseableIterator<CbclData> {

    private final byte[][] cachedTile;
    private final int[] cachedTilePosition;

    private CbclData queue = null;
    private Iterator<AbstractIlluminaPositionFileReader.PositionInfo> positionInfoIterator;

    private final CycleData[] cycleData;
    private final Map<Integer, File> filterFileMap;
    private final Map<Integer, boolean[]> cachedFilter = new HashMap<>();
    private final Map<Integer, Map<Integer, File>> surfaceToTileToCbclMap;

    private static final int INITIAL_HEADER_SIZE = 6;
    private static final Log log = Log.getInstance(CbclReader.class);

    public CbclReader(final List<File> cbcls, final Map<Integer, File> filterFileMap, final int[] outputLengths, int tileNum, List<AbstractIlluminaPositionFileReader.PositionInfo> locs) {
        super(outputLengths);
        surfaceToTileToCbclMap = sortCbcls(cbcls);
        this.filterFileMap = filterFileMap;
        cycleData = new CycleData[cycles];
        cachedTile = new byte[cycles][];
        cachedTilePosition = new int[cycles];
        readSurfaceTile(tileNum, locs);
        close();
    }

    private void readSurfaceTile(int tileNum, List<AbstractIlluminaPositionFileReader.PositionInfo> locs) {
        log.info("Processing tile " + tileNum);
        try {
            for (Map.Entry<Integer, Map<Integer, File>> entry : surfaceToTileToCbclMap.entrySet()) {
                Map<Integer, File> cycleMap = entry.getValue();
                for (int i = 0; i < cycles; i++) {
                    //cycleMap is 1 indexed
                    final ByteBuffer byteBuffer = ByteBuffer.allocate(INITIAL_HEADER_SIZE);
                    byteBuffer.order(ByteOrder.LITTLE_ENDIAN);

                    final File bclFile = cycleMap.get(i + 1);
                    if (bclFile == null) {
                        throw new PicardException("Expected cbcl file for surface " + entry.getKey() + " cycle " + (i + 1) + " but it was not found.");
                    }

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
                        throw new PicardException(String.format("BCL %s has invalid header structure.", bclFile.getAbsoluteFile()));
                    }

                    byte bitsPerBasecall = headerBuffer.get();
                    byte bitsPerQualityScore = headerBuffer.get();

                    if (bitsPerBasecall != 2 && bitsPerBasecall != bitsPerQualityScore) {
                        close();
                        throw new PicardException("CBCL data not encoded in nibbles. (not currently supported)");
                    }

                    int numberOfBins = headerBuffer.getInt();

                    byte[] qualityBins = new byte[numberOfBins];
                    //each bin has a pair of 4 byte mappings
                    for (int j = 0; j < numberOfBins; j++) {
                        headerBuffer.getInt(); // first int is "from" value, which we don't need
                        int to = headerBuffer.getInt();
                        qualityBins[j] = (byte) to;
                    }
                    long filePos = 0;

                    int numTiles = headerBuffer.getInt();
                    TileData tileInfo = null;
                    for (int j = 0; j < numTiles; j++) {
                        int tile = headerBuffer.getInt();
                        int numClustersInTile = headerBuffer.getInt();
                        int uncompressedBlockSize = headerBuffer.getInt();
                        int compressedBlockSize = headerBuffer.getInt();
                        if (tile == tileNum) {
                            tileInfo = new TileData(tile, numClustersInTile, uncompressedBlockSize, compressedBlockSize, filePos);
                        }
                        filePos += compressedBlockSize;
                    }

                    boolean pfExcluded = headerBuffer.get() == 1;
                    //try the next surface if we didn't find the tile
                    if (tileInfo == null) continue;

                    cycleData[i] = new CycleData(version, headerSize, bitsPerBasecall, bitsPerQualityScore, numberOfBins, qualityBins, numTiles, tileInfo, pfExcluded);
                    this.streams[i] = stream;
                    this.streamFiles[i] = bclFile;
                    byteBuffer.clear();
                    headerBuffer.clear();
                }
            }
            int totalCycleCount = 0;

            if (cycleData[totalCycleCount].tileInfo == null) {
                throw new PicardException("Could not find tile " + cycleData[totalCycleCount].tileInfo);
            }

            for (int outputLength : outputLengths) {
                for (int cycle = 0; cycle < outputLength; cycle++) {
                    CycleData currentCycleData = cycleData[totalCycleCount];
                    try {
                        if (cachedTile[totalCycleCount] == null) {
                            if (!cachedFilter.containsKey(cycleData[totalCycleCount].tileInfo.tileNum)) {
                                cacheFilterAndLocs(cycleData[totalCycleCount].tileInfo, locs);
                            }
                            cacheTile(totalCycleCount, cycleData[totalCycleCount].tileInfo, currentCycleData);
                        }
                    } catch (IOException e) {
                        // when logging the error, increment cycle by 1, since totalCycleCount is zero-indexed but Illumina directories are 1-indexed.
                        throw new PicardException(String.format("Error while reading from BCL file for cycle %d. Offending file on disk is %s",
                                (totalCycleCount + 1), this.streamFiles[totalCycleCount].getAbsolutePath()), e);
                    }
                    totalCycleCount++;
                }
            }

        } catch (final IOException ioe) {
            throw new RuntimeIOException(ioe);
        }
    }

    private Map<Integer, Map<Integer, File>> sortCbcls(List<File> cbcls) {
        Map<Integer, Map<Integer, File>> sortedMap = new TreeMap<>();
        Pattern pattern = Pattern.compile("^.+C(\\d{1,4}).+L(\\d{1,3})_(\\d).cbcl$");
        for (File cbcl : cbcls) {
            Matcher matcher = pattern.matcher(cbcl.getAbsolutePath());
            if (!matcher.matches()) {
                throw new PicardException("CBCL File " + cbcl.getAbsolutePath() + " does not match expected pattern.");
            }
            Integer surface = Integer.valueOf(matcher.group(3));
            Integer cycle = Integer.valueOf(matcher.group(1));
            if (sortedMap.containsKey(surface)) {
                sortedMap.get(surface).put(cycle, cbcl);
            } else {
                Map<Integer, File> cycleMap = new HashMap<>();
                cycleMap.put(cycle, cbcl);
                sortedMap.put(surface, cycleMap);
            }
        }

        return sortedMap;
    }

    @Override
    public boolean hasNext() {
        if (queue == null) {
            advance();
        }
        return !(queue == null);
    }

    public CbclData next() {
        if (queue == null) {
            advance();
        }
        CbclData data = queue;
        queue = null;
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
        CbclData data = new CbclData(outputLengths, cycleData[totalCycleCount].tileInfo.tileNum);

        for (int read = 0; read < outputLengths.length; read++) {
            for (int cycle = 0; cycle < outputLengths[read]; cycle++) {
                CycleData currentCycleData = cycleData[totalCycleCount];

                if (cachedTilePosition[totalCycleCount] >= cachedTile[totalCycleCount].length) {
                    // end of tile
                    return;
                }

                int singleByte = cachedTile[totalCycleCount][cachedTilePosition[totalCycleCount]++];

                decodeQualityBinnedBasecall(data, read, cycle, singleByte, currentCycleData);

                totalCycleCount++;
            }
        }
        data.setPositionInfo(positionInfoIterator.next());
        this.queue = data;
    }

    private void cacheFilterAndLocs(TileData currentTileData, List<AbstractIlluminaPositionFileReader.PositionInfo> locs) {
        boolean[] filterValues = new boolean[currentTileData.numClustersInTile];
        FilterFileReader reader = new FilterFileReader(filterFileMap.get(currentTileData.tileNum));
        Iterator<AbstractIlluminaPositionFileReader.PositionInfo> positionInfoIterator = locs.iterator();
        int count = 0;
        while (reader.hasNext()) {
            filterValues[count] = reader.next();
            count++;
        }

        List<AbstractIlluminaPositionFileReader.PositionInfo> positions = new ArrayList<>();
        for (boolean filterValue : filterValues) {
            AbstractIlluminaPositionFileReader.PositionInfo info = positionInfoIterator.next();
            if (filterValue) {
                positions.add(info);
            }
        }
        this.positionInfoIterator = positions.iterator();
        cachedFilter.put(currentTileData.tileNum, filterValues);
    }

    private void cacheTile(int totalCycleCount, TileData tileData, CycleData currentCycleData) throws IOException {
        byte[] tileByteArray = new byte[tileData.compressedBlockSize];
        //we are going to explode the nibbles in to bytes to make PF filtering easier
        byte[] uncompressedByteArray = new byte[tileData.uncompressedBlockSize];
        // ByteBuffer uncompressedByteArray = ByteBuffer.allocate(tileData.uncompressedBlockSize);
        byte[] unNibbledByteArray = new byte[tileData.uncompressedBlockSize * 2];

        // Read the whole compressed block into a buffer, then sanity check the length
        InputStream stream = this.streams[totalCycleCount];
        long dataLeft = tileData.filePosition - stream.skip(tileData.filePosition);
        while (dataLeft > 0) {
            dataLeft -= stream.skip(dataLeft);
        }

        int readBytes = stream.read(tileByteArray);
        if (readBytes != tileData.compressedBlockSize) {
            throw new PicardException(String.format("Error while reading from BCL file for cycle %d. Offending file on disk is %s",
                    (totalCycleCount + 1), this.streamFiles[totalCycleCount].getAbsolutePath()));
        }

        // Uncompress the data from the buffer we just wrote - use gzip input stream to write to uncompressed buffer
        ByteArrayInputStream byteInputStream = new ByteArrayInputStream(Arrays.copyOfRange(tileByteArray, 0, readBytes));

        GZIPInputStream gzipInputStream = new GZIPInputStream(byteInputStream, uncompressedByteArray.length);
        int read;
        int totalRead = 0;
        try {
            while ((read = gzipInputStream.read(uncompressedByteArray, totalRead, uncompressedByteArray.length - totalRead)) != -1) {
                if (read == 0) break;
                totalRead += read;
            }
        } catch (EOFException eofException) {
            throw new PicardException("Unexpected end of file " + this.streamFiles[totalCycleCount].getAbsolutePath()
                    + " this file is likely corrupt or truncated. We have read "
                    + totalRead + "and were expecting to read "
                    + uncompressedByteArray.length);
        }
        if (totalRead != tileData.uncompressedBlockSize) {
            throw new PicardException(String.format("Error while decompressing from BCL file for cycle %d. Offending file on disk is %s",
                    (totalCycleCount + 1), this.streamFiles[totalCycleCount].getAbsolutePath()));
        }

        // Read uncompressed data from the buffer and expand each nibble into a full byte for ease of use
        int index = 0;
        for (byte singleByte : uncompressedByteArray) {
            unNibbledByteArray[index] = (byte) (singleByte & 0x0f);
            index++;
            unNibbledByteArray[index] = (byte) ((singleByte >> 4) & 0x0f);
            index++;
        }
        gzipInputStream.close();

        // Write buffer contents to cached tile array
        // if nonPF reads are included we need to strip them out
        if (!currentCycleData.pfExcluded) {
            boolean[] filterDatas = cachedFilter.get(tileData.tileNum);
            int sum = 0;
            for (boolean b : filterDatas) {
                sum += b ? 1 : 0;
            }
            byte[] filteredByteArray = new byte[sum];
            int filterIndex = 0;
            int basecallIndex = 0;
            for (boolean filterData : filterDatas) {
                byte readByte = unNibbledByteArray[filterIndex];
                if (filterData) {
                    filteredByteArray[basecallIndex] = readByte;
                    basecallIndex++;
                }
                filterIndex++;
            }
            cachedTile[totalCycleCount] = filteredByteArray;
        } else {
            cachedTile[totalCycleCount] = unNibbledByteArray;
        }
        cachedTilePosition[totalCycleCount] = 0;
    }
}
