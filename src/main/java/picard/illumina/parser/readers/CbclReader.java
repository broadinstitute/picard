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

    private byte[][] cachedTile;
    private final int[] cachedTilePosition;

    private CbclData queue = null;
    private Iterator<AbstractIlluminaPositionFileReader.PositionInfo> positionInfoIterator;
    private final CycleData[] cycleData;
    private final Map<Integer, File> filterFileMap;
    private final Map<Integer, List<Boolean>> cachedFilter = new HashMap<>();
    private final Map<Integer, Map<Integer, File>> surfaceToTileToCbclMap;
    private int headerSize;
    private final Map<Integer, List<TileData>> allTiles = new HashMap<>();
    private final int[] outputCycles;

    private static final int INITIAL_HEADER_SIZE = 6;
    private static final Log log = Log.getInstance(CbclReader.class);
    private static final Pattern PATTERN = Pattern.compile("^.+C(\\d{1,4}).+L(\\d{1,3})_(\\d).cbcl$");

    public CbclReader(final List<File> cbcls, final Map<Integer, File> filterFileMap, final int[] outputLengths,
                      final int tileNum, final List<AbstractIlluminaPositionFileReader.PositionInfo> locs, final int[] outputCycles, final boolean headerOnly) {
        super(outputLengths);
        if (!filterFileMap.containsKey(tileNum)) {
            throw new PicardException("Filter file for tile " + tileNum + " does not exist.");
        }
        this.outputCycles = outputCycles;

        surfaceToTileToCbclMap = sortCbcls(cbcls);
        this.filterFileMap = filterFileMap;
        cycleData = new CycleData[cycles];
        cachedTile = new byte[cycles][];
        cachedTilePosition = new int[cycles];
        for (int i = 1; i <= cycles; i++) {
            allTiles.put(i, new ArrayList<>());
        }
        try {
            readSurfaceTile(tileNum, locs, headerOnly);
        } finally {
            close();
        }
    }

    private void readSurfaceTile(final int tileNum, final List<AbstractIlluminaPositionFileReader.PositionInfo> locs,
                                 final boolean headerOnly) {
        log.info("Processing tile " + tileNum);
        try {
            for (final Map.Entry<Integer, Map<Integer, File>> entry : surfaceToTileToCbclMap.entrySet()) {
                final Map<Integer, File> cycleMap = entry.getValue();
                for (int i = 0; i < cycles; i++) {
                    //cycleMap is 1 indexed
                    final ByteBuffer byteBuffer = ByteBuffer.allocate(INITIAL_HEADER_SIZE);
                    byteBuffer.order(ByteOrder.LITTLE_ENDIAN);

                    final File bclFile = cycleMap.get(outputCycles[i]);
                    if (bclFile == null) {
                        throw new PicardException("Expected cbcl file for surface " + entry.getKey() + " cycle " + (i + 1) + " but it was not found.");
                    }

                    final InputStream stream = open(bclFile, false, false, false);
                    int read = stream.read(byteBuffer.array());

                    //we need to read the first 6 bytes to determine the header size
                    if (read != INITIAL_HEADER_SIZE) {
                        throw new RuntimeIOException(String.format("BCL %s has invalid header structure.", bclFile.getAbsoluteFile()));
                    }

                    final short version = byteBuffer.getShort();
                    headerSize = byteBuffer.getInt();

                    final ByteBuffer headerBuffer = ByteBuffer.allocate(headerSize - INITIAL_HEADER_SIZE);
                    headerBuffer.order(ByteOrder.LITTLE_ENDIAN);

                    read = stream.read(headerBuffer.array());
                    if (read != headerSize - INITIAL_HEADER_SIZE) {
                        throw new PicardException(String.format("BCL %s has invalid header structure.", bclFile.getAbsoluteFile()));
                    }

                    final byte bitsPerBasecall = headerBuffer.get();
                    final byte bitsPerQualityScore = headerBuffer.get();

                    if (bitsPerBasecall != 2 && bitsPerBasecall != bitsPerQualityScore) {
                        throw new PicardException("CBCL data not encoded in nibbles. (not currently supported) bitsPerBasecall : "
                                + bitsPerBasecall + " bitsPerQualityScore : " + bitsPerQualityScore);
                    }

                    final int numberOfBins = headerBuffer.getInt();

                    final byte[] qualityBins = new byte[numberOfBins];
                    //each bin has a pair of 4 byte mappings
                    for (int j = 0; j < numberOfBins; j++) {
                        headerBuffer.getInt(); // first int is "from" value, which we don't need
                        final int to = headerBuffer.getInt();
                        qualityBins[j] = (byte) to;
                    }
                    long filePos = 0;

                    final int numTiles = headerBuffer.getInt();
                    TileData tileInfo = null;
                    for (int j = 0; j < numTiles; j++) {
                        final int tile = headerBuffer.getInt();
                        final int numClustersInTile = headerBuffer.getInt();
                        final int uncompressedBlockSize = headerBuffer.getInt();
                        final int compressedBlockSize = headerBuffer.getInt();
                        final TileData tileData = new TileData(tile, numClustersInTile, uncompressedBlockSize, compressedBlockSize, filePos);
                        allTiles.get(i + 1).add(tileData);
                        if (tile == tileNum) {
                            tileInfo = tileData;
                        }
                        filePos += compressedBlockSize;
                    }

                    final boolean pfExcluded = headerBuffer.get() == 1;
                    //try the next surface if we didn't find the tile
                    if (tileInfo == null) {
                        continue;
                    }

                    cycleData[i] = new CycleData(version, headerSize, bitsPerBasecall, bitsPerQualityScore, numberOfBins, qualityBins, numTiles, tileInfo, pfExcluded);
                    this.streams[i] = stream;
                    this.streamFiles[i] = bclFile;
                    byteBuffer.clear();
                    headerBuffer.clear();
                }
            }

            if (headerOnly) {
                return;
            }

            int totalCycleCount = 0;

            if (cycleData[totalCycleCount].tileInfo == null) {
                throw new PicardException("Could not find tile " + tileNum);
            }

            for (final int outputLength : outputLengths) {
                for (int cycle = 0; cycle < outputLength; cycle++) {
                    final CycleData currentCycleData = cycleData[totalCycleCount];
                    try {
                        if (cachedTile[totalCycleCount] == null) {
                            if (!cachedFilter.containsKey(cycleData[totalCycleCount].tileInfo.tileNum)) {
                                cacheFilterAndLocs(cycleData[totalCycleCount].tileInfo, locs);
                            }
                            cacheTile(totalCycleCount, cycleData[totalCycleCount].tileInfo, currentCycleData);
                        }
                    } catch (final IOException e) {
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

    private Map<Integer, Map<Integer, File>> sortCbcls(final List<File> cbcls) {
        final Map<Integer, Map<Integer, File>> sortedMap = new TreeMap<>();
        for (final File cbcl : cbcls) {
            final Matcher matcher = PATTERN.matcher(cbcl.getAbsolutePath());
            if (!matcher.matches()) {
                throw new PicardException("CBCL File " + cbcl.getAbsolutePath() + " does not match expected pattern.");
            }
            final Integer surface = Integer.valueOf(matcher.group(3));
            final Integer cycle = Integer.valueOf(matcher.group(1));
            if (sortedMap.containsKey(surface)) {
                sortedMap.get(surface).put(cycle, cbcl);
            } else {
                final Map<Integer, File> cycleMap = new HashMap<>();
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
        return queue != null;
    }

    public CbclData next() {
        if (queue == null) {
            advance();
        }
        final CbclData data = queue;
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
        final CbclData data = new CbclData(outputLengths, cycleData[totalCycleCount].tileInfo.tileNum);

        for (int read = 0; read < outputLengths.length; read++) {
            for (int cycle = 0; cycle < outputLengths[read]; cycle++) {
                final CycleData currentCycleData = cycleData[totalCycleCount];

                if (cachedTilePosition[totalCycleCount] >= cachedTile[totalCycleCount].length
                        || cachedTilePosition[totalCycleCount] >= cycleData[totalCycleCount].getTileInfo().getNumClustersInTile()) {
                    // end of tile
                    return;
                }

                final int singleByte = cachedTile[totalCycleCount][cachedTilePosition[totalCycleCount]++];

                decodeQualityBinnedBasecall(data, read, cycle, singleByte, currentCycleData);

                totalCycleCount++;
            }
        }
        data.setPositionInfo(positionInfoIterator.next());
        this.queue = data;
    }

    private void cacheFilterAndLocs(final TileData currentTileData, final List<AbstractIlluminaPositionFileReader.PositionInfo> locs) {
        final List<Boolean> filterValues = new ArrayList<>();
        final FilterFileReader reader = new FilterFileReader(filterFileMap.get(currentTileData.tileNum));
        final Iterator<AbstractIlluminaPositionFileReader.PositionInfo> positionInfoIterator = locs.iterator();

        while (reader.hasNext()) {
            filterValues.add(reader.next());
        }

        final List<AbstractIlluminaPositionFileReader.PositionInfo> positions = new ArrayList<>();
        for (final boolean filterValue : filterValues) {
            final AbstractIlluminaPositionFileReader.PositionInfo info = positionInfoIterator.next();
            if (filterValue) {
                positions.add(info);
            }
        }
        this.positionInfoIterator = positions.iterator();
        cachedFilter.put(currentTileData.tileNum, filterValues);
    }

    private void cacheTile(final int totalCycleCount, final TileData tileData, final CycleData currentCycleData) throws IOException {
        final byte[] tileByteArray = new byte[tileData.compressedBlockSize];

        // Read the whole compressed block into a buffer, then sanity check the length
        final InputStream stream = this.streams[totalCycleCount];
        long dataLeft = tileData.filePosition - stream.skip(tileData.filePosition);
        while (dataLeft > 0) {
            dataLeft -= stream.skip(dataLeft);
        }

        final int readBytes = stream.read(tileByteArray);
        if (readBytes != tileData.compressedBlockSize) {
            throw new PicardException(String.format("Error while reading from BCL file for cycle %d. Offending file on disk is %s",
                    (totalCycleCount + 1), this.streamFiles[totalCycleCount].getAbsolutePath()));
        }

        // Decompress the data from the buffer we just wrote - use gzip input stream to write to uncompressed buffer
        final ByteArrayInputStream byteInputStream = new ByteArrayInputStream(Arrays.copyOfRange(tileByteArray, 0, readBytes));
        byte[] decompressedByteArray = decompressTile(totalCycleCount, tileData, byteInputStream);

        // Read uncompressed data from the buffer and expand each nibble into a full byte for ease of use
        byte[] unNibbledByteArray = promoteNibblesToBytes(decompressedByteArray);
        cachedTile[totalCycleCount] = filterNonPfReads(tileData, currentCycleData, unNibbledByteArray);

        cachedTilePosition[totalCycleCount] = 0;
    }

    private byte[] filterNonPfReads(TileData tileData, CycleData currentCycleData, byte[] unNibbledByteArray) {
        // Write buffer contents to cached tile array
        // if nonPF reads are included we need to strip them out
        if (!currentCycleData.pfExcluded) {
            final List<Boolean> filterDatas = cachedFilter.get(tileData.tileNum);
            int sum = 0;
            for (final boolean b : filterDatas) {
                sum += b ? 1 : 0;
            }
            final byte[] filteredByteArray = new byte[sum];
            int filterIndex = 0;
            int basecallIndex = 0;
            for (final boolean filterData : filterDatas) {
                final byte readByte = unNibbledByteArray[filterIndex];
                if (filterData) {
                    filteredByteArray[basecallIndex] = readByte;
                    basecallIndex++;
                }
                filterIndex++;
            }
            return filteredByteArray;
        } else {
            return unNibbledByteArray;
        }
    }

    private byte[] promoteNibblesToBytes(byte[] decompressedByteArray) {
        //we are going to explode the nibbles in to bytes to make PF filtering easier
        final byte[] unNibbledByteArray = new byte[decompressedByteArray.length * 2];
        int index = 0;
        for (final byte singleByte : decompressedByteArray) {
            unNibbledByteArray[index] = (byte) (singleByte & 0x0f);
            index++;
            unNibbledByteArray[index] = (byte) ((singleByte >> 4) & 0x0f);
            index++;
        }
        return unNibbledByteArray;
    }

    private byte[] decompressTile(int totalCycleCount, TileData tileData, ByteArrayInputStream byteInputStream) throws IOException {
        final byte[] decompressedByteArray = new byte[tileData.uncompressedBlockSize];
        //only decompress the data if we are expecting data.
        if (decompressedByteArray.length == 0) {
            log.warn("Ignoring tile " + tileData.tileNum + " there are no PF reads.");
        } else {
            int read;
            int totalRead = 0;
            try (GZIPInputStream gzipInputStream = new GZIPInputStream(byteInputStream, decompressedByteArray.length)) {
                while ((read = gzipInputStream.read(decompressedByteArray, totalRead, decompressedByteArray.length - totalRead)) != -1) {
                    if (read == 0) {
                        break;
                    }
                    totalRead += read;
                }
            } catch (final EOFException eofException) {
                throw new PicardException("Unexpected end of file " + this.streamFiles[totalCycleCount].getAbsolutePath()
                        + " this file is likely corrupt or truncated. We have read "
                        + totalRead + " and were expecting to read "
                        + decompressedByteArray.length);
            }
            if (totalRead != tileData.uncompressedBlockSize) {
                throw new PicardException(String.format("Error while decompressing from BCL file for cycle %d. Offending file on disk is %s",
                        (totalCycleCount + 1), this.streamFiles[totalCycleCount].getAbsolutePath()));
            }
        }
        return decompressedByteArray;
    }

    public CycleData[] getCycleData() {
        return cycleData;
    }

    public int getHeaderSize() {
        return headerSize;
    }

    public List<File> getFilesForCycle(final int i) {
        final List<File> cbclFiles = new ArrayList<>();
        surfaceToTileToCbclMap.values().forEach(map -> {
            if (map.containsKey(i)) {
                cbclFiles.add(map.get(i));
            }
        });
        return cbclFiles;
    }

    public Map<Integer, List<TileData>> getAllTiles() {
        return allTiles;
    }

    public void clear() {
        cachedTile = null;
    }
}
