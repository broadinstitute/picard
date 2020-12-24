package picard.illumina.parser.readers;

import htsjdk.samtools.Defaults;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.IOUtil;
import picard.PicardException;
import picard.illumina.parser.BclData;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

public class BaseBclReader {
    private static final byte BASE_MASK = 0x0003;
    private static final byte[] BASE_LOOKUP = new byte[]{'A', 'C', 'G', 'T'};
    public static final byte NO_CALL_BASE = (byte) 'N';
    final InputStream[] streams;
    final File[] streamFiles;
    final int[] outputLengths;
    final int numReads;
    private final BclQualityEvaluationStrategy bclQualityEvaluationStrategy;
    final int[] numClustersPerCycle;
    final int cycles;

    /* Array of base values and quality values that are used to decode values from BCLs efficiently. */
    private static final byte[] BCL_BASE_LOOKUP = new byte[256];
    private static final byte[] BCL_QUAL_LOOKUP = new byte[256];

    static {
        BCL_BASE_LOOKUP[0] = NO_CALL_BASE;
        BCL_QUAL_LOOKUP[0] = BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY;

        for (int i=1; i<256; ++i) {
            // TODO: If we can remove the use of BclQualityEvaluationStrategy then in the lookup we
            // TODO: can just set the QUAL to max(2, (i >>> 2)) instead.
            BCL_BASE_LOOKUP[i] = BASE_LOOKUP[i & BASE_MASK];
            BCL_QUAL_LOOKUP[i] = (byte) Math.max(0, (i >>> 2));
        }
    }

    BaseBclReader(int[] outputLengths, BclQualityEvaluationStrategy bclQualityEvaluationStrategy) {
        this.outputLengths = outputLengths;
        this.numReads = outputLengths.length;
        this.bclQualityEvaluationStrategy = bclQualityEvaluationStrategy;

        int cycles = 0;
        for (final int outputLength : outputLengths) {
            cycles += outputLength;
        }

        this.cycles = cycles;
        this.streams = new InputStream[cycles];
        this.streamFiles = new File[cycles];
        this.numClustersPerCycle = new int[cycles];
    }

    BaseBclReader(int[] outputLengths) {
        this(outputLengths, null);
    }

    int getNumClusters() {
        return numClustersPerCycle[0];
    }

    InputStream open(final File file, final boolean seekable, final boolean isGzip, final boolean isBgzf) {
        final String filePath = file.getAbsolutePath();

        try {
            // Open up a buffered stream to read from the file and optionally wrap it in a gzip stream if necessary
            if (isBgzf) {
                // Only BlockCompressedInputStreams can seek, and only if they are fed a SeekableStream.
                return new BlockCompressedInputStream(IOUtil.maybeBufferedSeekableStream(file));
            } else if (isGzip) {
                if (seekable) {
                    throw new IllegalArgumentException(
                            String.format("Cannot create a seekable reader for gzip bcl: %s.", filePath)
                    );
                }
                return (IOUtil.maybeBufferInputStream(new GZIPInputStream(new FileInputStream(file), Defaults.BUFFER_SIZE / 2),
                        Defaults.BUFFER_SIZE / 2));
            } else {
                if (seekable) {
                    throw new IllegalArgumentException(
                            String.format("Cannot create a seekable reader for provided bcl: %s.", filePath)
                    );
                }
                return IOUtil.maybeBufferInputStream(new FileInputStream(file));
            }
        } catch (final FileNotFoundException fnfe) {
            throw new PicardException("File not found: (" + filePath + ")", fnfe);
        } catch (final IOException ioe) {
            throw new PicardException("Error reading file: (" + filePath + ")", ioe);
        }
    }

    final void decodeBasecall(final BclData bclData, final int read, final int cycle, final int byteToDecode) {
        bclData.bases[read][cycle] = BCL_BASE_LOOKUP[byteToDecode];
        bclData.qualities[read][cycle] = this.bclQualityEvaluationStrategy.reviseAndConditionallyLogQuality(BCL_QUAL_LOOKUP[byteToDecode]);
    }

    void decodeQualityBinnedBasecall(final BclData bclData, final int read, final int cycle, final int byteToDecode,
                                     final CycleData cycleData) {
        bclData.bases[read][cycle] = BCL_BASE_LOOKUP[byteToDecode];
        if (byteToDecode == 0) {
            bclData.qualities[read][cycle] = BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY;
        } else {
            bclData.qualities[read][cycle] = cycleData.qualityBins[byteToDecode >>> 2];
        }
    }

    public class CycleData {
        final short version;
        final int headerSize;

        final byte bitsPerBasecall;
        final byte bitsPerQualityScore;
        final int numberOfBins;
        final byte[] qualityBins;
        final int numTiles;
        final TileData tileInfo;
        final boolean pfExcluded;

        CycleData(final short version, final int headerSize, final byte bitsPerBasecall,
                  final byte bitsPerQualityScore, final int numberOfBins, final byte[] qualityBins,
                  final int numTiles, final TileData tileInfo, final boolean pfExcluded) {
            this.version = version;
            this.headerSize = headerSize;
            this.bitsPerBasecall = bitsPerBasecall;
            this.bitsPerQualityScore = bitsPerQualityScore;
            this.numberOfBins = numberOfBins;
            this.qualityBins = qualityBins;
            this.numTiles = numTiles;
            this.tileInfo = tileInfo;
            this.pfExcluded = pfExcluded;
        }

        public TileData getTileInfo() {
            return tileInfo;
        }
    }


    public class TileData {
        final int tileNum;
        final int numClustersInTile;
        final int uncompressedBlockSize;
        final int compressedBlockSize;
        final long filePosition;

        TileData(final int tileNum, final int numClustersInTile, final int uncompressedBlockSize,
                 final int compressedBlockSize, long filePosition) {
            this.tileNum = tileNum;
            this.numClustersInTile = numClustersInTile;
            this.uncompressedBlockSize = uncompressedBlockSize;
            this.compressedBlockSize = compressedBlockSize;
            this.filePosition = filePosition;
        }

        public int getTileNum() {
            return tileNum;
        }

        public int getCompressedBlockSize() {
            return compressedBlockSize;
        }

        int getNumClustersInTile() {
            return numClustersInTile;
        }

    }
}
