/*
* The MIT License
*
* Copyright (c) 2012 The Broad Institute
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
* THE SOFTWARE.
*/
package picard.illumina.parser.readers;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import picard.PicardException;
import picard.illumina.parser.BclData;
import picard.illumina.parser.TileIndex;
import picard.util.UnsignedTypeUtil;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayDeque;
import java.util.List;
import java.util.Queue;
import java.util.zip.GZIPInputStream;

/**
 * BCL Files are base call and quality score binary files containing a (base,quality) pair for successive clusters.
 * The file is structured as followed:
 * Bytes 1-4 : unsigned int numClusters
 * Bytes 5-numClusters + 5 : 1 byte base/quality score
 * <p/>
 * The base/quality scores are organized as follows (with one exception, SEE BELOW):
 * The right 2 most bits (these are the LEAST significant bits) indicate the base, where
 * A=00(0x00), C=01(0x01), G=10(0x02), and T=11(0x03)
 * <p/>
 * The remaining bytes compose the quality score which is an unsigned int.
 * <p/>
 * EXCEPTION: If a byte is entirely 0 (e.g. byteRead == 0) then it is a no call, the base
 * becomes '.' and the Quality becomes 2, the default illumina masking value
 * <p/>
 * (E.g. if we get a value in binary of 10001011 it gets transformed as follows:
 * <p/>
 * Value read: 10001011(0x8B)
 * <p/>
 * Quality     Base
 * <p/>
 * 100010      11
 * 00100010    0x03
 * 0x22        T
 * 34          T
 * <p/>
 * So the output base/quality will be a (T/34)
 */
public class BclReader extends BaseBclReader implements CloseableIterator<BclData> {
    private static final int HEADER_SIZE = 4;
    private static final int QUEUE_SIZE = 128;
    private final Queue<BclData> queue  = new ArrayDeque<>(QUEUE_SIZE);
    private final byte[] buffer         = new byte[QUEUE_SIZE];

    public BclReader(final List<File> bclsForOneTile, final int[] outputLengths,
                     final BclQualityEvaluationStrategy bclQualityEvaluationStrategy, final boolean seekable) {
        super(outputLengths, bclQualityEvaluationStrategy);
        try {

            final ByteBuffer byteBuffer = ByteBuffer.allocate(HEADER_SIZE);
            byteBuffer.order(ByteOrder.LITTLE_ENDIAN);

            for (int i = 0; i < cycles; ++i) {
                final File bclFile = bclsForOneTile.get(i);
                if (bclFile == null) {
                    close();
                    throw new RuntimeIOException(String.format("Could not find BCL file for cycle %d", i));
                }
                final String filePath = bclFile.getName();
                final boolean isGzip = filePath.endsWith(".gz");
                final boolean isBgzf = filePath.endsWith(".bgzf");
                final InputStream stream = open(bclFile, seekable, isGzip, isBgzf);
                final int read = stream.read(byteBuffer.array());
                if (read != HEADER_SIZE) {
                    close();
                    throw new RuntimeIOException(String.format("BCL %s has invalid header structure.", bclFile.getAbsoluteFile()));
                }
                numClustersPerCycle[i] = byteBuffer.getInt();
                if (!isBgzf && !isGzip) {
                    assertProperFileStructure(bclFile, numClustersPerCycle[i], stream);
                }
                this.streams[i] = stream;
                this.streamFiles[i] =  bclFile;
                byteBuffer.clear();
            }
        } catch (final IOException ioe) {
            throw new RuntimeIOException(ioe);
        }
    }

    public static boolean isGzipped(final File file) {
        return file.getAbsolutePath().endsWith(".gz");
    }

    public static boolean isBlockGzipped(final File file) {
        return file.getAbsolutePath().endsWith(".bgzf");
    }

    public static long getNumberOfClusters(final File file) {
        InputStream stream = null;
        try {
            if (isBlockGzipped(file)) stream = new BlockCompressedInputStream(IOUtil.maybeBufferedSeekableStream(file));
            else if (isGzipped(file)) stream = new GZIPInputStream(IOUtil.maybeBufferInputStream(new FileInputStream(file)));
            else stream = IOUtil.maybeBufferInputStream(new FileInputStream(file));

            return getNumberOfClusters(file.getAbsolutePath(), stream);

        } catch (final IOException ioe) {
            throw new PicardException("Could not open file " + file.getAbsolutePath() + " to get its cluster count: " + ioe.getMessage(), ioe);
        } finally {
            CloserUtil.close(stream);
        }
    }

    private static long getNumberOfClusters(final String filePath, final InputStream inputStream) {
        final byte[] header = new byte[HEADER_SIZE];

        try {
            final int headerBytesRead = inputStream.read(header);
            if (headerBytesRead != HEADER_SIZE) {
                throw new PicardException("Malformed file, expected header of size " + HEADER_SIZE + " but received " + headerBytesRead);
            }
        } catch (final IOException ioe) {
            throw new PicardException("Unable to read header for file (" + filePath + ")", ioe);
        }

        final ByteBuffer headerBuf = ByteBuffer.wrap(header);
        headerBuf.order(ByteOrder.LITTLE_ENDIAN);
        return UnsignedTypeUtil.uIntToLong(headerBuf.getInt());
    }


    public BclReader(final File bclFile, final BclQualityEvaluationStrategy bclQualityEvaluationStrategy, final boolean seekable) {
        super(new int[]{1}, bclQualityEvaluationStrategy);
        try {
            final ByteBuffer byteBuffer = ByteBuffer.allocate(HEADER_SIZE);
            final String filePath = bclFile.getName();
            final boolean isGzip = filePath.endsWith(".gz");
            final boolean isBgzf = filePath.endsWith(".bgzf");
            final InputStream stream = open(bclFile, seekable, isGzip, isBgzf);
            final int read = stream.read(byteBuffer.array());

            if (read != HEADER_SIZE) {
                throw new RuntimeIOException(String.format("BCL %s has invalid header structure.", bclFile.getAbsoluteFile()));
            }

            byteBuffer.order(ByteOrder.LITTLE_ENDIAN);
            this.numClustersPerCycle[0] = byteBuffer.getInt();
            if (!isBgzf && !isGzip) {
                assertProperFileStructure(bclFile, this.getNumClusters(), stream);
            }
            this.streams[0] = stream;
            this.streamFiles[0] = bclFile;
        } catch (final IOException ioe) {
            throw new PicardException("IOException opening file " + bclFile.getAbsoluteFile(), ioe);
        }
    }

    void assertProperFileStructure(final File file, final int numClusters, final InputStream stream) {
        final long elementsInFile = file.length() - HEADER_SIZE;
        if (numClusters != elementsInFile) {
            CloserUtil.close(stream);
            throw new PicardException("Expected " + numClusters + " in file " + file.getAbsolutePath() + " but found " + elementsInFile);
        }
    }

    public void close() {
        for (final InputStream stream : this.streams) {
            CloserUtil.close(stream);
        }
    }

    @Override
    public boolean hasNext() {
        if (queue.isEmpty()) {
            advance();
        }

        return !queue.isEmpty();
    }

    protected void assertProperFileStructure(final File file) {
        final long elementsInFile = file.length() - HEADER_SIZE;
        if (getNumClusters() != elementsInFile) {
            throw new PicardException("Expected " + getNumClusters() + " in file " + file.getAbsolutePath() + " but found " + elementsInFile);

        }
    }

    public BclData next() {
        if (!hasNext()) throw new IllegalStateException("next() called on BclReader that has no more items.");
        return queue.remove();
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    void advance() {
        int totalCycleCount = 0;

        try {
            // See how many clusters we can read and then make BclData objects for them
            final int clustersRead = this.streams[0].read(buffer);
            if (clustersRead == -1) return;

            final BclData[] bclDatas = new BclData[clustersRead];
            for (int i=0; i<clustersRead; ++i) {
                bclDatas[i] = new BclData(outputLengths);
            }

            // Process the data from the first cycle since we had to read it to know how many clusters we'd get
            updateClusterBclDatas(bclDatas, 0, 0);
            totalCycleCount += 1;

            for (int read = 0; read < numReads; ++read) {
                final int readLen = this.outputLengths[read];
                final int firstCycle = (read == 0) ? 1 : 0;  // For the first read we already did the first cycle above

                for (int cycle = firstCycle; cycle < readLen; ++cycle) {
                    final int n = this.streams[totalCycleCount].read(buffer, 0, clustersRead);
                    assert n == clustersRead;
                    totalCycleCount += 1;

                    updateClusterBclDatas(bclDatas, read, cycle);
                }
            }

            for (final BclData data : bclDatas) {
              this.queue.add(data);
            }
        }
        catch (final IOException ioe) {
            throw new RuntimeIOException(String.format("Error while reading from BCL file for cycle %d. Offending file on disk is %s",
                                    (totalCycleCount+1), this.streamFiles[totalCycleCount].getAbsolutePath()), ioe);
        }
    }

    /** Inserts the bases and quals at `cycle` of `read` in all of the bclDatas, using data in this.buffer. */
    private void updateClusterBclDatas(final BclData[] bclDatas, final int read, final int cycle) {
        final int numClusters = bclDatas.length;
        for (int dataIdx = 0; dataIdx< numClusters; ++dataIdx) {
            final BclData data = bclDatas[dataIdx];
            final int b = Byte.toUnsignedInt(this.buffer[dataIdx]);
            decodeBasecall(data, read, cycle, b);
        }
    }

    public static BclReader makeSeekable(final List<File> files, final BclQualityEvaluationStrategy bclQualityEvaluationStrategy, final int[] outputLengths) {
        return new BclReader(files, outputLengths, bclQualityEvaluationStrategy, true);
    }

    public int seek(final List<File> files, final TileIndex tileIndex, final int currentTile) {
        int count = 0;
        int numClustersInTile = 0;
        for (final InputStream inputStream : streams) {
            final TileIndex.TileIndexRecord tileIndexRecord = tileIndex.findTile(currentTile);
            final BclIndexReader bclIndexReader = new BclIndexReader(files.get(count));
            final long virtualFilePointer = bclIndexReader.get(tileIndexRecord.getZeroBasedTileNumber());
            if (!(inputStream instanceof BlockCompressedInputStream)) {
                throw new UnsupportedOperationException("Seeking only allowed on bzgf");
            } else {
                try {
                    if (tileIndex.getNumTiles() != bclIndexReader.getNumTiles()) {
                        throw new PicardException(String.format("%s.getNumTiles(%d) != %s.getNumTiles(%d)",
                                tileIndex.getFile().getAbsolutePath(), tileIndex.getNumTiles(), bclIndexReader.getBciFile().getAbsolutePath(), bclIndexReader.getNumTiles()));
                    }
                    ((BlockCompressedInputStream) inputStream).seek(virtualFilePointer);
                    numClustersInTile = tileIndexRecord.getNumClustersInTile();
                } catch (final IOException e) {
                    throw new PicardException("Problem seeking to " + virtualFilePointer, e);
                }
            }
            count++;
        }
        return numClustersInTile;
    }
}
