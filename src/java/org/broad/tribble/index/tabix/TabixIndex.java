/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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
package org.broad.tribble.index.tabix;

import net.sf.samtools.Bin;
import net.sf.samtools.BinningIndexContent;
import net.sf.samtools.Chunk;
import net.sf.samtools.LinearIndex;
import net.sf.samtools.util.BlockCompressedInputStream;
import net.sf.samtools.util.BlockCompressedOutputStream;
import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.util.StringUtil;
import org.broad.tribble.TribbleException;
import org.broad.tribble.index.Block;
import org.broad.tribble.index.Index;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.LittleEndianOutputStream;
import org.broad.tribble.util.TabixUtils;

import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.*;

/**
 * This class represent a Tabix index that has been built in memory or read from a file.  It can be queried or
 * written to a file.
 */
public class TabixIndex implements Index {
    private static final byte[] MAGIC = {'T', 'B', 'I', 1};
    public static final int MAGIC_NUMBER;
    static {
        final ByteBuffer bb = ByteBuffer.allocate(MAGIC.length);
        bb.put(MAGIC);
        bb.flip();
        MAGIC_NUMBER = bb.order(ByteOrder.LITTLE_ENDIAN).getInt();
    }

    private final TabixFormat formatSpec;
    private final List<String> sequenceNames;
    private final BinningIndexContent[] indices;

    /**
     * @param formatSpec Information about how to interpret the file being indexed.  Unused by this class other than
     *                   written to an output file.
     * @param sequenceNames Sequences in the file being indexed, in the order they appear in the file.
     * @param indices One for each element of sequenceNames
     */
    public TabixIndex(final TabixFormat formatSpec, final List<String> sequenceNames, final BinningIndexContent[] indices) {
        if (sequenceNames.size() != indices.length) {
            throw new IllegalArgumentException("sequenceNames.size() != indices.length");
        }
        this.formatSpec = formatSpec.clone();
        this.sequenceNames = Collections.unmodifiableList(new ArrayList<String>(sequenceNames));
        this.indices = indices;
    }

    /**
     * @param inputStream This is expected to be buffered and be gzip-decompressing as appropriate.  Caller
     *                    should close input stream after ctor returns.
     */
    public TabixIndex(final InputStream inputStream) throws IOException {
        this(inputStream, false);
    }

    /**
     * Convenient ctor that opens the file, wraps with with BGZF reader, and closes after reading index.
     */
    public TabixIndex(final File tabixFile) throws IOException {
        this(new BlockCompressedInputStream(tabixFile), true);
    }

    private TabixIndex(final InputStream inputStream, final boolean closeInputStream) throws IOException {
        final LittleEndianInputStream dis = new LittleEndianInputStream(inputStream);
        if (dis.readInt() != MAGIC_NUMBER) {
            throw new TribbleException(String.format("Unexpected magic number 0x%x", MAGIC_NUMBER));
        }
        final int numSequences = dis.readInt();
        indices = new BinningIndexContent[numSequences];
        formatSpec = new TabixFormat();
        formatSpec.flags = dis.readInt();
        formatSpec.sequenceColumn = dis.readInt();
        formatSpec.startPositionColumn = dis.readInt();
        formatSpec.endPositionColumn = dis.readInt();
        formatSpec.metaCharacter = (char)dis.readInt();
        formatSpec.numHeaderLinesToSkip = dis.readInt();
        final int nameBlockSize = dis.readInt();
        final byte[] nameBlock = new byte[nameBlockSize];
        if (dis.read(nameBlock) != nameBlockSize) throw new EOFException("Premature end of file reading Tabix header");
        final List<String> sequenceNames = new ArrayList<String>(numSequences);
        int startPos = 0;
        for (int i = 0; i < numSequences; ++i) {
            int endPos = startPos;
            while (nameBlock[endPos] != '\0') ++endPos;
            sequenceNames.add(StringUtil.bytesToString(nameBlock, startPos, endPos - startPos));
            startPos = endPos + 1;
        }
        if (startPos != nameBlockSize) {
            throw new TribbleException("Tabix header format exception.  Sequence name block is longer than expected");
        }
        for (int i = 0; i < numSequences; ++i) {
            indices[i] = loadSequence(i, dis);
        }
        if (closeInputStream) CloserUtil.close(dis);
        this.sequenceNames = Collections.unmodifiableList(sequenceNames);
    }

    /**
     *
     * @param chr the chromosome
     * @param start the start position, one-based, inclusive.
     * @param end the end position, one-based, inclusive.
     * @return List of regions of file that are candidates for the given query.
     *
     * TODO: This method has not yet been tested, since the primary task is index writing.
     */
    @Override
    public List<Block> getBlocks(final String chr, final int start, final int end) {
        final int sequenceIndex = sequenceNames.indexOf(chr);
        if (sequenceIndex == -1 || indices[sequenceIndex] == null) {
            return Collections.EMPTY_LIST;
        }
        final List<Chunk> chunks = indices[sequenceIndex].getChunksOverlapping(start, end);
        final List<Block> ret = new ArrayList<Block>(chunks.size());
        for (final Chunk chunk : chunks) {
            ret.add(new Block(chunk.getChunkStart(), chunk.getChunkEnd() - chunk.getChunkStart()));
        }
        return ret;
    }

    @Override
    public boolean isCurrentVersion() {
        return true;
    }

    @Override
    public List<String> getSequenceNames() {
        return sequenceNames;
    }

    @Override
    public boolean containsChromosome(final String chr) {
        return sequenceNames.contains(chr);
    }

    /**
     *
     * No arbitrary properties in Tabix
     */
    @Override
    public Map<String, String> getProperties() {
        return null;
    }

    @Override
    public boolean equalsIgnoreProperties(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final TabixIndex that = (TabixIndex) o;

        if (!formatSpec.equals(that.formatSpec)) return false;
        if (!Arrays.equals(indices, that.indices)) return false;
        return sequenceNames.equals(that.sequenceNames);

    }

    public TabixFormat getFormatSpec() {
        return formatSpec;
    }

    /**
     * Writes the index with BGZF.
     * @param tabixFile Where to write the index.
     */
    public void write(final File tabixFile) {
        final LittleEndianOutputStream los = new LittleEndianOutputStream(new BlockCompressedOutputStream(tabixFile));
        try {
            write(los);
            los.close();
        } catch (final IOException e) {
            throw new TribbleException("Exception writing " + tabixFile.getAbsolutePath(), e);
        }
    }

    /**
     * Writes to a file with appropriate name and directory based on feature file.
     * @param featureFile File being indexed.
     */
    @Override
    public void writeBasedOnFeatureFile(final File featureFile) throws IOException {
        if (!featureFile.isFile()) return;
        write(new File(featureFile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION));
    }

    /**
     *
     * @param los It is assumes that caller has done appropriate buffering and BlockCompressedOutputStream wrapping.
     *            Caller should close output stream after invoking this method.
     * @throws IOException
     */
    @Override
    public void write(final LittleEndianOutputStream los) throws IOException {
        los.writeInt(MAGIC_NUMBER);
        los.writeInt(sequenceNames.size());
        los.writeInt(formatSpec.flags);
        los.writeInt(formatSpec.sequenceColumn);
        los.writeInt(formatSpec.startPositionColumn);
        los.writeInt(formatSpec.endPositionColumn);
        los.writeInt(formatSpec.metaCharacter);
        los.writeInt(formatSpec.numHeaderLinesToSkip);
        int nameBlockSize = sequenceNames.size(); // null terminators
        for (final String sequenceName : sequenceNames) nameBlockSize += sequenceName.length();
        los.writeInt(nameBlockSize);
        for (final String sequenceName : sequenceNames) {
            los.write(StringUtil.stringToBytes(sequenceName));
            los.write(0);
        }
        for (final BinningIndexContent index : indices) {
            writeSequence(index, los);
        }
    }

    private void writeSequence(final BinningIndexContent indexContent, final LittleEndianOutputStream los) throws IOException {
        if (indexContent == null) {
            los.writeInt(0);
        } else {
            final BinningIndexContent.BinList binList = indexContent.getBins();
            los.writeInt(binList.numberOfNonNullBins);
            for (final Bin bin : binList) {
                writeBin(bin, los);
            }
            writeLinearIndex(indexContent.getLinearIndex(), los);
        }
    }

    private void writeLinearIndex(final LinearIndex linearIndex, final LittleEndianOutputStream los) throws IOException {
        if (linearIndex.getIndexStart() != 0) {
            // This could be handled by writing zeroes, but it is not expected so just fail.
            throw new IllegalArgumentException("Non-zero linear index start");
        }
        final long[] entries = linearIndex.getIndexEntries();
        los.writeInt(entries.length);
        for (final long entry : entries) los.writeLong(entry);
    }

    private void writeBin(final Bin bin, final LittleEndianOutputStream los) throws IOException {
        los.writeInt(bin.getBinNumber());
        final List<Chunk> chunkList = bin.getChunkList();
        los.writeInt(chunkList.size());
        for (final Chunk chunk: chunkList) {
            los.writeLong(chunk.getChunkStart());
            los.writeLong(chunk.getChunkEnd());
        }
    }

    /**
     * Although this is probably identical to BAM index reading code, code does not exist there to load directly
     * into a BinningIndexContent object, so that is implemented here.
     * @param referenceSequenceIndex Merely for setting in the returned object, not for seeking into the file.
     * @param dis This method assumes that the current position is at the start of the reference.
     */
    private BinningIndexContent loadSequence(final int referenceSequenceIndex, final LittleEndianInputStream dis) throws IOException {
        final int numBins = dis.readInt();
        if (numBins == 0) return null;
        int nonNullBins = 0;
        final ArrayList<Bin> bins = new ArrayList<Bin>();
        for (int i = 0; i < numBins; ++i) {
            final Bin bin = loadBin(referenceSequenceIndex, dis);
            if (bin != null) {
                // File is not sparse, but array being produced is sparse, so grow array with nulls as appropriate
                // so that bin number == index into array.
                ++nonNullBins;
                if (bins.size() > bin.getBinNumber()) {
                    if (bins.get(bin.getBinNumber()) != null) {
                        throw new TribbleException("Bin " + bin.getBinNumber() + " appears more than once in file");
                    }
                    bins.set(bin.getBinNumber(),bin);
                } else {
                    // Grow bins array as needed.
                    bins.ensureCapacity(bin.getBinNumber() + 1);
                    while (bins.size() < bin.getBinNumber()) bins.add(null);
                    bins.add(bin);
                }
            }
        }
        final LinearIndex linearIndex = loadLinearIndex(referenceSequenceIndex, dis);
        return new BinningIndexContent(referenceSequenceIndex,
                new BinningIndexContent.BinList(bins.toArray(new Bin[bins.size()]), nonNullBins), linearIndex);
    }

    private LinearIndex loadLinearIndex(final int referenceSequenceIndex, final LittleEndianInputStream dis) throws IOException {
        final int numElements = dis.readInt();
        final long[] elements = new long[numElements];
        for (int i = 0; i < numElements; ++i) {
            elements[i] = dis.readLong();
        }
        return new LinearIndex(referenceSequenceIndex, 0, elements);
    }

    private Bin loadBin(final int referenceSequenceIndex, final LittleEndianInputStream dis) throws IOException {
        final int binNumber = dis.readInt();
        final Bin ret = new Bin(referenceSequenceIndex, binNumber);
        final int numChunks = dis.readInt();
        final List<Chunk> chunkList = new ArrayList<Chunk>(numChunks);
        for (int i = 0; i < numChunks; ++i) {
            chunkList.add(loadChunk(dis));
        }
        ret.setChunkList(chunkList);
        return ret;
    }

    private Chunk loadChunk(final LittleEndianInputStream dis) throws IOException {
        final long start = dis.readLong();
        final long end = dis.readLong();
        return new Chunk(start, end);
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final TabixIndex index = (TabixIndex) o;

        if (!formatSpec.equals(index.formatSpec)) return false;
        if (!Arrays.equals(indices, index.indices)) return false;
        if (!sequenceNames.equals(index.sequenceNames)) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = formatSpec.hashCode();
        result = 31 * result + sequenceNames.hashCode();
        result = 31 * result + Arrays.hashCode(indices);
        return result;
    }
}
