/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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
package net.sf.picard.illumina.parser;

import net.sf.picard.PicardException;

import java.util.Iterator;
import java.util.List;
import java.util.Arrays;
import java.io.File;

/**
 * Abstract base class for parser for Illumina RTA CIF raw intensity and CNF noise files.
 * These files are produced per cycle, with all the clusters
 * in a tile in a single file, so in order to get the values for a read, multiple files must be accessed.
 *
 * All that the concrete classes need to do is to specify the file type via ctor argument (cif or cnf), and
 * specify (via implementation of abstract method) which slot in IlluminaEndData are to receive
 * the FourChannelIntensityData.
 *
 * @author alecw@broadinstitute.org
 */
public abstract class IlluminaCycleFileSetParser implements IlluminaParser {

    // Describe relationship between cycle numbers and read ends and barcode
    private final ReadConfiguration readConfiguration;

    // Location of illumina output files to be parsed.  Typically this is Data/Intensities/L00<lane>
    private final File laneDirectory;
    private final int lane;
    // Iterator for all the files for the lane.
    private final CycleFileSetIterator cycleFileSetIterator;
    // List of readers for the current tile, in cycle order.
    private final ClusterIntensityFileReader[] readers;
    // iteration index
    private int cluster = 0;

    /**
     * Prepare to iterate over a set of ClusterIntensityFiles
     * @param readConfiguration describes the correspondence btw cycle number and read end & offset.
     * @param directory Data/Intensities directory containing the L00<lane> subdirectories.
     * @param lane lane of interest.
     * @param fileType cnf or cif.
     */
    public IlluminaCycleFileSetParser(final ReadConfiguration readConfiguration, final File directory, final int lane,
                                      final ClusterIntensityFileReader.FileType fileType,
                                      final List<Integer> tiles) {
        this.lane = lane;
        this.readConfiguration = readConfiguration;
        this.laneDirectory = new File(directory, "L00" + this.lane);
        this.cycleFileSetIterator = new CycleFileSetIterator(this.laneDirectory, this.lane,
                fileType, readConfiguration.getMaxCycleNumber(), tiles);
        readers = new ClusterIntensityFileReader[this.cycleFileSetIterator.getNumberOfCycleFiles()];
        if (!this.cycleFileSetIterator.hasNext()) {
            return;
        }
        getNextReaders();
    }

    /**
     * Jump so that the next record parsed will be from the specified tile.
     */
    @Override
    public void seekToTile(final int oneBasedTileNumber) {
        cycleFileSetIterator.seekToTile(oneBasedTileNumber);
        getNextReaders();
    }

    /**
     * Read the next read's set of data and set it into the provided data object.  The object must have
     * the appropriate IlluminaEndData objects set into it for first end, second end, barcode.
     *
     * Note that this implementation is coded based on the assumption that a ClusterIntensityFile has a single cycle
     * in it.  It will work if that is not the case, but could be more efficiently coded.
     */
    @Override
    public void next(final ClusterData data) {
        // Create the appropriate FourChannelIntensityData objects to receive the data being parsed.
        initializeFCIDs(data);
        for (final ClusterIntensityFileReader reader : readers) {
            for (int cycle = reader.getFirstCycle(); cycle < reader.getFirstCycle() + reader.getNumCycles(); ++cycle) {
                final ReadData readData = data.getEnd(readConfiguration.getEndTypeForCycle(cycle));
                final FourChannelIntensityData fcid = getFCID(readData);
                final int offset = readConfiguration.getOffsetForCycle(cycle);
                for (final IntensityChannel channel : IntensityChannel.values()) {
                    fcid.getChannel(channel)[offset] = reader.getValue(cluster, channel, cycle);
                }
            }
        }
        ++cluster;
        if (cluster == readers[0].getNumClusters() && cycleFileSetIterator.hasNext()) {
            cluster = 0;
            getNextReaders();
        }
    }

    @Override
    public boolean hasNext() {
        return readers.length > 0 && cluster < readers[0].getNumClusters();
    }

    /**
     * Get the set of readers for reading the next tile of reads.
     */
    private void getNextReaders() {
        final Iterator<TiledIlluminaFile> it = cycleFileSetIterator.next();
        for (int i = 0; i < readers.length; ++i) {
            final TiledIlluminaFile tiledFile = it.next();
            readers[i] = new ClusterIntensityFileReader(tiledFile.file);
        }
        if (it.hasNext()) {
            throw new PicardException("More cycle files than expected: " + it.next().file);
        }
    }

    /**
     * Create FourChannelIntensityData objects as appropriate for paired-end or single-end, barcoded or not,
     * and stuff into the right slots according to the concrete class.
     */
    private void initializeFCIDs(final ClusterData data) {
        setFCID(data.getFirstEnd(), new FourChannelIntensityData(readConfiguration.getFirstLength()));
        if (readConfiguration.isPairedEnd()) {
            setFCID(data.getSecondEnd(), new FourChannelIntensityData(readConfiguration.getSecondLength()));
        }
        if (readConfiguration.isBarcoded()) {
            setFCID(data.getBarcodeRead(), new FourChannelIntensityData(readConfiguration.getBarcodeLength()));
        }
    }

    // These two methods abstract out whether the object is to read noise or raw intensities

    /**
     * Stuff the FourChannelIntensityData into the right slot for the concrete class.
     * @param readData Object into which to stuff the FourChannelIntensityData.
     * @param fourChannelIntensityData thing to stuff into the right slot of end.
     */
    protected abstract void setFCID(ReadData readData, FourChannelIntensityData fourChannelIntensityData);

    /**
     * Get the FourChannelIntensityData from the right slot for the concrete class.
     * @param readData Object from which to get the FourChannelIntensityData.
     * @return object from the appropriate slot for the concrete class.
     */
    protected abstract FourChannelIntensityData getFCID(ReadData readData);

    /**
     * Determine if files exist for the given lane and tile.
     * @param directory
     * @param lane
     * @param tile
     * @return True if CIF files exist for the given intensity directory, lane and tile.  It is not necessarily
     * the case that all appropriats CIFs exist, rather that at least one exists.
     */
    protected static boolean filesExist(final File directory, final int lane,
                                        final ClusterIntensityFileReader.FileType fileType, final int tile) {
        try {
            CycleFileSetIterator fileSetIterator = new CycleFileSetIterator(new File(directory, "L00" + lane), lane,
                    fileType, 1, Arrays.asList(tile));
            return fileSetIterator.hasNext();
        } catch (IlluminaFileNotFoundException e) {
            return false;
        }
    }
}
