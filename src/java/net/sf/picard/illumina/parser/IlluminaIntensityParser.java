
package net.sf.picard.illumina.parser;

import net.sf.picard.PicardException;

import java.io.File;
import java.util.List;
import java.util.Map;

/**
 * Common parent class to CnfParser and CifParser.  Creates the relevant CycleFileParser to read Intensity files and
 * partially implements the makeData used to create either the RawIntensityData or NoiseData object that will be
 * returned in IlluminaCycledFileSetParser.next().
 *
 * @author jburke@broadinstitute.org
 */
abstract class IlluminaIntensityParser<T extends IlluminaData> extends PerTilePerCycleParser<T> {
    public IlluminaIntensityParser(final File directory, final int lane, final CycleIlluminaFileMap tilesToCycleFiles, final int [] readLengths) {
        super(directory, lane, tilesToCycleFiles, readLengths);
    }

    /**
     * Add the intensity for one cycle to the illData object at index for the given channel.
     * @param illData The IlluminaDataObject created in the makeData method for this cluster
     * @param index Where in the FCID arrays where intensity should be placed
     * @param channel A,C,G, or T
     * @param intensity The value read by the CycleFileParser for this cycle
     */
    protected abstract void addIntensityToIlluminaData(final T illData, final CompositeIndex index, final IntensityChannel channel, final short intensity);

    /**
     * Both CnfParser and CifParser produces an array of FCIDs, implement this method to convert them into the relevant IlluminaData subclass(RawIntensityData or NoiseData)
     * @param fcids The intensity data that should be returned by the IlluminaData subclass
     * @return Fcids wrapped by an IlluminaData class of type T
     */
    protected abstract T intensityToIlluminaData(final FourChannelIntensityData[] fcids);

    /**
     * Make an IlluminaData object of type T with FCIDs of the given lengths.
     * @param outputLengths The expected lengths of the output data
     * @return IlluminaData object of Type T
     */
    @Override
    protected T makeData(final int[] outputLengths) {
        final FourChannelIntensityData [] fcids = new FourChannelIntensityData[outputLengths.length];
        for(int i = 0; i < outputLengths.length; i++) {
            fcids[i] = new FourChannelIntensityData(outputLengths[i]);
        }

        return intensityToIlluminaData(fcids);
    }

    /**
     * Return an IntensityFileParser for the given file/cycle
     * @param file The file to parse
     * @param cycle The cycle that file represents
     * @return CycleFileParser
     */
    @Override
    protected CycleFileParser<T> makeCycleFileParser(final File file, final int cycle) {
        return new IntensityFileParser(file, cycle);
    }

    @Override
    public void verifyData(final ReadStructure readStructure, final List<Integer> tiles) {
        //Verify we have the correct number of cycles and for each cycle we have the correct number of tiles
        super.verifyData(readStructure, tiles);

        Long fileSize;
        Integer numClusters;
        Integer elementSize = null;

        for(final Map.Entry<Integer, CycleFilesIterator> tileToFiles : tilesToCycleFiles.entrySet()) {
            fileSize    = null;
            numClusters = null;

            for(final File intensityFile : tileToFiles.getValue()) {
                final ClusterIntensityFileReader.ClusterIntensityFileHeader header = ClusterIntensityFileReader.readHeaders(intensityFile); //catches a number of file structure errors
                if(fileSize == null) {
                    fileSize    = intensityFile.length();
                    numClusters = header.numClusters;
                    elementSize = header.elementSize;
                } else {
                    if(fileSize != intensityFile.length()) {
                        throw new PicardException("Intensity cycle files for tile(" + tileToFiles.getValue() + ") were not of the same size, offending file(" + intensityFile.getAbsolutePath() + ")");
                    }

                    if(numClusters != header.numClusters) {
                        throw new PicardException("Intensity cycle files for tile(" + tileToFiles.getValue() + ") did not have the same number of clusters(" + intensityFile.getAbsolutePath() + ")");
                    }

                    if(elementSize != header.elementSize) {
                        throw new PicardException("Intensity cycle files for tile(" + tileToFiles.getValue() + ") did not have the same element sizes(" + intensityFile.getAbsolutePath() + ")");
                    }

                    if(header.numCycles != 1) {
                        throw new PicardException("Intensity cycle file for tile(" + tileToFiles.getValue() + ") had more than one cycle per file num cycles found(" + header.numCycles + ")");
                    }
                }
            }
        }
    }

    protected class IntensityFileParser implements CycleFileParser<T>{
        private final ClusterIntensityFileReader reader;
        private final int cycle;
        private int cluster;
        private final CompositeIndex index;

        public IntensityFileParser(final File file, final int cycle) {
            reader = new ClusterIntensityFileReader(file);
            cluster = 0;
            index = getIndex(cycle);
            this.cycle = cycle;
        }

        public void close() {
        }

        public void next(final T ild) {
            for(final IntensityChannel channel : IntensityChannel.values()) {
                addIntensityToIlluminaData(ild, index, channel, reader.getValue(cluster, channel, cycle));
            }
            ++cluster;
        }

        public boolean hasNext() {
            return cluster < reader.getNumClusters();
        }
    }
}

