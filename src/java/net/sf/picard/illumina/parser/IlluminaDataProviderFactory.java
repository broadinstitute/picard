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

import java.io.File;
import java.util.*;

import net.sf.picard.PicardException;

/**
 * IlluminaDataProviderFactory accepts options for parsing Illumina data files for a lane and creates an
 * IlluminaDataProvider, an iterator over the ClusterData for that lane, which utilizes these options.
 *
 *
 * Note: Since we tend to use IlluminaDataProviderFactory in multithreaded environments (e.g. we call makeDataProvider
 * in a different thread per tile in IlluminaBasecallsToSam).  I've made it essentially immutable.  makeDataProvider/getTiles
 * are now idempotent (well as far as IlluminaDataProviderFactory is concerned, many file handles and other things are
 * opened when makeDataProvider is called).  We may in the future want dataTypes to be provided to the
 * makeDataProvider factory methods so configuration is not done multiple times for the same basecallDirectory in
 * client code.
 *
 * @author jburke@broadinstitute.org
 */
public class IlluminaDataProviderFactory {
    // The following properties must be specified by caller.
    /** basecallDirectory holds QSeqs or bcls **/
    private final File basecallDirectory;
    private final int lane;

    //rawIntensity directory is assumed to be the parent of basecallDirectory
    private final File rawIntensityDirectory;

    /**
     * BarcodeCycle and length are not available in output QSeqs and must be passed by the user,
     * these are used to determine the ReadStructure and will later be completely replaced by
     * ReadStructure.
     * */
    private final Integer barcodeCycle;
    private final Integer barcodeLength;

    /** The types of data that will be returned by any IlluminaDataProviders created by this factory.
     *
     * Note: In previous version, data of types not specified might be returned if a data type was specified
     * for data residing in QSeqs (since QSeqs span multiple data types).  This is no longer the case, you
     * MUST specify all data types that should be returned.*/
    private final Set<IlluminaDataType> dataTypes;

    /**
     * Soon to be removed, currently used in auto-detecting ReadStructure
     * numQSeq is the number of QSeq "end" files for this lane
     */
    private int numQSeq;

    /** readStructure is the computed readStructure past to individual parsers */
    private final ReadStructure readStructure;

    /**
     * Prepare to iterate over non-barcoded Illumina Basecall output.
     *
     * @param basecallDirectory Where the Illumina basecall output is.  Some files are found by looking relative to this
     *                          directory, at least by default.
     * @param lane              Which lane to iterate over.
     * @param readStructure     If ReadStructure is specified then the IlluminaDataProvider produced by this factory will
     *                          produce clusters conforming to readStructure, otherwise the readStructure is detected by the factory
     * @param dataTypes         Which data types to read.
     */
    public IlluminaDataProviderFactory(final File basecallDirectory, final int lane, final ReadStructure readStructure,
                                final IlluminaDataType... dataTypes) {
        this.basecallDirectory     = basecallDirectory;
        this.rawIntensityDirectory = basecallDirectory.getParentFile();
        
        this.lane = lane;
        this.barcodeCycle = null;
        this.barcodeLength = null;
        this.dataTypes = new HashSet<IlluminaDataType>(Arrays.asList(dataTypes));

        if(readStructure == null) {
            this.readStructure = computeReadStructure();
        } else {
            this.readStructure = readStructure;
        }
    }

    /**
     * Prepare to iterate over the Illumina basecall output with barcodes.
     *
     * @param basecallDirectory Where the basecall output is.  Some files are found by looking relative to this
     *                          directory, at least by default.
     * @param lane              Which lane to iterate over.
     * @param barcodeCycle 1-base cycle number where barcode starts
     * @param barcodeLength number of cycles in barcode
     * @param dataTypes         Which data types to read.
     */
    public IlluminaDataProviderFactory(final File basecallDirectory, final int lane,
                                       final int barcodeCycle, final int barcodeLength,
                                final IlluminaDataType... dataTypes) {
        this.basecallDirectory     = basecallDirectory;
        this.rawIntensityDirectory = basecallDirectory.getParentFile();

        this.lane = lane;
        if (barcodeCycle < 1) {
            throw new IllegalArgumentException("barcodeCycle is < 1: " + barcodeCycle);
        }
        if (barcodeLength < 1) {
            throw new IllegalArgumentException("barcodeLength is < 1: " + barcodeLength);
        }
        this.barcodeCycle = barcodeCycle;
        this.barcodeLength = barcodeLength;
        this.dataTypes = new HashSet<IlluminaDataType>(Arrays.asList(dataTypes));
        this.readStructure = computeReadStructure();
    }

    /**
     * Return the list of tiles available for this flowcell and lane.  These are in ascending numerical order.
     *
     * @return List of all tiles available for this flowcell and lane.
     */
    public List<Integer> getTiles() {
        return new ArrayList<Integer>(IlluminaFileUtil.getEndedIlluminaBasecallFiles(basecallDirectory, "qseq", lane, 1).keySet());
    }

    public ReadStructure readStructure() {
        return readStructure;
    }

    /**
     * Call this method to create an iterator over all clusters for all tiles in ascending numeric order.
     *
     * @return An iterator for reading the Illumina basecall output for the lane specified in the ctor.
     */
    public IlluminaDataProvider makeDataProvider() {
        return makeDataProvider(null);
    }
    /**
     * Call this method to create an iterator over the specified tiles.
     *
     * @param tiles The tiles to iterate over, in order, or null to get all tiles in ascending numerical order.
     * @return An iterator for reading the Illumina basecall output for the lane specified in the constructor.
     */
    public IlluminaDataProvider makeDataProvider(final List<Integer> tiles) {
        if (dataTypes.isEmpty()) {
            throw new PicardException("No data types have been specified for basecall output " + basecallDirectory +
            ", lane " + lane);
        }

        final int [] outputLengths = new int[readStructure.descriptors.size()];
        for(int i = 0; i < outputLengths.length; i++) {
            outputLengths[i] = readStructure.descriptors.get(i).length;
        }

        final int totalCycles = readStructure.totalCycles;
        if(!dataTypes.contains(IlluminaDataType.Position)) {
            boolean addPosition = false;
            for(final IlluminaDataType dt : dataTypes) {
                switch(dt) {
                    case BaseCalls:
                    case PF:
                    case QualityScores:
                        addPosition = true;
                        break;
                    default: //do nothing
                }
            }

            if(addPosition) {
                dataTypes.add(IlluminaDataType.Position);
            }
        }

        final Map<IlluminaParser, Set<IlluminaDataType>> parsersToDataType = new HashMap<IlluminaParser, Set<IlluminaDataType>>();
        final SortedSet<IlluminaDataType> toSupport = new TreeSet<IlluminaDataType>(dataTypes);
        while(!toSupport.isEmpty()) {
            final IlluminaParser parser = getPreferredParser(toSupport.first(), tiles, totalCycles, outputLengths);
            final Set<IlluminaDataType> providedTypes = new HashSet<IlluminaDataType>(parser.supportedTypes());
            providedTypes.retainAll(toSupport); //provide only those we want
            toSupport.removeAll(providedTypes); //remove the ones we want from the set still needing support
            parsersToDataType.put(parser, providedTypes);
        }

        for(final IlluminaParser parser : parsersToDataType.keySet()) {
            parser.verifyData(readStructure, tiles);
        }

        return new IlluminaDataProvider(readStructure, parsersToDataType, basecallDirectory, lane);
    }

    /**
     * In the future their may be multiple parsers for the same IlluminaDataType (e.g. BCLParser and QSeqParser).  Instantiate an instance of the preferred parser for
     * the given data type and return it.
     * @param dataType The type of data we want to parse
     * @param tiles The tiles over which we will be parsing data
     * @param totalCycles The total number of cycles per tiles
     * @param outputLengths The expected arrangement of output data
     * @return A parser that will parse dataType data over the given tiles and cycles and output it in groupings of the sizes specified in outputLengths
     */
    private IlluminaParser getPreferredParser(final IlluminaDataType dataType, final List<Integer> tiles, final int totalCycles, final int [] outputLengths) {
        final IlluminaParser parser;
        switch (dataType) {
            case Barcodes:
                parser = new BarcodeParser(lane, IlluminaFileUtil.getNonEndedIlluminaBasecallFiles(basecallDirectory, "barcode", lane, tiles));
                break;
            case Noise:
                final CycleIlluminaFileMap cnfFileMap = IlluminaFileUtil.getCyledIlluminaFiles(rawIntensityDirectory, "cnf", lane, tiles, totalCycles);
                parser = new CnfParser(rawIntensityDirectory, lane, cnfFileMap, outputLengths);
                break;
            case QualityScores:
            case PF:
            case BaseCalls:
            case Position:
                if(numQSeq == 0) {
                    numQSeq = IlluminaFileUtil.getNumberOfIlluminaEnds(basecallDirectory, "qseq", lane);
                }
                final List<IlluminaFileMap> readTileMap = new ArrayList<IlluminaFileMap>();
                for(int i = 0; i < numQSeq; i++) {
                       readTileMap.add(IlluminaFileUtil.getEndedIlluminaBasecallFiles(basecallDirectory, "qseq", lane, i+1, tiles));
                }
                parser = new QseqParser(lane, outputLengths, readTileMap);
                break;
            case RawIntensities:
                final CycleIlluminaFileMap cifFileMap = IlluminaFileUtil.getCyledIlluminaFiles(rawIntensityDirectory, "cif", lane, tiles, totalCycles);
                parser = new CifParser(rawIntensityDirectory, lane, cifFileMap, outputLengths);
                break;
            default:
                throw new PicardException("Unrecognized data type(" + dataType + ") found by IlluminaDataProviderFactory!");
        }
        return parser;
    }

    /**
     * Based on the number of QSeqs in basecallDirectory and the barcode cycle/length provided, create the ReadStructure
     * that will be passed to each IlluminaDataProvider created by this factory.
     * @return An ReadStructure that fits the parameters specified by this factory and the QSeqs in basecallDirectory
     */
    private ReadStructure computeReadStructure() {
       numQSeq = IlluminaFileUtil.getNumberOfIlluminaEnds(basecallDirectory, "qseq", lane);
        if(numQSeq == 0) {
            throw new PicardException("Zero Qseqs found for lane " + lane);
        }

        final int[] qSeqLengths = new int[numQSeq];

        int totalReadLength = 0;
        for(int i = 0; i < numQSeq; i++) {
            final IlluminaFileMap files = IlluminaFileUtil.getEndedIlluminaBasecallFiles(basecallDirectory, "qseq", lane, i+1);
            qSeqLengths[i] = QseqReadParser.getReadLength(files.firstEntry().getValue());
            totalReadLength += qSeqLengths[i];
        }

        final List<ReadDescriptor> readDescriptors = new ArrayList<ReadDescriptor>();
        if(barcodeCycle != null) {
            if(barcodeCycle <= 0) {
                throw new PicardException("Barcode cycle( " + barcodeCycle + ") must be greater than 0.");
            } else {
                if(barcodeCycle != 1) {
                    readDescriptors.add(new ReadDescriptor(barcodeCycle - 1, ReadType.Template));
                }
                
                readDescriptors.add(new ReadDescriptor(barcodeLength, ReadType.Barcode));

                final int remainingLength = totalReadLength - (barcodeCycle - 1) - barcodeLength;
                if(remainingLength > 0) {
                    readDescriptors.add(new ReadDescriptor(remainingLength, ReadType.Template));
                }
            }
        } else {
            if(numQSeq > 2) {
                throw new PicardException("More than 2 qseqs(" + numQSeq + ") for non-indexed run!");
            } else {
                for(int i = 0; i < qSeqLengths.length; i++) {
                    readDescriptors.add(new ReadDescriptor(qSeqLengths[i], ReadType.Template));
                }
            }
        }

        return new ReadStructure(readDescriptors);
    }
}
