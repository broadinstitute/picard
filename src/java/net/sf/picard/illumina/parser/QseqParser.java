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
package net.sf.picard.illumina.parser;

import net.sf.picard.util.*;
import net.sf.picard.PicardException;

import java.io.File;
import java.util.*;
import java.util.Map.Entry;
import static net.sf.picard.illumina.parser.IlluminaDataType.*;

/**
 * QSeqParser takes a List<IlluminaFileMaps>(1 for each "end" in order) and provides an iterator over all clusters in the
 * for all tiles found in the maps.  "Clusters" in this sense refers to the chunks delineated by outputLengths.  Position, PF, Base, and Quality data
 * is output in chunks specified by these lengths.
 *
 * Each line contains the following fields:
 *     lane
 *     tile
 *     x-coordinate
 *     y-coordinate
 *     machine
 *     run
 *     passing filter
 *     bases
 *     quals
 *
 * Bases, quals and PF are extracted from this file.
 *
 *
 * @author jburke@broadinstitute.org
 */

class QseqParser implements IlluminaParser<QseqReadData> {
    private final int [] outputLengths; //Expected lengths of the bases/qualities arrays in the output QSeqReadData
    private final QseqReadParser[] parsers; //each parser parsers iterates over all tiles for 1 read
    private final List<IlluminaFileMap> tileMapByReadNumber;
    
    public static final Set<IlluminaDataType> SUPPORTED_TYPES = Collections.unmodifiableSet(CollectionUtil.makeSet(Position, BaseCalls, QualityScores, PF));

    public QseqParser(final int lane, final List<IlluminaFileMap> tileMapByReadNumber, final OutputMapping outputMapping) {
        this.outputLengths       = outputMapping.getOutputReadLengths();
        this.tileMapByReadNumber = tileMapByReadNumber;
        List<QseqReadParser> parsersList = makeReadParserList(lane, outputMapping);
        parsers = parsersList.toArray(new QseqReadParser[parsersList.size()]);
    }

    /**
     * QseqParser delegates the read/write operations of Qseq data providing to individual QseqReadParsers (1 per read that doesn't entirely
     * consist of skips).  To do this QseqParser makes a list of QseqReadParsers that each handles all the tile files for ONE specific READ of the run.
     * (E.g. if in lane 1 there are tiles 1,2,3 in the run and reads 1 and 2 then there will be ONE QseqReadParser that handles
     *  the tile files s_1_1_0001_qseq.txt, s_1_1_0002_qseq.txt, s_1_1_0003_qseq.txt and ONE QseqReadParser that handles
     *  the tile files s_1_2_0001_qseq.txt, s_1_2_0002_qseq.txt, s_1_2_0003_qseq.txt).
     *
     * A QseqReadParser reads the bases/qualities and other info from the qseq file of the current tile being handled
     * and maps them into an QseqReadData object.  However, the output bases/quals do not necessarily
     * correspond 1 to 1 with the bases/quals found in the Qseq files.  Therefore we create a mapping from
     * input ranges over the bases of a qseq file (RR below) to the location in the output arrays (LO) that the
     * bases in the range RR are to be copied to.
     *
     *
     * Assuming the reads below represent the structure of bases in the set of QSeq files in tileMapByReadNumber and
     * outputLengths = [Read1.length, Read2.length, Read3.length] - The lengths found in qseqs
     * outputRanges  = [(3,5),(8,13),(18,20)]                     - The Indexes(Cycle#-1) of desired output cycles
     *
     *          Read 1                                 Read2                      Read 3
     *     [B0, B1, B2, B3, B4, B5, B6, B7]     [B0, B1, B2, B3]     [B0, B1, B2, B3, B4,  B5,  B6,  B7]
     *      S   S  [O0, O1, O2]  S  S  [O0,      O1, O2, O3, O4,      O5, 06]  S   S   S  [O0,  O1,  O2]
     *  OCI 0   1   2   3   4   5   6   7         8   9  10  11       12  13  14  15  16   17   18   19
     *
     *  S   = Cycles skipped in the output
     *  Ox  = cycles to output at array index x of the output array delineated by the []
     *  OCI = Index of cycles in the run as a whole (i.e. cycle# - 1)
     *
     *  R(x,y) = Range from x to y inclusive
     *  T(ar,ele) = 2D array index where, assuming you have an array int [][]arr; you would index arr[ar][ele]
     *  For the reads above:
     *  RR for Read 1 = [R(2,4),R(7,7)]    LO for Read 1 = [T(0,0), T(1,0)]
     *  RR for read 2 = [R(0,4)]           LO for Read 2 = [T(1,1)]
     *  RR for read 3 = [R(0,2), R(5,7)]   LO for Read 3 = [T(1,5), T(2,0)]
     *
     * For each cluster, each read does something similar to:
     *
     * for(int i = 0; i < RR.length; i++) {
     *      Range R = RR[i]
     *      TwoDArrayIndex L = LO[i]
     *      System.arrayCopy(qseqData, R.start, R.length, outputArray[L.ar], L.ele)
     * }
     *
     * @param lane The current lane being parsed
     * @param outputMapping Information on what cycles are being output to ClusterData Objects
     * @return A list of QseqReadParsers that do the actual reading of qseq data and populating cluster data objects
     */
    public List<QseqReadParser> makeReadParserList(final int lane, final OutputMapping outputMapping) {
        //Since there is not a one-to-one relationship with ReadDescriptors in the ReadStructure and QSeq reads,
        //we have to inspect all reads for there sizes and exclude those we don't need in the parser itself
        final int [] readLengths = new int[tileMapByReadNumber.size()];
        for(int i = 0; i < readLengths.length; i++) {
            readLengths[i] = QseqReadParser.getReadLength(tileMapByReadNumber.get(i).firstEntry().getValue());
        }

        //Parsers required to fill the output cycles specified by OutputMapping
        final List<QseqReadParser> parsersList = new ArrayList<QseqReadParser>(readLengths.length);

        //split expected output ranges so that no one output range spans multiple Qseq files but all cycles in outputMapping
        //are still in at least one range
        final List<Range> splitOutputRanges         = splitOutputRangesOnQseqBoundaries(readLengths, outputMapping);

        //create a list of ranges(1 per read) of input indexes into the base/qual section of the qseq for each outputRages
        final List<List<Range>> inputRangesPerRead  = outputRangesToInputRanges(readLengths, splitOutputRanges);

        //just a list of the number of input/output range pairs there are for each read
        final List<Integer> outputRangesPerRead = new ArrayList<Integer>(inputRangesPerRead.size());
        for(final List<Range> readRange : inputRangesPerRead) {
            outputRangesPerRead.add(readRange.size());
        }

        //convert output ranges from cycle# to output indices, i.e find the exact indices into the 2D output array (on the cluster data object) that signifies
        // the start of that output range
        final List<List<OutputMapping.TwoDIndex>> outputTargets = outputRangesTo2DTargetsPerRead(outputRangesPerRead, splitOutputRanges, outputMapping);

        //create a QseqReadParser for each read that has at least one inputRange to outputTarget pair
        for(int i = 0; i < inputRangesPerRead.size(); i++) {
            if(outputRangesPerRead.get(i) > 0) {
                final List<Range> inputRanges = inputRangesPerRead.get(i);
                final List<OutputMapping.TwoDIndex> outputTarget = outputTargets.get(i);

                parsersList.add(new QseqReadParser(lane, tileMapByReadNumber.get(i),
                                                   inputRanges.toArray(new Range[inputRanges.size()]),
                                                   outputTarget.toArray(new OutputMapping.TwoDIndex[outputTarget.size()]), outputMapping));
            }
        }

        if(parsersList.size() == 0) {
            throw new PicardException("0 Qseq \"ends\" were found to parse!");
        }

        return parsersList;
    }

    /**
     * Given a list of outputRanges in which NO RANGE SPANS MULTIPLE QSEQ READS, convert these ranges into
     * input indexes for the given reads
     *(E.g.)
     *                     Read 1                        Read 2                     Read3
     * QSeq Reads[B0, B1, B2, B3, B4, B5, B6, B7]    [B0, B1, B2, B3]    [B0, B1, B2, B3, B4, B5, B6, B7]
     * Output Ranges  [ Range(1, 4)] [Range(5,7)]    [  Range(8,11) ]  [R(12,13)]    [   Range(15,19)   ]
     * ReadRanges =  List(
     *                  List(Range(1, 4), Range(5,7)),
     *                  List(Range(0, 4)),
     *                  List(Range(0, 1), Range(3,7))
     *               )
     *
     * @param readLengths The length of qseq reads
     * @param outputRanges Ranges of  cycle index numbers indicating cycles to output
     * @return A list of Ranges of indices into qseq reads , one per each outputRange segregated by read
     */
    public static List<List<Range>> outputRangesToInputRanges(final int [] readLengths, final List<Range> outputRanges) {
        final List<Range> scratchOutputRanges = new ArrayList<Range>(outputRanges);
        final List<List<Range>> readRanges = new ArrayList<List<Range>>(readLengths.length);
        for(int i = 0; i < readLengths.length; i++) {
            readRanges.add(new LinkedList<Range>());
        }

        int currentRead = 0;
        int startingCycle = 0;
        int endCycle = readLengths[startingCycle] - 1;
        List<Range> currentReadList = readRanges.get(currentRead);


        while(!scratchOutputRanges.isEmpty()) {
            final Range range = scratchOutputRanges.remove(0);
            while(range.start > endCycle) {
                //advance to next read
                currentRead    += 1;
                startingCycle   = endCycle + 1;
                endCycle        = startingCycle + readLengths[currentRead] - 1;
                currentReadList = readRanges.get(currentRead);
            }

            final int outRangeEnd;
            if(endCycle < range.end) {
                scratchOutputRanges.add(0, new Range(endCycle + 1, range.end));
                outRangeEnd = endCycle;
            } else{
                outRangeEnd = range.end;
            }

            currentReadList.add(new Range(range.start - startingCycle, outRangeEnd - startingCycle));

        }

        return readRanges;
    }

    /**
     * For each outputRange in cycle indices (cycle #s -1), translate the range into a starting "target" in
     * the 2D output array of bases/clusters in the output array
     * @param rangesPerRead An ordered list with one element per read where each element is the number of ranges for the corresponding read
     * @param splitOutputRanges A flat list of output ranges where no single range spans multiple reads
     * @param outputMapping Information on what cycles should be output and how they should be arranged
     * @return An ordered list with one element per read where each element is a list of starting points
     * to copy input ranges to
     */ //TODO: Add tests
    public static List<List<OutputMapping.TwoDIndex>> outputRangesTo2DTargetsPerRead(final List<Integer> rangesPerRead, final List<Range> splitOutputRanges, OutputMapping outputMapping) {
        List<List<OutputMapping.TwoDIndex>> outputTargets = new ArrayList<List<OutputMapping.TwoDIndex>>(rangesPerRead.size());
        for(int i = 0; i < rangesPerRead.size(); i++) {
            outputTargets.add(new ArrayList<OutputMapping.TwoDIndex>(rangesPerRead.get(i)));
        }

        int count = 0;
        int readIndex = 0;
        int numRangesForRead = rangesPerRead.get(readIndex);
        List<OutputMapping.TwoDIndex> currentReadTargets = outputTargets.get(readIndex);

        for(final Range or : splitOutputRanges) {

            while(count >= numRangesForRead) {
                count      = 0;
                readIndex += 1;
                numRangesForRead = rangesPerRead.get(readIndex);
                currentReadTargets = outputTargets.get(readIndex);
            }

            currentReadTargets.add(outputMapping.getOutputIndexForCycle(or.start+1));
            count += 1;
        }
        return outputTargets;
    }

    /**
     * Given a List of output cycle ranges with some that may range over multiple qseq read files, return an identical list except all
     * ranges that span multiple qseq read files will be split up and replaced with ranges that span over only 1 qseq file.
     *
     *                     Read 1                        Read 2                     Read3
     * QSeq Reads[B0, B1, B2, B3, B4, B5, B6, B7]    [B0, B1, B2, B3]    [B0, B1, B2, B3, B4, B5, B6, B7]
     * Ranges         [ Range(1, 4)] [                 Range(5,13)             ]     [   Range(15,19)   ]
     * Output Ranges  [ Range(1, 4)] [Range(5,7)]    [  Range(8,11) ]  [R(12,13)]    [   Range(15,19)   ]
     * @param qseqReadLengths
     * @param outputMapping
     * @return
     */
    public static List<Range> splitOutputRangesOnQseqBoundaries(int[] qseqReadLengths, OutputMapping outputMapping) {
        List<Range> outputCycleIndexRanges = new LinkedList<Range>(Arrays.asList(outputMapping.getCycleIndexRanges()));
        List<Range> outputRanges           = new LinkedList<Range>();

        int qseqArray = 0;
        int startCycleIndex = 0;
        while(outputCycleIndexRanges.size() > 0) {
            final Range curRange = outputCycleIndexRanges.remove(0);
            final int endCycleIndex = startCycleIndex + qseqReadLengths[qseqArray] - 1;
            if(curRange.start <= endCycleIndex) {
                if(curRange.end > endCycleIndex) {
                    outputRanges.add(new Range(curRange.start, endCycleIndex));
                    outputCycleIndexRanges.add(0, new Range(endCycleIndex+1, curRange.end));
                } else {
                    outputRanges.add(curRange);
                }
            } else {
                outputCycleIndexRanges.add(0, curRange);
                startCycleIndex = endCycleIndex + 1;
                ++qseqArray;
            }
        }

        return outputRanges;
    }

    /**
     * Assert that all ends have the same number of tiles and that the tiles for a given end have the same number of bases (checking only the first cluster
     * per file).  In the binary formats it will be easier to check the number of clusters etc but for Qseqs any problems with number of clusters will be
     * discovered when we iterate past the end of one of the reads.
     */
    @Override
    public void verifyData(final List<Integer> tiles, final int [] cycles) { //in this case, runConfiguration related information was verified in construction
        Integer numTiles = null;
        int end = 1;
        for(final IlluminaFileMap tilesToFiles : tileMapByReadNumber) {
            if(tiles != null && !tilesToFiles.keySet().containsAll(tiles)) {     //now with properly util classes this shouldn't happen but since this is only done once up front, test it anyway
                TreeSet<Integer> missingTiles = new TreeSet<Integer>(tiles);
                missingTiles.removeAll(tilesToFiles.keySet());
                
                String missing = missingTiles.first().toString();
                missingTiles.remove(missingTiles.first());
                for(final Integer tile : missingTiles) {
                    missing += ", " + tile;
                }
                throw new PicardException("IlluminaFileMap for \"end\" number " + end + " is missing tiles: " + missing);
            }

            if(numTiles == null) {  //make sure all ends have the same number of tiles
                numTiles = tilesToFiles.size();
            } else if(numTiles != tilesToFiles.size()) {
                throw new PicardException("Qseq \"end\" files do not have the same number of tiles expected(" + numTiles + ") found(" + tilesToFiles.size() + ") on end (" + end + ")");
            }

            //make sure this end has files with the same read lengths
            Integer numBases = null;
            for(final Entry<Integer, File> tileNoToFile : tilesToFiles.entrySet()) {
                final int readLength = QseqReadParser.getReadLength(tileNoToFile.getValue());
                if(numBases == null) {
                    numBases = readLength;
                } else if(numBases != readLength) {
                    throw new PicardException("Qseq \"end\" (" + end + ") has tiles with different numbers of bases per read.  Found on Tile(" + tileNoToFile.getKey() + ") File(" + tileNoToFile.getValue().getAbsolutePath() + ")");
                }
            }
            ++end;
        }
    }



    @Override
    public void seekToTile(final int oneBasedTileNumber) {
        for(final QseqReadParser parser : parsers) {
            parser.seekToTile(oneBasedTileNumber);
        }
    }

    @Override
    public QseqReadData next() {
        final QseqReadData qseqRd = new QseqReadData(outputLengths);
        for(final QseqReadParser parser : parsers) {
            parser.next(qseqRd);
        }
        return qseqRd;
    }

    public void remove() {
        throw new UnsupportedOperationException("Remove is not supported by " + QseqParser.class.getName());
    }

    @Override
    public boolean hasNext() {
        return parsers[0].hasNext();
    }

    @Override
    public Set<IlluminaDataType> supportedTypes() {
        return SUPPORTED_TYPES;
    }
}

/**
 * QSeqReadFilesParser parsers all tiles for 1 particular read (e.g. first of pair, second of pair, barcode).
 */
class QseqReadParser {
    public static final int MACHINE_COLUMN = 0;
    public static final int RUN__COLUMN = 1;
    public static final int LANE_COLUMN = 2;
    public static final int TILE_COLUMN = 3;
    public static final int X_COLUMN = 4;
    public static final int Y_COLUMN = 5;
    public static final int PF_COLUMN = 10;
    public static final int BASES_COLUMN = 8;
    public static final int QUALS_COLUMN = 9;

    private final IlluminaTextIterator textParser; //Parses multiple illumina tabbed delimited files
    private final FormatUtil formatter;
    private final Range[] srcRanges;
    private final OutputMapping.TwoDIndex[] dstStarts;
    private final int [] copyLengths;
    protected final int readLength;
    private static final SolexaQualityConverter qualConverter = SolexaQualityConverter.getSingleton();

    public QseqReadParser(final int lane, final IlluminaFileMap tilesToReadFiles, final Range [] srcRanges, final OutputMapping.TwoDIndex[] dstStarts, final OutputMapping outputMapping)
    {
        this.textParser = new IlluminaTextIterator(lane, tilesToReadFiles);
        this.formatter = new FormatUtil();
        this.readLength = QseqReadParser.getReadLength(tilesToReadFiles.firstEntry().getValue());

        this.srcRanges = srcRanges;
        this.dstStarts = dstStarts;

        copyLengths = new int[srcRanges.length];
        for(int i = 0; i < copyLengths.length; i++) {
            copyLengths[i] = srcRanges[i].length;
        }
    }

    public void seekToTile(int oneBasedTileNumber) {
        textParser.seekToTile(oneBasedTileNumber);
    }

    public boolean hasNext() {
        return textParser.hasNext();
    }

    public void next(final QseqReadData readData) {
        final String [] fields = textParser.next();

        final int lane = formatter.parseInt(fields[LANE_COLUMN]);
        textParser.validateLane(lane);

        final int tile = formatter.parseInt(fields[TILE_COLUMN]);
        final int x = formatter.parseInt(fields[X_COLUMN]);
        final int y = formatter.parseInt(fields[Y_COLUMN]);

        // Currently unused
        //final String machine = fields[MACHINE_COLUMN];
        //final int run = getFormatter().parseInt(fields[RUN__COLUMN]);
        final boolean pf = formatter.parseInt(fields[PF_COLUMN]) == 1;

        final String baseString = fields[BASES_COLUMN];
        final String qualString = fields[QUALS_COLUMN];
        if (baseString.length() != qualString.length()) {
            throw new PicardException("Length of bases and quals don't match in " + textParser.getCurrentFilename());
        }

        readData.setOrCheckLane(lane);
        readData.setOrCheckTile(tile);
        readData.setOrCheckXCoordinate(x);
        readData.setOrCheckYCoordinate(y);
        readData.setOrCheckPf(pf);

        stringToBytes(baseString, srcRanges, dstStarts, readData.getBases());
        stringToQuals(qualString, srcRanges, dstStarts, copyLengths, readData.getQualities());
    }

    private static void stringToBytes(final String s, final Range[] sourceRanges, final OutputMapping.TwoDIndex[] destRanges, byte[][] outputBuffers) {
        for(int i = 0; i < sourceRanges.length; i++) {
            s.getBytes(sourceRanges[i].start, sourceRanges[i].end + 1, outputBuffers[destRanges[i].majorIndex], destRanges[i].minorIndex);
        }
    }

    private static void stringToQuals(final String s, final Range[] sourceRanges, final OutputMapping.TwoDIndex[] destRanges, final int [] copyLengths, final byte [][] outputBuffers ) {
        stringToBytes(s, sourceRanges, destRanges, outputBuffers);

        for(int i = 0; i < destRanges.length; i++) {
            qualConverter.convertSolexa_1_3_QualityCharsToPhredBinary(destRanges[i].minorIndex, copyLengths[i], outputBuffers[destRanges[i].majorIndex]);
        }
    }

    public static int getReadLength(final File qseqFile) {
        final BasicInputParser parser = new BasicInputParser(true, qseqFile);
        if (!parser.hasNext()) {
            throw new PicardException("Unexpected empty qseq file: " + qseqFile);
        }
        final String[] fields = parser.next();
        return fields[BASES_COLUMN].length();
    }
}


class QseqReadData implements PositionalData, BaseData, QualityData, PfData {
    private final byte [][] bases;
    private final byte [][] qualities;
    private Boolean pf;
    private Integer xCoord;
    private Integer yCoord;
    private Integer lane;
    private Integer tile;

    public QseqReadData(final int [] bufferSizes) {
        bases = new byte[bufferSizes.length][];
        qualities = new byte[bufferSizes.length][];
        for(int i = 0; i < bufferSizes.length; i++) {
            bases[i] = new byte[bufferSizes[i]];
            qualities[i] = new byte[bufferSizes[i]];
        }

        pf = null;
        xCoord = null;
        yCoord = null;
        lane = null;
        tile = null;
    }

    @Override
    public byte[][] getBases() {
        return bases;
    }

    @Override
    public byte[][] getQualities() {
        return qualities;
    }

    public void setOrCheckPf(final boolean pf) {
        if(this.pf == null) {
            this.pf = pf;
        } else {
            assertEquals(this.pf, pf, "pf");
        }
    }

    @Override
    public boolean isPf() {
        return pf;
    }

    @Override
    public int getXCoordinate() {
        return xCoord;
    }

    public void setOrCheckXCoordinate(final int xCoord) {
        if(this.xCoord == null) {
            this.xCoord = xCoord;
        } else {
            assertEquals(this.xCoord, xCoord, "xCoord");
        }
    }

    @Override
    public int getYCoordinate() {
        return yCoord;
    }

    public void setOrCheckYCoordinate(final int yCoord) {
        if(this.yCoord == null) {
            this.yCoord = yCoord;
        } else {
            assertEquals(this.yCoord, yCoord, "yCoord");
        }
    }

    @Override
    public int getLane() {
        return lane;
    }

    public void setOrCheckLane(final int lane) {
        if(this.lane == null) {
            this.lane = lane;
        } else {
            assertEquals(this.lane, lane, "lane");
        }
    }

    @Override
    public int getTile() {
        return tile;
    }

    public void setOrCheckTile(final int tile) {
        if(this.tile == null) {
            this.tile = tile;
        } else {
            assertEquals(this.tile, tile, "tile");
        }
    }

    public <T> void assertEquals(T current, T newValue, String valueName) {
        if(!current.equals(newValue)) {
            throw new PicardException(valueName + " values don't match: original(" + current + ") new(" + newValue + ")");
        }
    }
}
