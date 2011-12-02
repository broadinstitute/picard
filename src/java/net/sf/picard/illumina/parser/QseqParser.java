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
    
    private static final Set<IlluminaDataType> SupportedTypes = Collections.unmodifiableSet(CollectionUtil.makeSet(Position, BaseCalls, QualityScores, PF));

    public QseqParser(final int lane, final int [] outputLengths, final List<IlluminaFileMap> tileMapByReadNumber) {
        this.outputLengths       = outputLengths;
        this.tileMapByReadNumber = tileMapByReadNumber;

        int totalExpectedLength = 0;
        for(int i = 0; i < outputLengths.length; i++) {
            totalExpectedLength += outputLengths[i];
        }

        int parsersIndex = 0;
        int writeOffset = 0;
        parsers = new QseqReadParser[tileMapByReadNumber.size()];
        for(final IlluminaFileMap tileMap : tileMapByReadNumber) {
            final QseqReadParser parser = new QseqReadParser(lane, tileMap, writeOffset, outputLengths);
            parsers[parsersIndex++] = parser;
            writeOffset += parser.readLength;
        }

        if(parsers.length == 0) {
            throw new PicardException("0 Qseq \"ends\" were found to parse!");
        }

        //This check is handled up front because we've already calculated the relevant values, verifyData will calculate a few other checks
        //In general, each parser should try to keep as much verification in verifyData as possible in because (as with the last iteration of this code
        //it tends to end up in code sprawl
        if(totalExpectedLength != writeOffset) {
            throw new PicardException("Requested read lengths(" + Arrays.toString(outputLengths) + ") span multiple clusters!  Specified read lengths do not land on Qseq cluster boundaries. " +
                    "Total expected length(" + totalExpectedLength + ") and actual length(" + writeOffset + ")");
        }
    }

    /**
     * Assert that all ends have the same number of tiles and that the tiles for a given end have the same number of bases (checking only the first cluster
     * per file).  In the binary formats it will be easier to check the number of clusters etc but for Qseqs any problems with number of clusters will be
     * discovered when we iterate past the end of one of the reads.
     *
     * @param readStructure read structure against which to verify this run
     */
    @Override
    public void verifyData(final ReadStructure readStructure, final List<Integer> tiles) { //in this case, runConfiguration related information was verified in construction
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
        return SupportedTypes;
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
    private final Range[] sourceRanges;
    private final CompositeIndex [] destRanges;
    private final int [] copyLengths;
    protected final int readLength;

    public QseqReadParser(final int lane, final IlluminaFileMap tilesToReadFiles,  int writeOffset, final int [] outputLengths)
    {
        this.textParser = new IlluminaTextIterator(lane, tilesToReadFiles);
        this.formatter = new FormatUtil();
        this.readLength = QseqReadParser.getReadLength(tilesToReadFiles.firstEntry().getValue());

        int bufferIndex = 0;
        while(writeOffset >= outputLengths[bufferIndex]) {
            writeOffset -= outputLengths[bufferIndex++];
        }

        int readOffset = 0;
        int coveredLength = 0;
        final List<Integer> lengthsList = new ArrayList<Integer>();
        final List<Range> srcRangesList = new ArrayList<Range>();
        final List<CompositeIndex> dstRangesList = new ArrayList<CompositeIndex>();
        while(coveredLength < readLength) { //TODO: Something is not being caught here if qseq length is actually > than runConfig length
            int length = Math.min(outputLengths[bufferIndex] - writeOffset, readLength - coveredLength);
            lengthsList.add(length);
            srcRangesList.add(new Range(readOffset, readOffset + length));
            dstRangesList.add(new CompositeIndex(bufferIndex, writeOffset));
            ++bufferIndex;          //The only way coveredLength < readLength is if the buffer - offset was to small to hold entire read length
            writeOffset = 0;  //so increment bufferIndex and appropriate writeOffset
            readOffset += length;
            coveredLength += length;
        }

        copyLengths = new int[lengthsList.size()];
        for(int i = 0; i < copyLengths.length; i++) {
            copyLengths[i] = lengthsList.get(i);
        }
        sourceRanges = srcRangesList.toArray(new Range[srcRangesList.size()]);
        destRanges   = dstRangesList.toArray(new CompositeIndex[dstRangesList.size()]);
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

        stringToBases(baseString, sourceRanges, destRanges, readData.getBases());
        stringToQuals(qualString, sourceRanges, destRanges, copyLengths, readData.getQualities());
    }

    private static void stringToBases(final String s, final Range[] sourceRanges, final CompositeIndex[] destRanges, byte [][] outputBuffers ) {
        for(int i = 0; i < sourceRanges.length; i++) {
            s.getBytes(sourceRanges[i].start, sourceRanges[i].end, outputBuffers[destRanges[i].arrayIndex], destRanges[i].elementIndex);
        }
    }

    private static void stringToQuals(final String s, final Range[] sourceRanges, final CompositeIndex[] destRanges, final int [] copyLengths, final byte [][] outputBuffers ) {
        stringToBases(s, sourceRanges, destRanges, outputBuffers);

        for(int i = 0; i < destRanges.length; i++) {
            SolexaQualityConverter.getSingleton().convertSolexa_1_3_QualityCharsToPhredBinary(destRanges[i].elementIndex, copyLengths[i], outputBuffers[destRanges[i].arrayIndex]);
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
