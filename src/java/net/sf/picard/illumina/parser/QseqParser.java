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

import net.sf.picard.util.BasicInputParser;
import net.sf.picard.PicardException;
import net.sf.picard.util.IlluminaUtil;
import net.sf.samtools.util.StringUtil;

import java.io.File;
import java.util.List;

/**
 * Parser for files of the form s_<lane>_<end>_<tile>_qseq.txt(.gz)* produced by Illumina 1.3 and later.
 * Note that there is a separate file for each end, but for a barcoded run, the barcode data is pre or postpended
 * to one of the ends.  For Illumina 1.4 and later, the barcode will be in a separate qseq file.  In that case,
 * barcode could be in any of s_<lane>_[123]_<tile>_qseq.txt depending on when in the run the barcode is read.
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
 * @author alecw@broadinstitute.org
 */
public class QseqParser extends AbstractIlluminaTextParser {

    private static final int LANE_COLUMN = 2;
    private static final int TILE_COLUMN = 3;
    private static final int X_COLUMN = 4;
    private static final int Y_COLUMN = 5;
    private static final int MACHINE_COLUMN = 0;
    private static final int RUN__COLUMN = 1;
    private static final int PF_COLUMN = 10;
    private static final int BASES_COLUMN = 8;
    private static final int QUALS_COLUMN = 9;

    private final ReadType readType;

    /**
     * Prepare to parse qseq files
     * @param readConfiguration Used to determine if the end being processed contains the barcode read, and where
     * it is.
     * @param directory Where to find qseq files.
     * @param lane
     * @param readType Only FIRST and SECOND allowed.  For barcoded run with barcode-aware basecall files, use the
     * ctor where file number can be specified explicitly, because end # and file # do not necessarily correspond.
     * @param tiles List of tiles to be processed, in order, or null to process all tiles.
     */
    public QseqParser(final ReadConfiguration readConfiguration, final File directory, final int lane,
                      final ReadType readType, final List<Integer> tiles) {
        this(readConfiguration, directory, lane, readType, readType.getDefaultFileNumber(), tiles);
        if (readType != ReadType.FIRST && readType != ReadType.SECOND) {
            // Don't use this ctor for barcode-aware basecall files.
            throw new IllegalArgumentException("Invalid EndType: " + readType);
        }
    }

    /**
     * Prepare to parse qseq files
     * @param readConfiguration Used to determine if the end being processed contains the barcode read, and where
     * it is.
     * @param directory Where to find qseq files.
     * @param lane
     * @param readType Whether this is the parser for first end, second end or barcode end.
     * @param fileNumber number for this end, i.e. s__<lane>_<fileNumber>_<tile>_qseq.txt.
     * @param tiles List of tiles to be processed, in order, or null to process all tiles.
     */
    public QseqParser(final ReadConfiguration readConfiguration, final File directory, final int lane,
                      final ReadType readType, final int fileNumber, final List<Integer> tiles) {
        super(readConfiguration, lane, directory);
        if (fileNumber < 1 || fileNumber > 4) {
            throw new IllegalArgumentException("Invalid fileNumber: " + fileNumber);
        }
        this.readType = readType;
        setFiles(IlluminaFileUtil.getEndedIlluminaBasecallFiles(directory, "qseq", lane, fileNumber, tiles));
        initializeParser(0);
    }

    /**
     * Process a line of qseq input.  Set the PF value on the ClusterData.  If barcoded, and this
     * is the end containing the barcode, split the bases and quals.  Stuff the bases and quals into the
     * appropriate end (as determined by the oneBasedEnd property), and into the barcode end if appropriate.
     *
     * @param data Parsed input is stored in this object, which has already been set up appropriately
     * for paired or single end read, and for barcode or not.
     * @param fields Input line, split on whitespace.
     */
    @Override
    protected void processLine(final ClusterData data, final String[] fields) {
        final int lane = getFormatter().parseInt(fields[LANE_COLUMN]);
        validateLane(lane);
        final int tile = getFormatter().parseInt(fields[TILE_COLUMN]);
        final int x = getFormatter().parseInt(fields[X_COLUMN]);
        final int y = getFormatter().parseInt(fields[Y_COLUMN]);

        // Currently unused
        //final String machine = fields[MACHINE_COLUMN];
        //final int run = getFormatter().parseInt(fields[RUN__COLUMN]);
        final boolean pf = getFormatter().parseInt(fields[PF_COLUMN]) == 1;

        final String baseString = fields[BASES_COLUMN];
        final String qualString = fields[QUALS_COLUMN];
        if (baseString.length() != qualString.length()) {
            throw new PicardException("Length of bases and quals don't match in " + getCurrentFilename());
        }
        final int expectedLength;
        if (readType == ReadType.FIRST) {
            expectedLength = getReadConfiguration().getFirstLength();
        } else if (readType == ReadType.SECOND) {
            expectedLength = getReadConfiguration().getSecondLength();
        } else {
            expectedLength = getReadConfiguration().getBarcodeLength();
        }
        int expectedLengthIncludingBarcode = expectedLength;
        final boolean containsBarcode = getReadConfiguration().isBarcoded() &&
                getReadConfiguration().getBarcodeRead() == readType && !isBarcodeQseq();
        if (containsBarcode) {
            // Barcode is on this end
            expectedLengthIncludingBarcode += getReadConfiguration().getBarcodeLength();
        }
        if (expectedLengthIncludingBarcode != baseString.length()) {
            if (isBarcodeQseq() && getReadConfiguration().getBarcodeReads() > 1) {
                // That's ok then
            }
            else {
                throw new PicardException("Length of bases does not match expected in " + getCurrentFilename());
            }
        }

        data.setOrCheckLane(lane);
        data.setOrCheckTile(tile);
        data.setOrCheckX(x);
        data.setOrCheckY(y);
        data.setOrCheckPf(pf);

        final ReadData end = (isBarcodeQseq()? data.getBarcodeRead(): data.getEnd(readType));

        if (containsBarcode) {
            final int barcodeOffset = getReadConfiguration().getOffsetOfBarcodeInRead();
            final int readOffset = getReadConfiguration().getOffsetOfNonBarcodeInRead();
            if (expectedLength > 0) {
                // PFS-113 2nd end consists solely of barcode, so end == null
                final byte[] readBases = StringUtil.stringToBytes(baseString, readOffset, expectedLength);
                final byte[] readQuals = IlluminaUtil.makePhredBinaryFromSolexaQualityAscii_1_3(qualString, readOffset, expectedLength);
                end.setBases(readBases);
                end.setQualities(readQuals);
            }
            final byte[] barcodeBases = StringUtil.stringToBytes(baseString, barcodeOffset, getReadConfiguration().getBarcodeLength());
            final byte[] barcodeQuals = IlluminaUtil.makePhredBinaryFromSolexaQualityAscii_1_3(qualString, barcodeOffset,
                    getReadConfiguration().getBarcodeLength());
            data.getBarcodeRead().setBases(barcodeBases);
            data.getBarcodeRead().setQualities(barcodeQuals);
        } else {
            final byte[] bases = StringUtil.stringToBytes(baseString);
            final byte[] quals = IlluminaUtil.makePhredBinaryFromSolexaQualityAscii_1_3(qualString);

            if (!isBarcodeQseq() || end.getBases() == null) {
                end.setBases(bases);
                end.setQualities(quals);
            }
            else {
                end.setBases(concat(end.getBases(), bases));
                end.setQualities(concat(end.getQualities(), quals));
            }
        }
    }

    // Concatenate two arrays
    private byte[] concat(final byte[] b1, final byte[] b2) {
        final byte[] retval = new byte[b1.length + b2.length];
        System.arraycopy(b1, 0, retval, 0, b1.length);
        System.arraycopy(b2, 0, retval, b1.length, b2.length);
        return retval;
    }

    /**
     * True if this is a qseq that contains only the barcode.
     */
    private boolean isBarcodeQseq() {
        return readType == ReadType.BARCODE;
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
