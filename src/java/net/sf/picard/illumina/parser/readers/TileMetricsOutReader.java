package net.sf.picard.illumina.parser.readers;

import net.sf.picard.PicardException;
import net.sf.picard.util.UnsignedTypeUtil;

import java.io.File;
import java.nio.ByteBuffer;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Reads a TileMetricsOut file commonly found in the InterOp directory of an Illumina Run Folder.  This
 * reader DOES NOT try to interpret the metrics code or metrics value but instead returns them in what
 * is essentially a struct.
 *
 * File Format:
 * byte 0 (unsigned byte) = The version number which MUST be 2 or an exception will be thrown
 * byte 1 (unsigned byte) = The record size which must be 10 or an exception will be thrown
 * bytes 3 + (current_record * 10) to (current_record * 10 + 10) (TileMetrics Record) = The actual records each of size 10 that
 *          get converted into IlluminaTileMetrics objects
 *
 * TileMetrics Record Format:
 * Each 10 byte record is of the following format:
 * byte 0-1 (unsigned short) = lane number
 * byte 2-3 (unsigned short) = tile number
 * byte 4-5 (unisgned short) = metrics code, see Theory of RTA document by Illumina for definition
 * byte 6-9 (float)          = metrics value, see Theory of RTA document by Illumina for definition
 */
public class TileMetricsOutReader implements Iterator<TileMetricsOutReader.IlluminaTileMetrics> {
    private static final int HEADER_SIZE = 2;
    private static final int EXPECTED_RECORD_SIZE = 10;
    private static final int EXPECTED_VERSION = 2;

    private final BinaryFileIterator<ByteBuffer> bbIterator;

    /**
     * Return a TileMetricsOutReader for the specified file
     * @param tileMetricsOutFile The file to read
     */
    public TileMetricsOutReader(final File tileMetricsOutFile) {
        bbIterator = MMapBackedIteratorFactory.getByteBufferIterator(HEADER_SIZE, EXPECTED_RECORD_SIZE, tileMetricsOutFile);

        final ByteBuffer header = bbIterator.getHeaderBytes();

        //Get the version, should be EXPECTED_VERSION, which is 2
        final int actualVersion = UnsignedTypeUtil.uByteToInt(header.get());
        if(actualVersion != EXPECTED_VERSION) {
            throw new PicardException("TileMetricsOutReader expects the version number to be " + EXPECTED_VERSION + ".  Actual Version in Header( " + actualVersion + ")" );
        }

        final int actualRecordSize = UnsignedTypeUtil.uByteToInt(header.get());
        if(EXPECTED_RECORD_SIZE != actualRecordSize) {
            throw new PicardException("TileMetricsOutReader expects the record size to be " + EXPECTED_RECORD_SIZE + ".  Actual Record Size in Header( " + actualRecordSize + ")" );
        }
    }

    public boolean hasNext() {
        return bbIterator.hasNext();
    }

    public IlluminaTileMetrics next() {
        if(!hasNext()) {
            throw new NoSuchElementException();
        }
        return new IlluminaTileMetrics(bbIterator.next());
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }

    /**
     * IlluminaTileMetrics corresponds to a single record in a TileMetricsOut file
     */
    public class IlluminaTileMetrics {
        private final int laneNumber;
        private final int tileNumber;
        private final int metricCode;
        private final float metricValue;

        public IlluminaTileMetrics(final ByteBuffer bb) {
            laneNumber  = UnsignedTypeUtil.uShortToInt(bb.getShort());
            tileNumber  = UnsignedTypeUtil.uShortToInt(bb.getShort());
            metricCode  = UnsignedTypeUtil.uShortToInt(bb.getShort());
            metricValue = bb.getFloat();
        }

        public int getLaneNumber() {
            return laneNumber;
        }

        public int getTileNumber() {
            return tileNumber;
        }

        public int getMetricCode() {
            return metricCode;
        }

        public float getMetricValue() {
            return metricValue;
        }
    }
}
