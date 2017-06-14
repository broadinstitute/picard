package picard.illumina.parser.readers;

import picard.util.UnsignedTypeUtil;

import java.io.File;
import java.nio.ByteBuffer;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class EmpiricalPhasingMetricsOutReader implements Iterator<EmpiricalPhasingMetricsOutReader.IlluminaPhasingMetrics> {
    private static final int HEADER_SIZE = 2;
    private static final int RECORD_SIZE = 16;

    private final BinaryFileIterator<ByteBuffer> bbIterator;

    public EmpiricalPhasingMetricsOutReader(File phasingMetricsOutFile) {
        bbIterator = MMapBackedIteratorFactory.getByteBufferIterator(HEADER_SIZE, RECORD_SIZE, phasingMetricsOutFile);
        bbIterator.getHeaderBytes();
    }

    @Override
    public boolean hasNext() {
        return bbIterator.hasNext();
    }

    @Override
    public IlluminaPhasingMetrics next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        return new IlluminaPhasingMetrics(bbIterator.next());
    }

    public static class IlluminaPhasingMetrics {
        public final TileMetricsOutReader.IlluminaLaneTileCode laneTileCode;
        public final int cycle;
        public final float phasingWeight;
        public final float prephasingWeight;

        IlluminaPhasingMetrics(final ByteBuffer bb) {
            this.laneTileCode = new TileMetricsOutReader.IlluminaLaneTileCode(UnsignedTypeUtil.uShortToInt(bb.getShort()), bb.getInt(), 0);
            this.cycle = UnsignedTypeUtil.uShortToInt(bb.getShort());
            this.phasingWeight = bb.getFloat();
            this.prephasingWeight = bb.getFloat();
        }
    }
}
