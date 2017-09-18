package picard.illumina;

import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;
import picard.PicardException;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.IlluminaDataProviderFactory;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;

import java.io.File;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

abstract class BasecallsConverter<CLUSTER_OUTPUT_RECORD> {
    private static final Log log = Log.getInstance(BasecallsConverter.class);

    final Comparator<CLUSTER_OUTPUT_RECORD> outputRecordComparator;
    final int maxReadsInRamPerTile;
    final boolean demultiplex;
    final List<File> tmpDirs;
    final boolean ignoreUnexpectedBarcodes;
    final SortingCollection.Codec<CLUSTER_OUTPUT_RECORD> codecPrototype;
    final Class<CLUSTER_OUTPUT_RECORD> outputRecordClass;
    final Map<String, ? extends ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap;
    final ProgressLogger readProgressLogger = new ProgressLogger(log, 1000000, "Read");
    final ProgressLogger writeProgressLogger = new ProgressLogger(log, 1000000, "Write");
    int numThreads;
    ClusterDataConverter<CLUSTER_OUTPUT_RECORD> converter = null;

    protected final BclQualityEvaluationStrategy bclQualityEvaluationStrategy;
    protected List<Integer> tiles;
    protected final IlluminaDataProviderFactory factory;

    /**
     * @param barcodeRecordWriterMap   Map from barcode to CLUSTER_OUTPUT_RECORD writer.  If demultiplex is false, must contain
     *                                 one writer stored with key=null.
     * @param demultiplex              If true, output is split by barcode, otherwise all are written to the same output stream.
     * @param maxReadsInRamPerTile     Configures number of reads each tile will store in RAM before spilling to disk.
     * @param tmpDirs                  For SortingCollection spilling.
     * @param numProcessors            Controls number of threads.  If <= 0, the number of threads allocated is
     *                                 available cores - numProcessors.
     * @param outputRecordComparator   For sorting output records within a single tile.
     * @param codecPrototype           For spilling output records to disk.
     * @param outputRecordClass        Inconveniently needed to create SortingCollections.
     * @param ignoreUnexpectedBarcodes If true, will ignore reads whose called barcode is not found in barcodeRecordWriterMap,
     */
    BasecallsConverter(final Map<String, ? extends ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap,
                       final int maxReadsInRamPerTile,
                       final List<File> tmpDirs,
                       final SortingCollection.Codec<CLUSTER_OUTPUT_RECORD> codecPrototype,
                       final boolean ignoreUnexpectedBarcodes,
                       final boolean demultiplex,
                       final Comparator<CLUSTER_OUTPUT_RECORD> outputRecordComparator,
                       final BclQualityEvaluationStrategy bclQualityEvaluationStrategy,
                       final Class<CLUSTER_OUTPUT_RECORD> outputRecordClass,
                       final int numProcessors,
                       final IlluminaDataProviderFactory factory) {

        this.barcodeRecordWriterMap = barcodeRecordWriterMap;
        this.maxReadsInRamPerTile = maxReadsInRamPerTile;
        this.tmpDirs = tmpDirs;
        this.codecPrototype = codecPrototype;
        this.ignoreUnexpectedBarcodes = ignoreUnexpectedBarcodes;
        this.demultiplex = demultiplex;
        this.outputRecordComparator = outputRecordComparator;
        this.bclQualityEvaluationStrategy = bclQualityEvaluationStrategy;
        this.outputRecordClass = outputRecordClass;
        this.factory = factory;


        if (numProcessors == 0) {
            this.numThreads = Runtime.getRuntime().availableProcessors();
        } else if (numProcessors < 0) {
            this.numThreads = Runtime.getRuntime().availableProcessors() + numProcessors;
        } else {
            this.numThreads = numProcessors;
        }
    }

    public IlluminaDataProviderFactory getFactory() {
        return factory;
    }

    public abstract void doTileProcessing();

    /**
     * Must be called before doTileProcessing.  This is not passed in the ctor because often the
     * IlluminaDataProviderFactory is needed in order to construct the converter.
     *
     * @param converter Converts ClusterData to CLUSTER_OUTPUT_RECORD
     */
    void setConverter(final ClusterDataConverter<CLUSTER_OUTPUT_RECORD> converter) {
        this.converter = converter;
    }

    void setTileLimits(final Integer firstTile, final Integer tileLimit) {
        if (firstTile != null) {
            int i;
            for (i = 0; i < tiles.size(); ++i) {
                if (tiles.get(i).intValue() == firstTile.intValue()) {
                    tiles = tiles.subList(i, tiles.size());
                    break;
                }
            }
            if (tiles.get(0).intValue() != firstTile.intValue()) {
                throw new PicardException("firstTile=" + firstTile + ", but that tile was not found.");
            }
        }
        if (tileLimit != null && tiles.size() > tileLimit) {
            tiles = tiles.subList(0, tileLimit);
        }
    }

    interface ClusterDataConverter<OUTPUT_RECORD> {
        /**
         * Creates the OUTPUT_RECORDs from the cluster
         */
        OUTPUT_RECORD convertClusterToOutputRecord(final ClusterData cluster);
    }

    interface ConvertedClusterDataWriter<OUTPUT_RECORD> {
        void write(final OUTPUT_RECORD rec);

        void close();
    }

    /**
     * A comparator for tile numbers, which are not necessarily ordered by the number's value.
     */
    public static final Comparator<Integer> TILE_NUMBER_COMPARATOR = (integer1, integer2) -> {
        final String s1 = integer1.toString();
        final String s2 = integer2.toString();
        // Because a the tile number is followed by a colon, a tile number that
        // is a prefix of another tile number should sort after. (e.g. 10 sorts after 100).
        if (s1.length() < s2.length()) {
            if (s2.startsWith(s1)) {
                return 1;
            }
        } else if (s2.length() < s1.length() && s1.startsWith(s2)) {
            return -1;
        }
        return s1.compareTo(s2);
    };
}
