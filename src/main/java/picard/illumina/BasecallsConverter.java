package picard.illumina;

import htsjdk.samtools.util.IOUtil;
import picard.PicardException;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.IlluminaDataProviderFactory;
import picard.illumina.parser.IlluminaDataType;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;

import java.io.File;
import java.util.*;
import java.util.regex.Pattern;

/**
 * BasecallsConverter utilizes an underlying IlluminaDataProvider to convert parsed and decoded sequencing data
 * from standard Illumina formats to specific output records (FASTA records/SAM records).
 * <p>
 * The underlying IlluminaDataProvider apply several optional transformations that can include EAMSS filtering,
 * non-PF read filtering and quality score recoding using a BclQualityEvaluationStrategy.
 * <p>
 * The converter can also limit the scope of data that is converted from the data provider by setting the
 * tile to start on (firstTile) and the total number of tiles to process (tileLimit).
 * <p>
 * Additionally, BasecallsConverter can optionally demultiplex reads by outputting barcode specific reads to
 * their associated writers..
 *
 * @param <CLUSTER_OUTPUT_RECORD> The type of record that this converter will convert to.
 */
public abstract class BasecallsConverter<CLUSTER_OUTPUT_RECORD> {
    public static final IlluminaDataType[] DATA_TYPES_WITH_BARCODE = {
            IlluminaDataType.BaseCalls,
            IlluminaDataType.QualityScores,
            IlluminaDataType.Position,
            IlluminaDataType.PF,
            IlluminaDataType.Barcodes
    };
    public static final IlluminaDataType[] DATA_TYPES_WITHOUT_BARCODE =
            Arrays.copyOfRange(DATA_TYPES_WITH_BARCODE, 0, DATA_TYPES_WITH_BARCODE.length - 1);

    protected final IlluminaDataProviderFactory factory;
    protected final boolean demultiplex;
    protected final boolean ignoreUnexpectedBarcodes;
    protected final Map<String, ? extends ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap;
    protected final boolean includeNonPfReads;
    protected List<Integer> tiles;
    protected ClusterDataConverter<CLUSTER_OUTPUT_RECORD> converter = null;
    protected int numThreads;

    /**
     * Constructs a new BasecallsConverter object.
     *
     * @param basecallsDir                 Where to read basecalls from.
     * @param barcodesDir                  Where to read barcodes from (optional; use basecallsDir if not specified).
     * @param lane                         What lane to process.
     * @param readStructure                How to interpret each cluster.
     * @param barcodeRecordWriterMap       Map from barcode to CLUSTER_OUTPUT_RECORD writer.  If demultiplex is false, must contain
     *                                     one writer stored with key=null.
     * @param demultiplex                  If true, output is split by barcode, otherwise all are written to the same output stream.
     * @param numThreads                   Controls number of threads.
     * @param firstTile                    (For debugging) If non-null, start processing at this tile.
     * @param tileLimit                    (For debugging) If non-null, process no more than this many tiles.
     * @param bclQualityEvaluationStrategy The basecall quality evaluation strategy that is applyed to decoded base calls.
     * @param ignoreUnexpectedBarcodes     If true, will ignore reads whose called barcode is not found in barcodeRecordWriterMap.
     * @param applyEamssFiltering          If true, apply EAMSS filtering if parsing BCLs for bases and quality scores.
     * @param includeNonPfReads            If true, will include ALL reads (including those which do not have PF set).
     *                                     This option does nothing for instruments that output cbcls (Novaseqs)
     */
    public BasecallsConverter(
            final File basecallsDir,
            final File barcodesDir,
            final int lane,
            final ReadStructure readStructure,
            final Map<String, ? extends ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap,
            final boolean demultiplex,
            final int numThreads,
            final Integer firstTile,
            final Integer tileLimit,
            final BclQualityEvaluationStrategy bclQualityEvaluationStrategy,
            final boolean ignoreUnexpectedBarcodes,
            final boolean applyEamssFiltering,
            final boolean includeNonPfReads
    ) {
        this.barcodeRecordWriterMap = barcodeRecordWriterMap;
        this.ignoreUnexpectedBarcodes = ignoreUnexpectedBarcodes;
        this.demultiplex = demultiplex;
        this.numThreads = numThreads;

        this.factory = new IlluminaDataProviderFactory(basecallsDir,
                barcodesDir, lane, readStructure, bclQualityEvaluationStrategy, getDataTypesFromReadStructure(readStructure, demultiplex));
        this.factory.setApplyEamssFiltering(applyEamssFiltering);
        this.includeNonPfReads = includeNonPfReads;
        this.tiles = factory.getAvailableTiles();
        tiles.sort(TILE_NUMBER_COMPARATOR);
        setTileLimits(firstTile, tileLimit);
    }

    /**
     * Abstract method for processing tiles of data and outputting records for each barcode.
     *
     * @param barcodes The barcodes used optionally for demultiplexing. Must contain at least a single null value if
     *                 no demultiplexing is being done.
     */
    public abstract void processTilesAndWritePerSampleOutputs(final Set<String> barcodes);

    /**
     * Interface that defines a converter that takes ClusterData and returns OUTPUT_RECORD type objects.
     *
     * @param <OUTPUT_RECORD> The recode type to convert to.
     */
    protected interface ClusterDataConverter<OUTPUT_RECORD> {
        /**
         * Creates the OUTPUT_RECORDs from the cluster
         */
        OUTPUT_RECORD convertClusterToOutputRecord(final ClusterData cluster);
    }

    /**
     * Interface that defines a writer that will write out OUTPUT_RECORD type objects.
     *
     * @param <OUTPUT_RECORD> The recode type to convert to.
     */
    protected interface ConvertedClusterDataWriter<OUTPUT_RECORD> {
        /**
         * Write out a single record of type OUTPUT_RECORD.
         *
         * @param rec The record to write.
         */
        void write(final OUTPUT_RECORD rec);

        /**
         * Closes the writer.
         */
        void close();
    }

    /**
     * A comparator used to sort Illumina tiles in their proper order.
     */
    public static final Comparator<Integer> TILE_NUMBER_COMPARATOR = (integer1, integer2) -> {
        final String s1 = integer1.toString();
        final String s2 = integer2.toString();
        // Because the tile number is followed by a colon, a tile number that
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

    /**
     * Applies an lane and tile based regex to return all files matching that regex for each tile.
     *
     * @param baseDirectory The directory to search for tiled files.
     * @param pattern       The pattern used to match files.
     * @return A file array of all of the tile based files that match the regex pattern.
     */
    public static File[] getTiledFiles(final File baseDirectory, final Pattern pattern) {
        return IOUtil.getFilesMatchingRegexp(baseDirectory, pattern);
    }

    /**
     * Given a read structure return the data types that need to be parsed for this run
     *
     * @param readStructure The read structure that defines how the read is set up.
     * @param demultiplex   If true, output is split by barcode, otherwise all are written to the same output stream.
     * @return A data type array for each piece of data needed to satisfy the read structure.
     */
    protected static IlluminaDataType[] getDataTypesFromReadStructure(final ReadStructure readStructure,
                                                                      final boolean demultiplex) {
        if (!readStructure.hasSampleBarcode() || !demultiplex) {
            return DATA_TYPES_WITHOUT_BARCODE;
        } else {
            return DATA_TYPES_WITH_BARCODE;
        }
    }

    /**
     * Gets the data provider factory used to create the underlying data provider.
     *
     * @return A factory used for create the underlying data provider.
     */
    protected IlluminaDataProviderFactory getFactory() {
        return factory;
    }

    /**
     * Must be called before doTileProcessing.  This is not passed in the ctor because often the
     * IlluminaDataProviderFactory is needed in order to construct the converter.
     *
     * @param converter Converts ClusterData to CLUSTER_OUTPUT_RECORD
     */
    protected void setConverter(final ClusterDataConverter<CLUSTER_OUTPUT_RECORD> converter) {
        this.converter = converter;
    }

    /**
     * Uses the firstTile and tileLimit parameters to set which tiles will be processed. The processor will
     * start with firstTile and continue to process tiles in order until it has processed at most tileLimit tiles.
     *
     * @param firstTile The tile to begin processing at.
     * @param tileLimit The maximum number of tiles to process.
     */
    protected void setTileLimits(final Integer firstTile, final Integer tileLimit) {
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
}
