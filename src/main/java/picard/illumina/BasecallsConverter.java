package picard.illumina;

import htsjdk.io.AsyncWriterPool;
import htsjdk.io.Writer;
import htsjdk.samtools.util.IOUtil;
import picard.PicardException;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.IlluminaDataProviderFactory;
import picard.illumina.parser.IlluminaDataType;
import picard.illumina.parser.ReadData;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.util.ThreadPoolExecutorWithExceptions;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

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
    public static final Set<IlluminaDataType> DATA_TYPES_WITH_BARCODE = new HashSet<>(Arrays.asList(
            IlluminaDataType.BaseCalls,
            IlluminaDataType.QualityScores,
            IlluminaDataType.Position,
            IlluminaDataType.PF,
            IlluminaDataType.Barcodes));

    public static final Set<IlluminaDataType> DATA_TYPES_WITHOUT_BARCODE = new HashSet<>(Arrays.asList(
            IlluminaDataType.BaseCalls,
            IlluminaDataType.QualityScores,
            IlluminaDataType.Position,
            IlluminaDataType.PF));

    protected final IlluminaDataProviderFactory[] laneFactories;
    protected final boolean demultiplex;
    protected final boolean ignoreUnexpectedBarcodes;
    protected final Map<String, ? extends Writer<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap;
    protected final boolean includeNonPfReads;
    protected final AsyncWriterPool writerPool;
    protected ClusterDataConverter<CLUSTER_OUTPUT_RECORD> converter = null;
    protected List<Integer> tiles;
    protected BarcodeExtractor barcodeExtractor;

    /**
     * Constructs a new BasecallsConverter object.
     *
     * @param basecallsDir                 Where to read basecalls from.
     * @param barcodesDir                  Where to read barcodes from (optional; use basecallsDir if not specified).
     * @param lanes                        What lanes to process.
     * @param readStructure                How to interpret each cluster.
     * @param barcodeRecordWriterMap       Map from barcode to CLUSTER_OUTPUT_RECORD writer.  If demultiplex is false, must contain
     *                                     one writer stored with key=null.
     * @param demultiplex                  If true, output is split by barcode, otherwise all are written to the same output stream.
     * @param firstTile                    (For debugging) If non-null, start processing at this tile.
     * @param tileLimit                    (For debugging) If non-null, process no more than this many tiles.
     * @param bclQualityEvaluationStrategy The basecall quality evaluation strategy that is applyed to decoded base calls.
     * @param ignoreUnexpectedBarcodes     If true, will ignore reads whose called barcode is not found in barcodeRecordWriterMap.
     * @param applyEamssFiltering          If true, apply EAMSS filtering if parsing BCLs for bases and quality scores.
     * @param includeNonPfReads            If true, will include ALL reads (including those which do not have PF set).
     *                                     This option does nothing for instruments that output cbcls (Novaseqs)
     * @param barcodeExtractor             The `BarcodeExtractor` used to do inline barcode matching.
     */
    public BasecallsConverter(
            final File basecallsDir,
            final File barcodesDir,
            final int[] lanes,
            final ReadStructure readStructure,
            final Map<String, ? extends Writer<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap,
            final boolean demultiplex,
            final Integer firstTile,
            final Integer tileLimit,
            final BclQualityEvaluationStrategy bclQualityEvaluationStrategy,
            final boolean ignoreUnexpectedBarcodes,
            final boolean applyEamssFiltering,
            final boolean includeNonPfReads,
            final AsyncWriterPool writerPool,
            final BarcodeExtractor barcodeExtractor
    ) {
        this.barcodeRecordWriterMap = barcodeRecordWriterMap;
        this.ignoreUnexpectedBarcodes = ignoreUnexpectedBarcodes;
        this.demultiplex = demultiplex;
        this.barcodeExtractor = barcodeExtractor;
        this.writerPool = writerPool;
        this.laneFactories = new IlluminaDataProviderFactory[lanes.length];
        for(int i = 0; i < lanes.length; i++) {
            this.laneFactories[i] = new IlluminaDataProviderFactory(basecallsDir,
                    barcodesDir, lanes[i], readStructure, bclQualityEvaluationStrategy, getDataTypesFromReadStructure(readStructure, demultiplex, barcodesDir));
            this.laneFactories[i].setApplyEamssFiltering(applyEamssFiltering);
        }
        this.includeNonPfReads = includeNonPfReads;
        Set<Integer> allTiles = new TreeSet<>(TILE_NUMBER_COMPARATOR);
        for(IlluminaDataProviderFactory laneFactory: laneFactories) {
            allTiles.addAll(laneFactory.getAvailableTiles());
        }
        this.tiles = new ArrayList<>(allTiles);
        setTileLimits(firstTile, tileLimit);
    }

    /**
     * Abstract method for processing tiles of data and outputting records for each barcode.
     *
     * @param barcodes The barcodes used optionally for demultiplexing. Must contain at least a single null value if
     *                 no demultiplexing is being done.
     */
    public abstract void processTilesAndWritePerSampleOutputs(final Set<String> barcodes) throws IOException;

    /**
     * Closes all writers. If an AsycnWriterPool is used call close on that, otherwise iterate each writer and close it.
     *
     * @throws IOException throw if there is an error closing the writer.
     */
    public void closeWriters() throws IOException {
        if (writerPool != null) {
            writerPool.close();
        } else {
            for (Writer<CLUSTER_OUTPUT_RECORD> writer : barcodeRecordWriterMap.values()) {
                writer.close();
            }
        }
    }

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
    protected interface ConvertedClusterDataWriter<OUTPUT_RECORD> extends Writer<OUTPUT_RECORD> {
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
     * Because the tile number is followed by a colon, a tile number that is a prefix of another tile number
     * should sort after. (e.g. 10 sorts after 100). Tile numbers with the same number of digits are sorted numerically.
     */
    public static final Comparator<Integer> TILE_NUMBER_COMPARATOR = (integer1, integer2) -> {
        final String s1 = integer1.toString();
        final String s2 = integer2.toString();

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
     * @param barcodesDir   The barcodes dir that contains barcode files.
     * @return A data type array for each piece of data needed to satisfy the read structure.
     */
    protected static Set<IlluminaDataType> getDataTypesFromReadStructure(final ReadStructure readStructure,
                                                                         final boolean demultiplex,
                                                                         File barcodesDir) {
        if (!readStructure.hasSampleBarcode() || !demultiplex || barcodesDir == null) {
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
    protected IlluminaDataProviderFactory[] getLaneFactories() {
        return laneFactories;
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

    /**
     * If we are demultiplexing and a barcodeExtractor is defined then this method will perform on-the-fly demuxing.
     * Otherwise it will just return the pre-demuxed barcode from `ExtractIlluminaBarcodes`.
     * @param cluster The cluster data to demux
     * @param metrics The metrics object that will store the demux metrics.
     * @param noMatch A no-match metric object to store metrice for any read that doesn't demux
     * @param outputReadStructure The output `ReadStructure` for this cluster
     * @return The matched barcode or null if no barcode was matched.
     */
    protected String maybeDemultiplex(ClusterData cluster,
                                      Map<String, BarcodeMetric> metrics,
                                      BarcodeMetric noMatch,
                                      ReadStructure outputReadStructure) {
        String barcode = null;
        if (demultiplex) {
            // If a barcode extractor was provided for on-the-fly demux use it
            if (barcodeExtractor != null) {
                int[] barcodeIndices = outputReadStructure.sampleBarcodes.getIndices();
                byte[][] readSubsequences = new byte[barcodeIndices.length][];
                byte[][] qualityScores = new byte[barcodeIndices.length][];
                for (int i = 0; i < barcodeIndices.length; i++) {
                    ReadData barcodeRead = cluster.getRead(barcodeIndices[i]);
                    readSubsequences[i] = barcodeRead.getBases();
                    qualityScores[i] = barcodeRead.getQualities();
                }
                BarcodeExtractor.BarcodeMatch match = barcodeExtractor.findBestBarcode(readSubsequences,
                        qualityScores, true);

                BarcodeExtractor.updateMetrics(match, cluster.isPf(), metrics, noMatch);

                if (match.isMatched()) {
                  barcode = match.getBarcode();
                }
                cluster.setMatchedBarcode(barcode);
            } else {
                barcode = cluster.getMatchedBarcode();
            }
        }
        return barcode;
    }

    protected void interruptAndShutdownExecutors(ThreadPoolExecutorWithExceptions... executors) {
        final int tasksRunning = Arrays.stream(executors).mapToInt(test -> test.shutdownNow().size()).sum();
        final String errorMessages = Arrays.stream(executors).map(e -> {
            if (e.exception != null) {
                return e.exception.toString();
            } else {
                return "";
            }
        }).collect(Collectors.joining(","));
        throw new PicardException("Exceptions in tile processing. There were " + tasksRunning
                + " tasks still running or queued and they have been cancelled. Errors: " + errorMessages);
    }

    protected synchronized void updateMetrics(Map<String, BarcodeMetric> metrics, BarcodeMetric noMatch) {
        if(barcodeExtractor != null) {
            for (final String key : metrics.keySet()) {
                barcodeExtractor.getMetrics().get(key).merge(metrics.get(key));
            }
            barcodeExtractor.getNoMatchMetric().merge(noMatch);
        }
    }
}
