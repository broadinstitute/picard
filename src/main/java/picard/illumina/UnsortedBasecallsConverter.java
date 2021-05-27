package picard.illumina;

import htsjdk.io.AsyncWriterPool;
import htsjdk.io.Writer;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.illumina.parser.BaseIlluminaDataProvider;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.IlluminaDataProviderFactory;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.Set;

/**
 * UnortedBasecallsConverter utilizes an underlying IlluminaDataProvider to convert parsed and decoded sequencing data
 * from standard Illumina formats to specific output records (FASTA records/SAM records). This data is processed
 * on a tile by tile basis.
 * <p>
 * The underlying IlluminaDataProvider applies several optional transformations that can include EAMSS filtering,
 * non-PF read filtering and quality score recoding using a BclQualityEvaluationStrategy.
 * <p>
 * The converter can also limit the scope of data that is converted from the data provider by setting the
 * tile to start on (firstTile) and the total number of tiles to process (tileLimit).
 * <p>
 * Additionally, BasecallsConverter can optionally demultiplex reads by outputting barcode specific reads to
 * their associated writers.
 */
public class UnsortedBasecallsConverter<CLUSTER_OUTPUT_RECORD> extends BasecallsConverter<CLUSTER_OUTPUT_RECORD> {
    private static final Log log = Log.getInstance(UnsortedBasecallsConverter.class);
    private final ProgressLogger progressLogger = new ProgressLogger(log, 1000000, "Processed");

    /**
     * Constructs a new BasecallsConverter object.
     *
     * @param basecallsDir                 Where to read basecalls from.
     * @param barcodesDir                  Where to read barcodes from (optional; use basecallsDir if not specified).
     * @param lanes                        What lane to process.
     * @param readStructure                How to interpret each cluster.
     * @param barcodeRecordWriterMap       Map from barcode to CLUSTER_OUTPUT_RECORD writer.  If demultiplex is false, must contain
     *                                     one writer stored with key=null.
     * @param demultiplex                  If true, output is split by barcode, otherwise all are written to the same output stream.
     *                                     available cores - numProcessors.
     * @param firstTile                    (For debugging) If non-null, start processing at this tile.
     * @param tileLimit                    (For debugging) If non-null, process no more than this many tiles.
     * @param bclQualityEvaluationStrategy The basecall quality evaluation strategy that is applyed to decoded base calls.
     * @param ignoreUnexpectedBarcodes     If true, will ignore reads whose called barcode is not found in barcodeRecordWriterMap.
     * @param applyEamssFiltering          If true, apply EAMSS filtering if parsing BCLs for bases and quality scores.
     * @param includeNonPfReads            If true, will include ALL reads (including those which do not have PF set).
     */
    protected UnsortedBasecallsConverter(
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
            final AsyncWriterPool writerPool
    ) {
        super(basecallsDir, barcodesDir, lanes, readStructure, barcodeRecordWriterMap, demultiplex,
                firstTile, tileLimit, bclQualityEvaluationStrategy,
                ignoreUnexpectedBarcodes, applyEamssFiltering, includeNonPfReads, writerPool);
    }

    /**
     * Set up tile processing and record writing threads for this converter.  This creates a tile reading thread
     * pool of size 4. The tile processing threads notify the completed work checking thread when they are
     * done processing a thread. The completed work checking thread will then dispatch the record writing for tiles
     * in order.
     *
     * @param barcodes The barcodes used for demultiplexing. When there is no demultiplexing done this should be a Set
     *                 containing a single null value.
     */
    @Override
    public void processTilesAndWritePerSampleOutputs(final Set<String> barcodes) throws IOException {
        for(IlluminaDataProviderFactory laneFactory : laneFactories) {
            for (Integer tileNum : tiles) {
                if (laneFactory.getAvailableTiles().contains(tileNum)) {
                    final BaseIlluminaDataProvider dataProvider = laneFactory.makeDataProvider(tileNum);

                    while (dataProvider.hasNext()) {
                        final ClusterData cluster = dataProvider.next();
                        if (includeNonPfReads || cluster.isPf()) {
                            barcodeRecordWriterMap.get(cluster.getMatchedBarcode()).write(converter.convertClusterToOutputRecord(cluster));
                            progressLogger.record(null, 0);
                        }
                    }
                    dataProvider.close();
                }
            }
        }
        closeWriters();
    }
}
