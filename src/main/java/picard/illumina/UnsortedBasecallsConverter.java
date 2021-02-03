package picard.illumina;

import picard.PicardException;
import picard.illumina.parser.BaseIlluminaDataProvider;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;

import java.io.File;
import java.util.*;

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
    private int tilesProcessing = 0;

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
     * @param numThreads                   Controls number of threads.  If <= 0, the number of threads allocated is
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
        super(basecallsDir, barcodesDir, lane, readStructure, barcodeRecordWriterMap, demultiplex,
                numThreads, firstTile, tileLimit, bclQualityEvaluationStrategy,
                ignoreUnexpectedBarcodes, applyEamssFiltering, includeNonPfReads, 1);
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
    public void processTilesAndWritePerSampleOutputs(final Set<String> barcodes) {
        int MAX_TILES_IN_CACHE = 4;

        int tilesSubmitted = 0;
        while (tilesSubmitted < tiles.size()) {
            if (tilesProcessing < MAX_TILES_IN_CACHE) {
                int tile = tiles.get(tilesSubmitted);
                tileReadExecutor.submit(new TileReadProcessor(tile));
                tilesProcessing++;
                tilesSubmitted++;
            } else {
                try {
                    Thread.sleep(100);
                } catch (InterruptedException e) {
                    throw new PicardException("Tile processing thread interrupted: " + e.getMessage());
                }
            }
        }
        awaitTileProcessingCompletion();
    }

    /**
     * TileProcessor is a Runnable that process all records for a given tile. It uses the underlying
     * IlluminaDataProvider to iterate over cluster data for a specific tile. Records are added to a
     * cache as they are read.
     * <p>
     * After the tile processing is complete it notifies the CompletedWorkChecker that data is ready
     * for writing.
     */
    private class TileReadProcessor implements Runnable {
        private final int tileNum;

        TileReadProcessor(final int tileNum) {
            this.tileNum = tileNum;
        }

        @Override
        public void run() {
            Queue<ClusterData> queue = new ArrayDeque<>();
            final BaseIlluminaDataProvider dataProvider = factory.makeDataProvider(tileNum);

            while (dataProvider.hasNext()) {
                final ClusterData cluster = dataProvider.next();
                readProgressLogger.record(null, 0);
                if (includeNonPfReads || cluster.isPf()) {
                    queue.add(cluster);
                }
            }

            dataProvider.close();
            notifyWorkComplete(tileNum, Collections.singletonList(new RecordToWriterPump(queue)));
        }
    }

    private class RecordToWriterPump implements Runnable {
        private final Queue<ClusterData> clusterData;

        RecordToWriterPump(final Queue<ClusterData> clusterData) {
            this.clusterData = clusterData;
        }

        @Override
        public void run() {
            ClusterData cluster;
            while ((cluster = clusterData.poll()) != null) {
                if (includeNonPfReads || cluster.isPf()) {
                    barcodeRecordWriterMap.get(cluster.getMatchedBarcode()).write(converter.convertClusterToOutputRecord(cluster));
                    writeProgressLogger.record(null, 0);
                }
            }
            tilesProcessing--;
        }
    }
}
