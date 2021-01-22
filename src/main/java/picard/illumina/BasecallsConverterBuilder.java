package picard.illumina;

import htsjdk.samtools.SAMFileWriterImpl;
import htsjdk.samtools.util.SortingCollection;
import picard.PicardException;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;

import java.io.File;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

/**
 * BasecallsConverterBuilder creates and configures BasecallsConverter objects. It builds and returns a
 * BasecallsConverter that includes record sorting.
 *
 * @param <CLUSTER_OUTPUT_RECORD> The type of records that the built converter will output.
 */
public class BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> {
    private final File basecallsDir;
    private final int lane;
    private final ReadStructure readStructure;
    private final Map<String, ? extends BasecallsConverter.ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap;
    private Comparator<CLUSTER_OUTPUT_RECORD> outputRecordComparator;
    private SortingCollection.Codec<CLUSTER_OUTPUT_RECORD> codecPrototype;
    private Class<CLUSTER_OUTPUT_RECORD> outputRecordClass;
    private int maxReadsInRamPerThread = SAMFileWriterImpl.getDefaultMaxRecordsInRam();
    private List<File> tmpDirs = Collections.singletonList(new File(System.getProperty("java.io.tmpdir")));
    private File barcodesDir;
    private boolean demultiplex = false;
    private int numThreads = Runtime.getRuntime().availableProcessors();
    private Integer firstTile = null;
    private Integer tileLimit = null;
    private BclQualityEvaluationStrategy bclQualityEvaluationStrategy = new BclQualityEvaluationStrategy(BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY);
    private boolean ignoreUnexpectedBarcodes = false;
    private boolean applyEamssFiltering = false;
    private boolean includeNonPfReads = false;

    /**
     * Constructs a new builder used for creating BasecallsConverter objects.
     *
     * @param basecallsDir           Where to read basecalls from.
     * @param lane                   What lane to process.
     * @param readStructure          How to interpret each cluster.
     * @param barcodeRecordWriterMap Map from barcode to CLUSTER_OUTPUT_RECORD writer.  If demultiplex is false, must contain
     *                               one writer stored with key=null.
     */
    public BasecallsConverterBuilder(final File basecallsDir, final Integer lane, final ReadStructure readStructure,
                                     Map<String, ? extends BasecallsConverter.ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap) {
        this.basecallsDir = basecallsDir;
        this.lane = lane;
        this.readStructure = readStructure;
        this.barcodeRecordWriterMap = barcodeRecordWriterMap;
    }

    /**
     * Builds a new sorting basecalls converter.
     *
     * @param outputRecordComparator For sorting output records within a single tile.
     * @param codecPrototype         For spilling output records to disk.
     * @param outputRecordClass      Class needed to create SortingCollections.
     * @param tmpDirs                For SortingCollection spilling.
     * @return A basecalls converter that will output sorted records.
     */
    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> withSorting(Comparator<CLUSTER_OUTPUT_RECORD> outputRecordComparator,
                                                                        SortingCollection.Codec<CLUSTER_OUTPUT_RECORD> codecPrototype,
                                                                        Class<CLUSTER_OUTPUT_RECORD> outputRecordClass,
                                                                        List<File> tmpDirs) {
        if(outputRecordComparator == null || codecPrototype == null || outputRecordClass == null){
            throw new PicardException("outputRecordComparator, codecPrototype and outputRecordClasse can not be null.");
        }
        this.outputRecordComparator = outputRecordComparator;
        this.codecPrototype = codecPrototype;
        this.outputRecordClass = outputRecordClass;
        this.tmpDirs = tmpDirs;
        return this;
    }

    /**
     * Builds a basecalls converter
     *
     * @return A basecalls converter that will output records according to the parameters set.
     */
    public BasecallsConverter<CLUSTER_OUTPUT_RECORD> build() {

        if (outputRecordComparator != null && codecPrototype != null && outputRecordClass != null && tmpDirs != null) {
            return new SortedBasecallsConverter<>(basecallsDir, barcodesDir, lane, readStructure,
                    barcodeRecordWriterMap, demultiplex, maxReadsInRamPerThread,
                    tmpDirs, numThreads,
                    firstTile, tileLimit, outputRecordComparator,
                    codecPrototype,
                    outputRecordClass, bclQualityEvaluationStrategy, ignoreUnexpectedBarcodes, applyEamssFiltering, includeNonPfReads);
        } else {
            return new UnsortedBasecallsConverter<>(basecallsDir, barcodesDir, lane, readStructure,
                    barcodeRecordWriterMap, demultiplex, numThreads, firstTile, tileLimit,
                    bclQualityEvaluationStrategy, ignoreUnexpectedBarcodes, applyEamssFiltering, includeNonPfReads);
        }
    }

    /**
     * Configures whether or not the converter will ignore unexpected barcodes or throw an exception if one is found.
     *
     * @param ignoreUnexpectedBarcodes If true, will ignore reads whose called barcode is not found in barcodeRecordWriterMap.
     * @return A builder that will create a converter with the ignoreUnexpectedBarcodes boolean set.
     */
    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> withIgnoreUnexpectedBarcodes(boolean ignoreUnexpectedBarcodes) {
        this.ignoreUnexpectedBarcodes = ignoreUnexpectedBarcodes;
        return this;
    }

    /**
     * Configures whether or not the converter will apply EAMSS filtering.
     *
     * @param applyEamssFiltering If true, apply EAMSS filtering if parsing BCLs for bases and quality scores.
     * @return A builder that will create a converter with the applyEamssFiltering boolean set.
     */
    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> withApplyEamssFiltering(boolean applyEamssFiltering) {
        this.applyEamssFiltering = applyEamssFiltering;
        return this;
    }

    /**
     * Configures whether or not the converter will ignore non-PF reads.
     *
     * @param includeNonPfReads If true, will include ALL reads (including those which do not have PF set).
     *                          This option does nothing for instruments that output cbcls (Novaseqs)
     * @return A builder that will create a converter with the includeNonPfReads boolean set.
     */
    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> withIncludeNonPfReads(boolean includeNonPfReads) {
        this.includeNonPfReads = includeNonPfReads;
        return this;
    }

    /**
     * Configures the bcl quality evaluation strategy that the converter will apply.
     *
     * @param bclQualityEvaluationStrategy The mechanism for revising and evaluating qualities read from a BCL file
     * @return A builder that will create a converter with the bclQualityEvaluationStrategy set.
     */
    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> withBclQualityEvaluationStrategy(BclQualityEvaluationStrategy bclQualityEvaluationStrategy) {
        this.bclQualityEvaluationStrategy = bclQualityEvaluationStrategy;
        return this;
    }

    /**
     * Configures the total number of tiles that the converter will process.
     *
     * @param tileLimit If non-null, process no more than this many tiles.
     * @return A builder that will create a converter with tileLimit set.
     */
    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> tileLimit(Integer tileLimit) {
        this.tileLimit = tileLimit;
        return this;
    }

    /**
     * Configures the first tile that the converter will begin processing at.
     *
     * @param firstTile If non-null, start processing at this tile.
     * @return A builder that will create a converter with firstTile set.
     */
    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> firstTile(Integer firstTile) {
        this.firstTile = firstTile;
        return this;
    }

    /**
     * Configures how many processors this converter will use.
     *
     * @param numProcessors Controls number of threads.  If <= 0, the number of threads allocated is
     *                      available cores - numProcessors.
     * @return A builder that will create a converter with numProcessors set.
     */
    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> numProcessors(int numProcessors) {
    	if (numProcessors == 0) {
            this.numThreads = Runtime.getRuntime().availableProcessors();
        } else if (numProcessors < 0) {
            this.numThreads = Runtime.getRuntime().availableProcessors() + numProcessors;
        } else {
           this.numThreads = numProcessors;
        }

        return this;
    }

    /**
     * Configures whether or not the converter will demultiplex reads by barcode.
     *
     * @param demultiplex If true, output is split by barcode, otherwise all are written to the same output stream.
     * @return A builder that will create a converter with demultiplex set.
     */
    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> withDemultiplex(boolean demultiplex) {
        this.demultiplex = demultiplex;
        return this;
    }

    /**
     * Configure the director used to find barcode files created by ExtractIlluminaBarcodes. These files are
     * used to demultiplex reads.
     *
     * @param barcodesDir Where to read barcodes from (optional; use basecallsDir if not specified).
     * @return A builder that will create a converter with barcodesDir set.
     */
    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> barcodesDir(File barcodesDir) {
        this.barcodesDir = (barcodesDir == null) ? basecallsDir : barcodesDir;
        return this;
    }

    /**
     * Specify the max number of records in RAM. This is divided by the number of tile processing threads and
     * sets the max number of records in RAM per tile.
     *
     * @param maxReadsInRam The maximum number of records in RAM to store for each tile before spilling to disk.
     * @return A builder that will create a converter with the maximum records in RAM set to `maxReadsInRam/numThreads`
     */
    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> withMaxRecordsInRam(int maxReadsInRam) {
        this.maxReadsInRamPerThread = Math.max(1, maxReadsInRam / this.numThreads);
        return this;
    }
}
