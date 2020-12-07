package picard.illumina;

import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.util.SortingCollection;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;

import java.io.File;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

public class BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> {

    File basecallsDir;
    File barcodesDir;
    int lane;
    ReadStructure readStructure;

    Comparator<CLUSTER_OUTPUT_RECORD> outputRecordComparator;
    SortingCollection.Codec<CLUSTER_OUTPUT_RECORD> codecPrototype;
    Class<CLUSTER_OUTPUT_RECORD> outputRecordClass;

    Map<String, ? extends BasecallsConverter.ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap;
    boolean demultiplex = false;
    int maxReadsInRamPerTile = 1200000;
    List<File> tmpDirs = Collections.singletonList(new File(System.getProperty("java.io.tmpdir")));
    int numProcessors = Runtime.getRuntime().availableProcessors();
    Integer firstTile = null;
    Integer tileLimit = null;
    BclQualityEvaluationStrategy bclQualityEvaluationStrategy =
            new BclQualityEvaluationStrategy(BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY);;
    boolean ignoreUnexpectedBarcodes = false;
    boolean applyEamssFiltering = false;
    boolean includeNonPfReads = false;

    public BasecallsConverterBuilder(final File basecallsDir, final Integer lane, final ReadStructure readStructure,
                                     Map<String, ? extends BasecallsConverter.ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap) {
        this.basecallsDir = basecallsDir;
        this.lane = lane;
        this.readStructure = readStructure;
        this.barcodeRecordWriterMap = barcodeRecordWriterMap;
    }

    public BasecallsConverter<CLUSTER_OUTPUT_RECORD> build() {
        return new BasecallsConverter<>(basecallsDir, barcodesDir, lane, readStructure,
                barcodeRecordWriterMap, demultiplex, maxReadsInRamPerTile,
                tmpDirs, numProcessors,
                firstTile, tileLimit, outputRecordComparator,
                codecPrototype,
                outputRecordClass, bclQualityEvaluationStrategy, ignoreUnexpectedBarcodes, applyEamssFiltering, includeNonPfReads);
    }

    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> withIgnoreUnexpectedBarcodes(boolean ignoreUnexpectedBarcodes) {
        this.ignoreUnexpectedBarcodes = ignoreUnexpectedBarcodes;
        return this;
    }

    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> withApplyEamssFiltering(boolean applyEamssFiltering) {
        this.applyEamssFiltering = applyEamssFiltering;
        return this;
    }

    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> withIncludeNonPfReads(boolean includeNonPfReads) {
        this.includeNonPfReads = includeNonPfReads;
        return this;
    }

    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> bclQualityEvaluationStrategy(BclQualityEvaluationStrategy bclQualityEvaluationStrategy) {
        this.bclQualityEvaluationStrategy = bclQualityEvaluationStrategy;
        return this;
    }

    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> withSorting(
            Comparator<CLUSTER_OUTPUT_RECORD> outputRecordComparator,
            SortingCollection.Codec<CLUSTER_OUTPUT_RECORD> codecPrototype,
            Class<CLUSTER_OUTPUT_RECORD> outputRecordClass,
            Integer maxReadsInRamPerTile,
            List<File> tmpDirs) {
        this.outputRecordComparator = outputRecordComparator;
        this.codecPrototype = codecPrototype;
        this.outputRecordClass = outputRecordClass;
        this.maxReadsInRamPerTile = maxReadsInRamPerTile;
        this.tmpDirs = tmpDirs;
        return this;
    }

    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> tileLimit(Integer tileLimit) {
        this.tileLimit = tileLimit;
        return this;
    }

    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> firstTile(Integer firstTile) {
        this.firstTile = firstTile;
        return this;
    }

    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> numProcessors(Integer numProcessors) {
        this.numProcessors = numProcessors;
        return this;
    }

    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> withDemultiplex(boolean demultiplex) {
        this.demultiplex = demultiplex;
        return this;
    }

    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> barcodeRecordWriterMap(
            Map<String, ? extends BasecallsConverter.ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap) {
        this.barcodeRecordWriterMap = barcodeRecordWriterMap;
        return this;
    }

    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> barcodesDir(File barcodesDir) {
        this.barcodesDir = (barcodesDir == null) ? basecallsDir : barcodesDir;
        return this;
    }
}
