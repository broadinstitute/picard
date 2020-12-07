package picard.illumina;

import htsjdk.samtools.util.SortingCollection;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;

import java.io.File;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

public class BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> {
    File basecallsDir;
    File barcodesDir;
    int lane;
    ReadStructure readStructure;
    Map<String, ? extends BasecallsConverter.ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap;
    boolean demultiplex;
    int maxReadsInRamPerTile;
    List<File> tmpDirs;
    int numProcessors;
    Integer firstTile;
    Integer tileLimit;
    Comparator<CLUSTER_OUTPUT_RECORD> outputRecordComparator;
    SortingCollection.Codec<CLUSTER_OUTPUT_RECORD> codecPrototype;
    Class<CLUSTER_OUTPUT_RECORD> outputRecordClass;
    BclQualityEvaluationStrategy bclQualityEvaluationStrategy;
    boolean ignoreUnexpectedBarcodes = false;
    boolean applyEamssFiltering = false;
    boolean includeNonPfReads = false;

    public BasecallsConverterBuilder(final File basecallsDir) {
        this.basecallsDir = basecallsDir;
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

    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> outputRecordComparator(Comparator<CLUSTER_OUTPUT_RECORD> outputRecordComparator) {
        this.outputRecordComparator = outputRecordComparator;
        return this;
    }

    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> codecPrototype(SortingCollection.Codec<CLUSTER_OUTPUT_RECORD> codecPrototype) {
        this.codecPrototype = codecPrototype;
        return this;
    }

    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> outputRecordClass(Class<CLUSTER_OUTPUT_RECORD> outputRecordClass) {
        this.outputRecordClass = outputRecordClass;
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

    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> tmpDirs(List<File> tmpDirs) {
        this.tmpDirs = tmpDirs;
        return this;
    }

    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> maxReadsInRamPerTile(Integer maxReadsInRamPerTile) {
        this.maxReadsInRamPerTile = maxReadsInRamPerTile;
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

    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> readStructure(ReadStructure readStructure) {
        this.readStructure = readStructure;
        return this;
    }

    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> lane(Integer lane) {
        this.lane = lane;
        return this;
    }

    public BasecallsConverterBuilder<CLUSTER_OUTPUT_RECORD> barcodesDir(File barcodesDir) {
        this.barcodesDir = (barcodesDir == null) ? basecallsDir : barcodesDir;
        return this;
    }


}
