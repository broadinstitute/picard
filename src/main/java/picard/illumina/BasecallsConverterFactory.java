package picard.illumina;

import htsjdk.samtools.util.SortingCollection;
import picard.illumina.parser.IlluminaFileUtil;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;

import java.io.File;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

public class BasecallsConverterFactory<CLUSTER_OUTPUT_RECORD> {
    public BasecallsConverter<CLUSTER_OUTPUT_RECORD> getConverter(
            final File basecallsDir, File barcodesDir, final int lane,
            final ReadStructure readStructure,
            final Map<String, ? extends BasecallsConverter.ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap,
            final boolean demultiplex,
            final int maxReadsInRamPerTile,
            final List<File> tmpDirs,
            final int numProcessors,
            final Integer firstTile,
            final Integer tileLimit,
            final Comparator<CLUSTER_OUTPUT_RECORD> outputRecordComparator,
            final SortingCollection.Codec<CLUSTER_OUTPUT_RECORD> codecPrototype,
            final Class<CLUSTER_OUTPUT_RECORD> outputRecordClass,
            final BclQualityEvaluationStrategy bclQualityEvaluationStrategy,
            final boolean ignoreUnexpectedBarcodes,
            final boolean applyEamssFiltering,
            final boolean includeNonPfReads,
            final boolean forceGc
    ) {
        if (IlluminaFileUtil.hasCbcls(basecallsDir, lane)) {
            if (barcodesDir == null) barcodesDir = basecallsDir;
            return new NewIlluminaBasecallsConverter<>(basecallsDir, barcodesDir, lane, readStructure,
                    barcodeRecordWriterMap, demultiplex, maxReadsInRamPerTile,
                    tmpDirs, numProcessors,
                    firstTile, tileLimit, outputRecordComparator,
                    codecPrototype,
                    outputRecordClass, bclQualityEvaluationStrategy, ignoreUnexpectedBarcodes);
        } else {
            return new IlluminaBasecallsConverter<>(basecallsDir, barcodesDir, lane, readStructure,
                    barcodeRecordWriterMap, demultiplex, maxReadsInRamPerTile, tmpDirs, numProcessors,
                    forceGc, firstTile, tileLimit, outputRecordComparator,
                    codecPrototype,
                    outputRecordClass, bclQualityEvaluationStrategy,
                    applyEamssFiltering, includeNonPfReads, ignoreUnexpectedBarcodes);
        }
    }
}
