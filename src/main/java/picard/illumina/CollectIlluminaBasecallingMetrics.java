/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.illumina;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Illumina;
import picard.illumina.parser.BaseIlluminaDataProvider;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.IlluminaDataProviderFactory;
import picard.illumina.parser.IlluminaDataType;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.text.DecimalFormat;
import java.util.SortedMap;
import java.util.TreeMap;

/***
        - *  A Command line tool to collect Illumina Basecalling metrics for a sequencing run
        - *  Requires a Lane and an input file of Barcodes to expect.
        - *  Outputs metrics:
        - *    *  Mean Clusters Per Tile
        - *    *  Standard Deviation of Clusters Per Tile
        - *    *  Mean Pf Clusters Per Tile
        - *    *  Standard Deviation of Pf Clusters Per Tile
        - *    *  Mean Percentage of Pf Clusters Per Tile
        - *    *  Standard Deviation of Percentage of Pf Clusters Per Tile
        - */

@CommandLineProgramProperties(
        usage = CollectIlluminaBasecallingMetrics.USAGE_SUMMARY + CollectIlluminaBasecallingMetrics.USAGE_DETAILS,
        usageShort = CollectIlluminaBasecallingMetrics.USAGE_SUMMARY,
        programGroup = Illumina.class
)
public class CollectIlluminaBasecallingMetrics extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Collects Illumina Basecalling metrics for a sequencing run.  ";
    static final String USAGE_DETAILS = "<p>This tool will produce per-barcode and per-lane basecall metrics for each sequencing run.  " +
            "Mean values for each metric are determined using data from all of the tiles.  This tool requires the following data, LANE(#), " +
            "BASECALLS_DIR, READ_STRUCTURE, and an input file listing the sample barcodes.  " +
            "Program will provide metrics including: the total numbers of bases, reads, and clusters, as well as the fractions of each " +
            "bases, reads, and clusters that passed Illumina quality filters (PF) both per barcode and per lane.  " +
            "For additional information on Illumina's PF quality metric, please see the corresponding " +
            "<a href='https://www.broadinstitute.org/gatk/guide/article?id=6329'>GATK Dictionary entry</a>.</p> " +
            "<p>The input barcode_list.txt file is a file containing all of the sample and molecular barcodes and can be obtained from the " +
            "<a href='http://broadinstitute.github.io/picard/command-line-overview.html#ExtractIlluminaBarcodes'>ExtractIlluminaBarcodes</a> " +
            "tool.  </p>" +
            ""   +
            "Note: Metrics labeled as percentages are actually expressed as fractions!  " +
            "" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CollectIlluminaBasecallingMetrics \\<br />" +
            "      BASECALLS_DIR=/BaseCalls/ \\<br />" +
            "      LANE=001 \\<br />" +
            "      READ_STRUCTURE=25T8B25T \\<br />" +
            "      INPUT=barcode_list.txt " +
            "</pre>" +

            "<p>Please see the CollectIlluminaBasecallingMetrics " +
            "<a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#IlluminaBasecallingMetrics'>definitions</a> " +
            "for a complete description of the metrics produced by this tool.  </p>" +

            "<hr />"
    ;
    //Command Line Arguments

    @Option(doc="The Illumina basecalls output directory from which data are read", shortName="B")
    public File BASECALLS_DIR;

    @Option(doc = "The barcodes directory with _barcode.txt files (generated by ExtractIlluminaBarcodes). If not set, use BASECALLS_DIR. ", shortName = "BCD", optional = true)
    public File BARCODES_DIR;

    @Option(doc="The lane whose data will be read", shortName = StandardOptionDefinitions.LANE_SHORT_NAME)
    public Integer LANE;

    // TODO: No longer optional after old workflows are through
    @Option(doc="The file containing barcodes to expect from the run - barcodeData.#",shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, optional = true)
    public File INPUT;

    @Option(doc=ReadStructure.PARAMETER_DOC, shortName="RS")
    public String READ_STRUCTURE;

    @Option(doc="The file to which the collected metrics are written", shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional = true)
    public File OUTPUT;

    private int barcodeLength = 0;
    private String unmatched_barcode;
    private final SortedMap<String, IlluminaMetricCounts> barcodeToMetricCounts;

    private static final String BARCODE_NAME_COLUMN = "barcode_name";
    private static final String BARCODE_SEQUENCE_COLUMN_NAME_STUB = "barcode_sequence_";

    public CollectIlluminaBasecallingMetrics() {
        this.barcodeToMetricCounts = new TreeMap<String, IlluminaMetricCounts>();
    }

    @Override
    protected int doWork() {
        // File and Directory Validation
        IOUtil.assertDirectoryIsReadable(BASECALLS_DIR);
        if (OUTPUT == null) OUTPUT = new File(BASECALLS_DIR, String.format("LANE%s_basecalling_metrics", LANE));
        IOUtil.assertFileIsWritable(OUTPUT);

        final IlluminaDataProviderFactory factory;
        final ReadStructure readStructure = new ReadStructure(READ_STRUCTURE);
        final BclQualityEvaluationStrategy bclQualityEvaluationStrategy = new BclQualityEvaluationStrategy(BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY);

        if (INPUT == null) {
            // TODO: Legacy support. Remove when INPUT is required, after all old workflows are through
            factory = new IlluminaDataProviderFactory(BASECALLS_DIR, LANE, readStructure, bclQualityEvaluationStrategy,
                    IlluminaDataType.PF, IlluminaDataType.Position);
        } else {
            // Grab expected barcode data from barcodeData.<LANE>
            IOUtil.assertFileIsReadable(INPUT);
            final TabbedTextFileWithHeaderParser barcodesParser = new TabbedTextFileWithHeaderParser(INPUT);
            for (final TabbedTextFileWithHeaderParser.Row row : barcodesParser) {
                final String barcodeName = row.getField(BARCODE_NAME_COLUMN);
                final StringBuilder barcode = new StringBuilder();
                for (int i = 1; i <= readStructure.sampleBarcodes.length(); i++) {
                    barcode.append(row.getField(BARCODE_SEQUENCE_COLUMN_NAME_STUB + i));
                    if (barcodeLength == 0) barcodeLength = barcode.length();
                }

                // Only add the barcode to the hash if it has sequences. For libraries
                // that don't have barcodes this won't be set in the file.
                if (barcode.length() > 0) {
                    barcodeToMetricCounts.put(barcode.toString(), new IlluminaMetricCounts(barcode.toString(), barcodeName, LANE));
                }
            }

            factory = barcodeToMetricCounts.isEmpty()
                    ? new IlluminaDataProviderFactory(
                        BASECALLS_DIR,
                        BARCODES_DIR,
                        LANE,
                        readStructure,
                        bclQualityEvaluationStrategy,
                        IlluminaDataType.PF,
                        IlluminaDataType.Position)
                    : new IlluminaDataProviderFactory(
                        BASECALLS_DIR,
                        BARCODES_DIR,
                        LANE,
                        readStructure,
                        bclQualityEvaluationStrategy,
                        IlluminaDataType.PF,
                        IlluminaDataType.Position,
                        IlluminaDataType.Barcodes);
        }

        unmatched_barcode = StringUtil.repeatCharNTimes('N', barcodeLength);

        //Initialize data provider, iterate over clusters, and collect statistics
        final BaseIlluminaDataProvider provider = factory.makeDataProvider();

        while (provider.hasNext()) {
            final ClusterData cluster = provider.next();
            addCluster(cluster);
        }

        onComplete();
        return 0;
    }

    /***
     * Process new cluster of Illumina data - increment a running counter of data
     */
    private void addCluster(final ClusterData cluster) {
        //compute hash of Barcode and Lane for key
        String barcode = cluster.getMatchedBarcode();
        if (barcode == null) barcode = unmatched_barcode;

        //increment counts
        IlluminaMetricCounts counters =  barcodeToMetricCounts.get(barcode);
        if (counters == null) {
             counters = new IlluminaMetricCounts(barcode,null,LANE);
             barcodeToMetricCounts.put(barcode, counters);
        }
        final int tileNumber = cluster.getTile();
        counters.incrementClusterCount(tileNumber,cluster.isPf());
    }

    /**
     * Handles completion of metric collection. Metrics are computed from counts of data and written out to a file.
     */
    private void onComplete() {
        try {
            final MetricsFile<IlluminaBasecallingMetrics, Comparable<?>> file = getMetricsFile();
            final IlluminaMetricCounts allLaneCounts = new IlluminaMetricCounts(null, null, LANE);
            for (final String s : barcodeToMetricCounts.keySet()) {
                final IlluminaMetricCounts counts = barcodeToMetricCounts.get(s);
                counts.addMetricsToFile(file);
                allLaneCounts.addIlluminaMetricCounts(counts);
            }
            if (! barcodeToMetricCounts.keySet().contains("")) allLaneCounts.addMetricsToFile(file);  // detect non-indexed case
            file.write(OUTPUT);
        } catch (final Exception ex) {
            throw new PicardException("Error writing output file " + OUTPUT.getPath(), ex);
        }
    }

    public static void main(final String[] argv) {
        new CollectIlluminaBasecallingMetrics().instanceMainWithExit(argv);
    }

    /***
     * This class manages counts of Illumina Basecalling data on a Per Barcode Per Lane basis.  Cluster and PFCluster
     * counts are stored per tile number.
     */
    private class IlluminaMetricCounts {
        /*** Stores counts of clusters found for a specific Barcode-Lane combination across all tiles.  Key = Tile Number, Value = count of clusters***/
        private final Histogram<Integer> tileToClusterHistogram;
        /*** Stores counts of pf clusters found for a specific Barcode-Lane combination across all tiles.  Key = Tile Number, Value = count of clusters***/
        private final Histogram<Integer> tileToPfClusterHistogram;
        final IlluminaBasecallingMetrics metrics;

        public IlluminaMetricCounts(final String barcode, final String barcodeName, final Integer laneNumber) {
            this.tileToClusterHistogram = new Histogram<Integer>();
            this.tileToPfClusterHistogram = new Histogram<Integer>();
            this.metrics = new IlluminaBasecallingMetrics();
            this.metrics.MOLECULAR_BARCODE_SEQUENCE_1 = barcode;
            this.metrics.MOLECULAR_BARCODE_NAME = barcodeName;
            this.metrics.LANE = Integer.toString(laneNumber);
        }

        /*  Increments cluster count by 1 for a given tile number */
        public void incrementClusterCount(final int tileNumber, final boolean isPf) {
            incrementClusterCount(tileNumber,1d, isPf);
        }

        /*  Increments cluster count by an amount for a given tile number */
        public void incrementClusterCount(final int tileNumber, final double incrementAmount, final boolean isPf) {
            incrementClusterCount(tileNumber, incrementAmount, (isPf ? 1d : 0d));
        }

        /*  Increments cluster count by an amount for a given tile number */
        public void incrementClusterCount(final Integer tileNumber, final double incrementAmount, final double pfIncrementAmount) {
            tileToClusterHistogram.increment(tileNumber, incrementAmount);
            tileToPfClusterHistogram.increment(tileNumber, pfIncrementAmount);
        }

        /* Handles calculating final metrics and updating the metric object */
        private void onComplete() {
            final double meanClustersPerTile =  tileToClusterHistogram.getMeanBinSize();
            metrics.MEAN_CLUSTERS_PER_TILE = Math.round(meanClustersPerTile);
            metrics.SD_CLUSTERS_PER_TILE = Math.round(tileToClusterHistogram.getStandardDeviationBinSize(meanClustersPerTile));

            final double meanPfClustersPerTile =  tileToPfClusterHistogram.getMeanBinSize();
            metrics.MEAN_PF_CLUSTERS_PER_TILE = Math.round(meanPfClustersPerTile);
            metrics.SD_PF_CLUSTERS_PER_TILE = Math.round(tileToPfClusterHistogram.getStandardDeviationBinSize(meanPfClustersPerTile));

            final DecimalFormat decFormat = new DecimalFormat("#.##");
            final Histogram<Integer> laneToPctPfClusterHistogram = tileToPfClusterHistogram.divideByHistogram(tileToClusterHistogram);
            final double meanPctPfClustersPerTile = laneToPctPfClusterHistogram.getMeanBinSize();
            metrics.MEAN_PCT_PF_CLUSTERS_PER_TILE = (Double.isNaN(meanPctPfClustersPerTile) ?  0 : Double.valueOf(decFormat.format(meanPctPfClustersPerTile * 100)));
            metrics.SD_PCT_PF_CLUSTERS_PER_TILE = Double.valueOf(decFormat.format(laneToPctPfClusterHistogram.getStandardDeviationBinSize(meanPctPfClustersPerTile) * 100));

            metrics.TOTAL_CLUSTERS = (long) this.tileToClusterHistogram.getSumOfValues();
            metrics.PF_CLUSTERS    = (long) this.tileToPfClusterHistogram.getSumOfValues();
            
            final ReadStructure readStructure = new ReadStructure(READ_STRUCTURE);
            int templateBaseCountPerCluster = 0;
            for (int i = 0; i < readStructure.templates.length(); i++) templateBaseCountPerCluster += readStructure.templates.get(i).length;
            metrics.TOTAL_READS = metrics.TOTAL_CLUSTERS * readStructure.templates.length();
            metrics.PF_READS = metrics.PF_CLUSTERS * readStructure.templates.length();
            metrics.TOTAL_BASES = metrics.TOTAL_CLUSTERS * templateBaseCountPerCluster;
            metrics.PF_BASES = metrics.PF_CLUSTERS * templateBaseCountPerCluster;
            
        }

        /* Computes final metric based on data counts and writes to output metric file */
        public void addMetricsToFile(final MetricsFile<IlluminaBasecallingMetrics, Comparable<?>> file) {
            onComplete();
            file.addMetric(metrics);
        }

        /*  Merges data from another IlluminaMetricCount object into current one.*/
        public void addIlluminaMetricCounts(final IlluminaMetricCounts counts) {
            this.tileToClusterHistogram.addHistogram(counts.tileToClusterHistogram);
            this.tileToPfClusterHistogram.addHistogram(counts.tileToPfClusterHistogram);
        }
    }
}
