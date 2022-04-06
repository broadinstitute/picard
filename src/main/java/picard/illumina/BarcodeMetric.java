package picard.illumina;

import htsjdk.samtools.metrics.MetricBase;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.util.IlluminaUtil;
import picard.util.help.HelpConstants;

/**
 * Metrics produced by the ExtractIlluminaBarcodes program that is used to parse data in
 * the basecalls directory and determine to which barcode each read should be assigned.
 */
@DocumentedFeature(groupName = HelpConstants.DOC_CAT_METRICS, summary = HelpConstants.DOC_CAT_METRICS_SUMMARY)
public class BarcodeMetric extends MetricBase {
    /**
     * The barcode (from the set of expected barcodes) for which the following metrics apply.
     * Note that the "symbolic" barcode of NNNNNN is used to report metrics for all reads that
     * do not match a barcode.
     */
    public String BARCODE;

    public String BARCODE_WITHOUT_DELIMITER;
    /**
     * The barcode name.
     */
    public String BARCODE_NAME = "";
    /**
     * The name of the library
     */
    public String LIBRARY_NAME = "";
    /**
     * The total number of reads matching the barcode.
     */
    public long READS = 0;
    /**
     * The number of PF reads matching this barcode (always less than or equal to READS).
     */
    public long PF_READS = 0;
    /**
     * The number of all reads matching this barcode that matched with 0 errors or no-calls.
     */
    public long PERFECT_MATCHES = 0;
    /**
     * The number of PF reads matching this barcode that matched with 0 errors or no-calls.
     */
    public long PF_PERFECT_MATCHES = 0;
    /**
     * The number of all reads matching this barcode that matched with 1 error or no-call.
     */
    public long ONE_MISMATCH_MATCHES = 0;
    /**
     * The number of PF reads matching this barcode that matched with 1 error or no-call.
     */
    public long PF_ONE_MISMATCH_MATCHES = 0;
    /**
     * The fraction of all reads in the lane that matched to this barcode.
     */
    public double PCT_MATCHES = 0d;
    /**
     * The rate of all reads matching this barcode to all reads matching the most prevalent barcode. For the
     * most prevalent barcode this will be 1, for all others it will be less than 1 (except for the possible
     * exception of when there are more orphan reads than for any other barcode, in which case the value
     * may be arbitrarily large).  One over the lowest number in this column gives you the fold-difference
     * in representation between barcodes.
     */
    public double RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT = 0d;
    /**
     * The fraction of PF reads in the lane that matched to this barcode.
     */
    public double PF_PCT_MATCHES = 0d;

    /**
     * The rate of PF reads matching this barcode to PF reads matching the most prevalent barcode. For the
     * most prevalent barcode this will be 1, for all others it will be less than 1 (except for the possible
     * exception of when there are more orphan reads than for any other barcode, in which case the value
     * may be arbitrarily large).  One over the lowest number in this column gives you the fold-difference
     * in representation of PF reads between barcodes.
     */
    public double PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT = 0d;

    /**
     * The "normalized" matches to each barcode. This is calculated as the number of pf reads matching
     * this barcode over the sum of all pf reads matching any barcode (excluding orphans). If all barcodes
     * are represented equally this will be 1.
     */
    public double PF_NORMALIZED_MATCHES;

    protected byte[][] barcodeBytes;

    public BarcodeMetric(final String barcodeName, final String libraryName,
                         final String barcodeDisplay, final String[] barcodeSeqs) {

        this.BARCODE = barcodeDisplay;
        this.BARCODE_WITHOUT_DELIMITER = barcodeDisplay.replaceAll(IlluminaUtil.BARCODE_DELIMITER, "");
        this.BARCODE_NAME = barcodeName;
        this.LIBRARY_NAME = libraryName;
        this.barcodeBytes = new byte[barcodeSeqs.length][];
        for (int i = 0; i < barcodeSeqs.length; i++) {
            barcodeBytes[i] = htsjdk.samtools.util.StringUtil.stringToBytes(barcodeSeqs[i]);
        }
    }

    /**
     * This ctor is necessary for when reading metrics from file
     */
    public BarcodeMetric() {
        barcodeBytes = null;
    }

    /**
     * Creates a copy of metric initialized with only non-accumulated and non-calculated values set
     */
    public BarcodeMetric copy() {
        final BarcodeMetric result = new BarcodeMetric();
        result.BARCODE = this.BARCODE;
        result.BARCODE_WITHOUT_DELIMITER = this.BARCODE_WITHOUT_DELIMITER;
        result.BARCODE_NAME = this.BARCODE_NAME;
        result.LIBRARY_NAME = this.LIBRARY_NAME;
        result.barcodeBytes = this.barcodeBytes;
        return result;
    }

    /**
     * Adds the non-calculated
     */
    public void merge(final BarcodeMetric metric) {
        this.READS += metric.READS;
        this.PF_READS += metric.PF_READS;
        this.PERFECT_MATCHES += metric.PERFECT_MATCHES;
        this.PF_PERFECT_MATCHES += metric.PF_PERFECT_MATCHES;
        this.ONE_MISMATCH_MATCHES += metric.ONE_MISMATCH_MATCHES;
        this.PF_ONE_MISMATCH_MATCHES += metric.PF_ONE_MISMATCH_MATCHES;
    }

}
