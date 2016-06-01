package picard.sam.markduplicates;

/**
 *
 * The purpose of this class is to show that MarkDuplicates gives the same results when run on files that do not have a
 * molecular barcode tag, even if the code is trying to use the molecular barcode
 *
 */

public class MarkDuplicateWithMissingReadTwoBarcodeTest extends MarkDuplicateWithMissingBarcodeTest {

    @Override
    protected String getArgumentName() {
        return "READ_TWO_BARCODE_TAG";
    }

    @Override
    protected String getTagValue() {
        return "RX";
    }
}
