package picard.sam.markduplicates;

/**
 * The purpose of this class is to show that MarkDuplicates gives the same results when run on files that do not have a
 * molecular barcode tag, even if the code is trying to use the molecular barcode
 */

abstract public class MarkDuplicateWithMissingBarcodeTest extends MarkDuplicatesTest {

    protected AbstractMarkDuplicatesCommandLineProgramTester getTester() {
        return new MarkDuplicatesWithMissingBarcodesTester();
    }

    abstract protected String getArgumentName();

    abstract protected String getTagValue();

    private class MarkDuplicatesWithMissingBarcodesTester extends MarkDuplicatesTester {
        @Override
        public void runTest() {
            boolean hasRX = false;
            boolean isDuplex = false;
            for (final String argument : this.getArgs()) {
                if (argument.startsWith(getArgumentName())) {
                    hasRX = true;
                    break;
                }
                if (argument.startsWith("DUPLEX_UMI")) {
                    isDuplex = true;
                    break;
                }
            }

            // Tests that specify the RX tag, or the DUPLEX_UMI argument should not be used with missing barcodes.
            if (!hasRX && !isDuplex) {
                addArg(getArgumentName() + "=" + getTagValue());
            }

            super.runTest();
        }
    }
}
