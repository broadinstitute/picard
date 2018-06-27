package picard.sam.SamErrorMetric;

import org.broadinstitute.barclay.argparser.CommandLineParser;

import java.util.function.Supplier;

/**
 * An enum that is used to generate a {@link Supplier <BaseErrorCalculator>} from a string
 * To use this given a String 'str':
 * <p>
 * Errors.valueOf(str).getErrorSupplier()
 * <p>
 * This is used in {@link CollectSamErrorMetrics} to convert an input argument to a {@link BaseErrorAggregation}.
 */
enum ErrorType implements CommandLineParser.ClpEnum {
    ERROR(SimpleErrorCalculator::new, "Collects the average error at the bases provided."),
    OVERLAPPING_ERROR(OverlappingReadsErrorCalculator::new, "Only considers bases from the overlapping parts of reads from the same template. " +
            "For those bases, it calculates the error that can be attributable to pre-sequencing, versus during-sequencing.");

    private final Supplier<? extends BaseCalculator> errorSupplier;

    ErrorType(Supplier<? extends BaseCalculator> errorSupplier, final String docString) {
        this.errorSupplier = errorSupplier;
        this.docString = docString;
    }

    public Supplier<? extends BaseCalculator> getErrorSupplier() {
        return errorSupplier;
    }

    private final String docString;

    @Override
    public String getHelpDoc() {
        return docString + " Suffix is: '" + errorSupplier.get().getSuffix() + "'.";
    }
}
