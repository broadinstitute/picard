/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
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
    ERROR(SimpleErrorCalculator::new, "Collects the average (SNP) error at the bases provided."),
    OVERLAPPING_ERROR(OverlappingReadsErrorCalculator::new, "Only considers bases from the overlapping parts of reads from the same template. " +
            "For those bases, it calculates the error that can be attributable to pre-sequencing, versus during-sequencing."),
    INDEL_ERROR(IndelErrorCalculator::new, "Collects insertion and deletion errors at the bases provided.");
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
