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

import htsjdk.samtools.reference.SamLocusAndReferenceIterator.SAMLocusAndReference;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.SamLocusIterator;

import java.util.*;
import java.util.function.Supplier;

/**
 * An interface and implementations for classes that apply a {@link ReadBaseStratification.RecordAndOffsetStratifier RecordAndOffsetStratifier}
 * to put bases into various "bins" and then compute an {@link ErrorMetric} on these bases using a {@link BaseErrorCalculator}.
 *
 */

public class BaseErrorAggregation<CALCULATOR extends BaseCalculator> {
    private final Supplier<CALCULATOR> simpleAggregatorGenerator;
    private final ReadBaseStratification.RecordAndOffsetStratifier stratifier;
    private final Map<Object, CALCULATOR> strataAggregatorMap;

    public BaseErrorAggregation(final Supplier<CALCULATOR> simpleAggregatorGenerator,
                                final ReadBaseStratification.RecordAndOffsetStratifier stratifier) {
        this.stratifier = stratifier;
        this.simpleAggregatorGenerator = simpleAggregatorGenerator;
        this.strataAggregatorMap =  new CollectionUtil.DefaultingMap<>
                (ignored -> simpleAggregatorGenerator.get(),true);
    }

    public void addBase(final SamLocusIterator.RecordAndOffset recordAndOffset, final SAMLocusAndReference locusInfo) {
        final Object stratus = stratifier.stratify(recordAndOffset, locusInfo);
        // this assumes we do not want to aggregate null.
        if (stratus != null) {
            strataAggregatorMap.get(stratus).addBase(recordAndOffset, locusInfo);
        }
    }

    public String getSuffix() {
        return simpleAggregatorGenerator.get().getSuffix() + "_by_" + stratifier.getSuffix();
    }

    public ErrorMetric[] getMetrics() {
        final List<ErrorMetric> metrics = new ArrayList<>();

        // we do this to sort
        final TreeSet<Object> strata = new TreeSet<>(strataAggregatorMap.keySet());

        for (final Object stratum : strata) {
            final ErrorMetric metric = strataAggregatorMap.get(stratum).getMetric();
            metric.COVARIATE = stratum.toString();
            metrics.add(metric);
        }
        return metrics.toArray(new ErrorMetric[0]);
    }
}

