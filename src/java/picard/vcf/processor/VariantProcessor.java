/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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
package picard.vcf.processor;

import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.variantcontext.VariantContext;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * Describes an object that processes variants and produces a result.
 * <p/>
 * A consumer typically builds an instance of this class via {@link Builder}, providing it the appropriate {@link AccumulatorGenerator} and
 * {@link ResultMerger}, then calls {@link #process()} to obtain the {@link RESULT} of the processing.
 * <p/>
 * Future work...?
 * - Make more efficient for the single-thread case.
 * - A {@link VcfFileSegmentGenerator} that is based on an interval list, so that segments' span a constant-size total-base-count overlap with
 * the interval list (or something in that vein).
 *
 * @author mccowan
 */
public class VariantProcessor<RESULT, ACCUMULATOR extends VariantProcessor.Accumulator<RESULT>> {

    /**
     * Handles {@link VariantContext}s, and accumulates their data in some fashion internally.
     * A call to {@link #result()} produces an embodiment of the results of this processing (which may or may not be the accumulator itself).
     *
     * @author mccowan
     */
    public static interface Accumulator<RESULT> {
        void accumulate(final VariantContext vc);

        RESULT result();
    }

    /**
     * Generates instances of {@link Accumulator}s.
     *
     * @author mccowan
     */
    public static interface AccumulatorGenerator<ACCUMULATOR extends Accumulator<RESULT>, RESULT> {
        ACCUMULATOR build();
    }

    /**
     * Takes a collection of results produced by {@link Accumulator#result()} and merges them into a single {@link RESULT}.
     *
     * @author mccowan
     */
    public static interface ResultMerger<RESULT> {
        RESULT merge(final Collection<RESULT> resultsToReduce);
    }


    final ResultMerger<RESULT> merger;
    final VariantAccumulatorExecutor<ACCUMULATOR, RESULT> executor;

    VariantProcessor(
            final ResultMerger<RESULT> merger,
            final VariantAccumulatorExecutor<ACCUMULATOR, RESULT> executor) {
        this.merger = merger;
        this.executor = executor;
    }

    public RESULT process() {
        executor.start();
        try {
            executor.awaitCompletion();
        } catch (final InterruptedException e) {
            throw new RuntimeException(e);
        }

        final List<RESULT> results = new ArrayList<RESULT>();
        for (final ACCUMULATOR a : executor.accumulators()) {
            results.add(a.result());
        }
        return merger.merge(results);
    }

    /** Simple builder of {@link VariantProcessor}s. */
    public static class Builder<A extends Accumulator<R>, R> {
        final AccumulatorGenerator<A, R> accumulatorGenerator;
        ResultMerger<R> reducer = null;
        IntervalList intervals = null;
        final List<File> inputs = new ArrayList<File>();
        int threadCount = 1;

        Builder(final AccumulatorGenerator<A, R> accumulatorGenerator) {
            this.accumulatorGenerator = accumulatorGenerator;
        }

        public Builder<A, R> multithreadingBy(final int threadCount) {
            if (threadCount < 1) throw new IllegalArgumentException("Multithreading value must exceed 0.");
            this.threadCount = threadCount;
            return this;
        }

        public Builder<A, R> withInput(final File... vcfs) {
            Collections.addAll(inputs, vcfs);
            return this;
        }

        public Builder<A, R> limitingProcessedRegionsTo(final IntervalList intervals) {
            if (this.intervals != null) throw new IllegalStateException("Already provided an interval list.");
            this.intervals = IntervalList.copyOf(intervals);
            return this;
        }

        public Builder<A, R> combiningResultsBy(final ResultMerger<R> reducer) {
            if (this.reducer != null) throw new IllegalStateException("Already provided a reducer.");
            this.reducer = reducer;
            return this;
        }

        public static <A extends Accumulator<R>, R> Builder<A, R> generatingAccumulatorsBy(final AccumulatorGenerator<A, R> generator) {
            return new Builder<A, R>(generator);
        }

        public VariantProcessor<R, A> build() {
            if (inputs.isEmpty()) throw new IllegalStateException("You need to provided some inputs before building.");
            if (reducer == null) throw new IllegalStateException("You must provide a reducer before building.");

            return new VariantProcessor<R, A>(reducer, new VariantAccumulatorExecutor.MultiThreadedChunkBased<A, R>(
                    threadCount,
                    composeVcfIteratorProducerFromBuilderArguments(),
                    accumulatorGenerator
            ));
        }

        private VariantIteratorProducer composeVcfIteratorProducerFromBuilderArguments() {
            /**
             * Be careful; if we pick chunkings that are highly granular (e.g., a chunking based on each interval in an exome-like 
             * interval list), it will result in a {@link htsjdk.variant.vcf.VCFFileReader#query(String, int, int)} call
             * per tiny chunk, which is very non-performant due to some implementations of that method.
             */
            final VariantIteratorProducer ret;
            if (intervals == null) {
                ret = VariantIteratorProducer.byHundredMegabaseChunks(inputs);
            } else {
                ret = VariantIteratorProducer.byHundredMegabaseChunksWithOnTheFlyFilteringByInterval(inputs, intervals);
            }
            return ret;
        }
    }
}
