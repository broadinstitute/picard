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

import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.base.Optional;
import com.google.common.collect.FluentIterable;
import picard.util.AtomicIterator;
import picard.util.Iterators;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Describes the functionality for an executor that manages the delegation of work to {@link VariantProcessor.Accumulator}s.
 * 
 * @author mccowan
 */
public interface VariantAccumulatorExecutor<ACCUMULATOR extends VariantProcessor.Accumulator<RESULT>, RESULT> {
    /** Starts the work of the executor, returning immediately. */
    void start();

    /** Blocks until the work is complete. */
    void awaitCompletion() throws InterruptedException;
    
    /** Returns the {@link VariantProcessor.Accumulator}s associated with this executor. */
    Collection<ACCUMULATOR> accumulators();

    /**
     * A {@link VariantAccumulatorExecutor} that breaks down work into chunks described by the provided {@link VariantIteratorProducer} and
     * spreads them over the indicated number of threads.
     *
     * @author mccowan
     */
    class MultiThreadedChunkBased<A extends VariantProcessor.Accumulator<R>, R> implements VariantAccumulatorExecutor<A, R> {
        private static final Log LOG = Log.getInstance(MultiThreadedChunkBased.class);

        final AtomicIterator<CloseableIterator<VariantContext>> vcIterators;
        final ExecutorService executor;
        final Collection<A> accumulators = Collections.synchronizedCollection(new ArrayList<A>());

        /** Signals whether or not this executor is started. */
        volatile boolean started = false;
        final int numThreads;

        private final List<Throwable> childrenErrors = Collections.synchronizedList(new ArrayList<Throwable>());

        final VariantProcessor.AccumulatorGenerator<A, R> accumulatorGenerator;

        public MultiThreadedChunkBased(
                final int numThreads,
                final VariantIteratorProducer vcIteratorProducer,
                final VariantProcessor.AccumulatorGenerator<A, R> accumulatorGenerator
        ) {
            this.executor = Executors.newFixedThreadPool(numThreads);
            this.vcIterators = Iterators.atomicIteratorOf(vcIteratorProducer.iterators());
            this.numThreads = numThreads;
            this.accumulatorGenerator = accumulatorGenerator;
        }


        @Override
        public synchronized void start() {
            started = true;
            for (int i = 0; i < numThreads; i++) {
                final A accumulator = accumulatorGenerator.build();
                accumulators.add(accumulator);
                executor.submit(new Worker(accumulator));
            }
            executor.shutdown();
        }

        public synchronized Collection<A> accumulators() {
            return Collections.unmodifiableCollection(accumulators);
        }

        @Override
        public void awaitCompletion() throws InterruptedException {
            if (!started) {
                throw new IllegalStateException("This method can be called only after the executor has been started.");
            } else {
                executor.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
                if (!childrenErrors.isEmpty()) {
                    throw new MultiException(childrenErrors);
                }
            }
        }

        static class MultiException extends RuntimeException {
            final List<Throwable> childrenExceptions;

            public MultiException(final List<Throwable> childrenExceptions) {
                this.childrenExceptions = childrenExceptions;
            }

            @Override
            public String getMessage() {
                return "Children threads encountered exceptions:\n" + Joiner.on("\n\t").join(FluentIterable.from(childrenExceptions).transform
                        (new Function<Throwable, String>() {
                            @Override
                            public String apply(final Throwable throwable) {
                                return throwable.getMessage();
                            }
                        }));
            }
        }

        /** Continually requests and exhausts variant context iterators, delegating each to the child {@link Worker#processor}. */
        class Worker implements Runnable {
            final VariantProcessor.Accumulator processor;

            Worker(final VariantProcessor.Accumulator processor) {
                this.processor = processor;
            }

            @Override
            public void run() {
                try {
                    Optional<CloseableIterator<VariantContext>> readerMaybe;
                    while ((readerMaybe = vcIterators.next()).isPresent()) {
                        final CloseableIterator<VariantContext> reader = readerMaybe.get();
                        while (reader.hasNext()) processor.accumulate(reader.next());
                        reader.close();

                        if (!childrenErrors.isEmpty()) {
                            LOG.error(Thread.currentThread() + " aborting: observed error in another child thread.");
                            break;
                        }
                    }
                } catch (final Throwable e) {
                    childrenErrors.add(e);
                    LOG.error(e, "Unexpected exception encountered in child thread.");
                } finally {
                    LOG.debug(String.format("Thread %s is finishing.", Thread.currentThread()));
                }
            }
        }
    }
}
