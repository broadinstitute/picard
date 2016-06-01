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
import com.google.common.base.Predicate;
import com.google.common.collect.FluentIterable;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import picard.vcf.processor.util.PredicateFilterDecoratingClosableIterator;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * A mechanism for iterating over {@link CloseableIterator} of {@link VariantContext}s in in some fashion, given VCF files and optionally
 * an interval list.
 * 
 * The produced iterators may perform on-the-fly filtering of the produced {@link VariantContext}s.
 *
 * @author mccowan
 */
public abstract class VariantIteratorProducer {
    final static int ONE_HUNDRED_MILLION = (int) 100e6;
    /** 
     * Renders the embodied regions of the VCF files in the form of {@link htsjdk.samtools.util.CloseableIterator}s over
     * {@link VariantContext}s.  The iterator may perform on-the-fly filtering of these elements.
     */
    public abstract Iterable<CloseableIterator<VariantContext>> iterators();

    /** Closes any latent file handles that may have been opened by calls to {@link #iterators()}. */
    public abstract void close();

    /**
     * Produces a chunking with segments of size 100 megabases (or less if a contig boundary is reached), that also performs on-the-fly
     * filtering of {@link VariantContext}
     */
    public static VariantIteratorProducer byHundredMegabaseChunksWithOnTheFlyFilteringByInterval(final List<File> vcfs, final IntervalList intervalList) {
        return new Threadsafe(VcfFileSegmentGenerator.byWholeContigSubdividingWithWidth(ONE_HUNDRED_MILLION), vcfs, intervalList);
    }

    /** Produces a chunking with segments of size 100 megabases (or less if a contig boundary is reached). */
    public static VariantIteratorProducer byHundredMegabaseChunks(final List<File> vcfs) {
        return new Threadsafe(VcfFileSegmentGenerator.byWholeContigSubdividingWithWidth(ONE_HUNDRED_MILLION), vcfs, null);
    }

    /**
     * A {@link VariantIteratorProducer} that is based on a given {@link VcfFileSegmentGenerator} and a list of VCFs.  The chunks are ordered by VCF, and
     * then by whatever ordering of segments are produced by {@link VcfFileSegmentGenerator#forVcf(java.io.File)} for each of those VCFs.
     * <p/>
     * The iterators produced by this class are safe to share between multiple threads.
     * <p/>
     * If an {@link IntervalList} provided, the produced iterators will perform on-the-fly filtering and calls to {@link Iterator#next()} will
     * only include {@link VariantContext}s that fall within the regions described by that list.
     * <p/>
     * This class maintains a {@link ThreadLocal} of {@link VCFFileReader} to ensure that each thread has its own, and at most one, reader per
     * VCF.  It also guarantees that, so long as a thread closes its "queried into" readers after expiring them, there is only one extant
     * "queried into" iterator per thread per VCF.
     *
     * @author mccowan
     */
    static class Threadsafe extends VariantIteratorProducer {
        final static Log LOG = Log.getInstance(Threadsafe.class);

        /** A list of the segments for which the corresponding {@link VariantContext}s will be produced. */
        final List<VcfFileSegment> segments;
        final OverlapDetector<Interval> intervalsOfInterestDetector;

        /** Maps directly to {@link #segments}; useful for determining if a given variant falls into multiple segments (don't double-count!). */
        final Map<File, OverlapDetector<VcfFileSegment>> multiSegmentDetectorPerFile =
                new CollectionUtil.DefaultingMap<File,OverlapDetector<VcfFileSegment>>(new CollectionUtil.DefaultingMap.Factory<OverlapDetector<VcfFileSegment>, File>() {
                    @Override
                    public OverlapDetector<VcfFileSegment> make(final File f) {
                        return new OverlapDetector<VcfFileSegment>(0, 0);
                    }
                }, true);


        Threadsafe(final VcfFileSegmentGenerator segmenter, final List<File> vcfs) {
            this(segmenter, vcfs, null);
        }

        Threadsafe(final VcfFileSegmentGenerator segmenter, final List<File> vcfs, final IntervalList intervals) {
            if (intervals != null) {
                final List<Interval> uniques = intervals.getUniqueIntervals(false);
                this.intervalsOfInterestDetector = new OverlapDetector<Interval>(0, 0);
                intervalsOfInterestDetector.addAll(uniques, uniques);
            } else {
                intervalsOfInterestDetector = null;
            }

            /**
             * Prepare {@link #segments} and {@link multiSegmentDetectorPerFile}.  First, we only want to accumulate segments that are
             * interesting, e.g., only ones that do not lie completely outside of the interval list provided by the consumer (if it was
             * provided).
             *
             * Then, break up our {@link VcfFileSegment}s by file, and build a corresponding {@link OverlapDetector} for those segments.  This
             * will be used downstream to ensure that we don't emit a given VC from a VCF more than once.
             *
             * While we're at it, since this class expects the provided {@link VcfFileSegmentGenerator} produces non-overlapping segments, 
             * assert that this is true.
             */
            final VcfFileSegmentGenerator interestingSegmentSegmenter =
                    intervalsOfInterestDetector == null ? segmenter : VcfFileSegmentGenerator.excludingNonOverlaps(segmenter, intervalsOfInterestDetector);
            segments = new ArrayList<VcfFileSegment>();
            for (final File vcf : vcfs) {
                for (final VcfFileSegment segment : interestingSegmentSegmenter.forVcf(vcf)) {
                    segments.add(segment);
                }
            }
            for (final VcfFileSegment segment : segments) {
                final Interval segmentInterval = segment.correspondingInterval();
                final OverlapDetector<VcfFileSegment> vcfSpecificDetector = multiSegmentDetectorPerFile.get(segment.vcf());
                if (vcfSpecificDetector.getOverlaps(segmentInterval).isEmpty()) {
                    vcfSpecificDetector.addLhs(segment, new Interval(segment.contig(), segment.start(), segment.stop()));
                } else {
                    throw new IllegalArgumentException(String.format(
                            "Provided segmenting strategy produced overlapping intervals; %s overlaps with: %s",
                            segment,
                            Joiner.on(", ").join(vcfSpecificDetector.getOverlaps(segmentInterval))
                    ));
                }
            }
        }

        /**
         * All of the {@link VCFFileReader}s opened by this object, across all threads.  (We can't get this
         * value out of {@link #localVcfFileReaders} because of the way {@link java.lang.ThreadLocal} works.
         */
        final Collection<VCFFileReader> allReaders = Collections.synchronizedCollection(new ArrayList<VCFFileReader>());

        /**
         * A map maintained for each thread that contains the {@link htsjdk.variant.vcf.VCFFileReader}s that it has opened for a
         * given VCF.  (Threads never share {@link htsjdk.variant.vcf.VCFFileReader}s.
         * <p/>
         * This is a {@link CollectionUtil.DefaultingMap} to effectively produce readers on-the-fly.  Note that it also adds every produced
         * reader into {@link #allReaders}.
         */
        final ThreadLocal<CollectionUtil.DefaultingMap<File, VCFFileReader>> localVcfFileReaders =
                new ThreadLocal<CollectionUtil.DefaultingMap<File, VCFFileReader>>() {
                    @Override
                    protected CollectionUtil.DefaultingMap<File, VCFFileReader> initialValue() {
                        return new CollectionUtil.DefaultingMap<File, VCFFileReader>(new CollectionUtil.DefaultingMap.Factory<VCFFileReader, File>() {
                            @Override
                            public VCFFileReader make(final File file) {
                                final VCFFileReader reader = new VCFFileReader(file);
                                LOG.debug(String.format("Producing a reader of %s for %s.", file, Thread.currentThread()));
                                allReaders.add(reader);
                                return reader;
                            }
                        }, true);
                    }
                };

        /**
         * Converts a {@link VcfFileSegment} into a {@link VariantContext} iterator.  Applies filtering via {@link #intervalsOfInterestDetector}
         * if it is defined.
         */
        private CloseableIterator<VariantContext> iteratorForSegment(final VcfFileSegment segment) {
            final CloseableIterator<VariantContext> query =
                    localVcfFileReaders.get() // Get the collection of VCF file readers local to this thread
                            .get(segment.vcf()) // Get or generate the reader for this segment's VCF file
                            .query(segment.contig(), segment.start(), segment.stop()); // Query the segment

            // Then wrap the iterator in a on-the-fly interval-list based filter, if requested.
            final Collection<Predicate<VariantContext>> filters = new ArrayList<Predicate<VariantContext>>();
            if (intervalsOfInterestDetector != null) {
                filters.add(new OverlapsPredicate());
            }
            filters.add(new NonUniqueVariantPredicate(segment));
            return new PredicateFilterDecoratingClosableIterator<VariantContext>(query, filters);
        }

        @Override
        public Iterable<CloseableIterator<VariantContext>> iterators() {
            return FluentIterable.from(segments).transform(new Function<VcfFileSegment, CloseableIterator<VariantContext>>() {
                @Override
                public CloseableIterator<VariantContext> apply(final VcfFileSegment segment) {
                    return iteratorForSegment(segment);
                }
            });
        }

        @Override
        public void close() {
            final Iterator<VCFFileReader> i = allReaders.iterator();
            while (i.hasNext()) {
                i.next().close();
                i.remove();
            }
        }

        /**
         * A predicate that I had difficulty naming. The value of this predicate is that it ensures that no single variant is produced multiple
         * times from a call to {@link Threadsafe#iterators()}.  It works by asking each variant "Hey variant, which of the
         * {@link VcfFileSegment}s that this {@link Threadsafe} is producing variants from do you fall into (and thus,
         * would be produced by a corresponding call to {@link VCFFileReader#query(String, int, int)}?  If it's more than one, we're only going
         * to emit you from the query for the very first of those {@link VcfFileSegment}s."
         */
        final class NonUniqueVariantPredicate implements Predicate<VariantContext> {
            final VcfFileSegment sourceSegment;

            NonUniqueVariantPredicate(final VcfFileSegment sourceSegment) {
                this.sourceSegment = sourceSegment;
            }

            @Override
            public boolean apply(final VariantContext vc) {
                if (vc.getStart() == vc.getEnd()) {
                    // This isn't a multi-allelic segment, so it can only be produced by a single segment.
                    return true;
                }
                final Collection<VcfFileSegment> intersectingSegments =
                        multiSegmentDetectorPerFile.get(sourceSegment.vcf()).getOverlaps(new Interval(vc.getChr(), vc.getStart(), vc.getEnd()));
                if (intersectingSegments.size() < 2) {
                    // There's only one segment that produces this variant
                    return true;
                }

                // The convention is: only emit the VC if it is produced from the first segment that can produce it in the list.
                final int sourceSegmentIndex = segments.indexOf(sourceSegment);
                LOG.debug("Found wide variant spanning multiple source segments: ", vc);
                for (final VcfFileSegment intersectingSegment : intersectingSegments) {
                    if (segments.indexOf(intersectingSegment) < sourceSegmentIndex) {
                        // There is a segment that produces this variant earlier in the segment list, exclude it.
                        return false;
                    }
                }
                LOG.debug("Emitting wide variant because it belongs to first segment: ", vc);
                return true;
            }
        }

        /** A predicate answering if the provided {@link VariantContext} overlaps {@link #intervalsOfInterestDetector}. */
        final class OverlapsPredicate implements Predicate<VariantContext> {
            @Override
            public boolean apply(final VariantContext vc) {
                final boolean include = !intervalsOfInterestDetector.getOverlaps(new Interval(vc.getChr(), vc.getStart(), vc.getEnd())).isEmpty();
                if (!include) LOG.debug("Filtering variant at ", vc.getChr(), ":", vc.getStart(), "-", vc.getEnd());
                return include;
            }
        }
    }

}
