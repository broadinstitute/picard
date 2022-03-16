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

import com.google.common.primitives.Ints;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import java.nio.file.Path;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Describes a mechanism for producing {@link VcfPathSegment}s from a VCF file path.
 *
 */
public abstract class VcfPathSegmentGenerator {
    final static Log LOG = Log.getInstance(VcfPathSegmentGenerator.class);

    public abstract Iterable<VcfPathSegment> forVcf(final Path vcf);

    public static VcfPathSegmentGenerator byWholeContigSubdividingWithWidth(final long segmentWidth) {
        return WidthLimitingDecorator.wrapping(ByWholeContig.getInstance(), segmentWidth);
    }

    /**
     * Returns a decorated {@link VcfPathSegmentGenerator} that filters out {@link VcfPathSegment}s that have no overlap with the provided
     * {@link OverlapDetector}.
     */
    public static <T> VcfPathSegmentGenerator excludingNonOverlaps(final VcfPathSegmentGenerator strategy, final OverlapDetector<T> overlaps) {
        return new VcfPathSegmentGenerator() {
            @Override
            public Iterable<VcfPathSegment> forVcf(final Path vcf) {
                return StreamSupport.stream(strategy.forVcf(vcf).spliterator(), false).filter(segment -> {
                    final boolean keep = !overlaps.getOverlaps(new Interval(segment.contig(), segment.start(), segment.stop())).isEmpty();
                    if (!keep) {
                        LOG.debug(String.format("Ignoring segment because it does not overlap with detector, %s::%s:%s-%s",
                                segment.vcf().getFileName(), segment.contig(), segment.start(), segment.stop())
                        );
                    }
                    return keep;
                }).collect(Collectors.toList());
            }
        };
    }

    /**
     * A very simple {@link VcfPathSegmentGenerator} that breaks up the provided vcfs into contig-sized chunks.
     *
     * @author mccowan
     */
    static class ByWholeContig extends VcfPathSegmentGenerator {
        // Singleton!
        ByWholeContig() {
        }

        private static final ByWholeContig singleton = new ByWholeContig();

        public static ByWholeContig getInstance() {
            return singleton;
        }

        @Override
        public Iterable<VcfPathSegment> forVcf(final Path vcf) {
            final List<SAMSequenceRecord> samSequenceRecords = readSequences(vcf);
            return samSequenceRecords.stream().map(samSequenceRecord -> VcfPathSegment.ofWholeSequence(samSequenceRecord, vcf)).collect(Collectors.toList());
        }

        private static List<SAMSequenceRecord> readSequences(final Path vcf) {
            final VCFFileReader reader = new VCFFileReader(vcf);
            final VCFHeader header = reader.getFileHeader();
            final SAMSequenceDictionary dict = header.getSequenceDictionary();
            reader.close();
            return dict.getSequences();
        }
    }

    /**
     * Decorator to apply to other {@link VcfPathSegmentGenerator} to enforce that no segment is larger than the specified width.
     *
     * @author mccowan
     */
    static final class WidthLimitingDecorator extends VcfPathSegmentGenerator {
        final VcfPathSegmentGenerator underlyingStrategy;
        final long width;

        public static WidthLimitingDecorator wrapping(final VcfPathSegmentGenerator basis, final long maximumWidth) {
            return new WidthLimitingDecorator(basis, maximumWidth);
        }

        private WidthLimitingDecorator(final VcfPathSegmentGenerator underlyingStrategy, final long maximumWidth) {
            this.underlyingStrategy = underlyingStrategy;
            this.width = maximumWidth - 1;
        }

        /**
         * The thing that does the work; accepts a {@link VcfPathSegment} (produced by the parent {@link VcfPathSegmentGenerator}) and breaks
         * it down into subsegments.
         */
        private final class VcfPathSegmentSubdivider implements Iterable<VcfPathSegment> {
            final VcfPathSegment basis;

            private VcfPathSegmentSubdivider(final VcfPathSegment basis) {
                this.basis = basis;
            }

            @Override
            public Iterator<VcfPathSegment> iterator() {
                return new Iterator<VcfPathSegment>() {
                    int nextStart = basis.start();

                    @Override
                    public boolean hasNext() {
                        return nextStart <= basis.stop();
                    }

                    @Override
                    public VcfPathSegment next() {
                        final int start = nextStart;
                        final VcfPathSegment ret = new VcfPathSegment() {
                            @Override
                            public int start() {
                                return start;
                            }

                            @Override
                            public int stop() {
                                return Ints.checkedCast(Math.min(start + width, basis.stop()));
                            }

                            @Override
                            public String contig() {
                                return basis.contig();
                            }

                            @Override
                            public Path vcf() {
                                return basis.vcf();
                            }
                        };
                        nextStart += width + 1;
                        return ret;
                    }

                    @Override
                    public void remove() {
                        throw new UnsupportedOperationException();
                    }
                };
            }
        }

        @Override
        public Iterable<VcfPathSegment> forVcf(final Path vcf) {
            // Turn the VCF into segments, and then apply our 
            return StreamSupport.stream(underlyingStrategy.forVcf(vcf).spliterator(), false)
                    .flatMap(vcfFileSegment -> StreamSupport.stream(new VcfPathSegmentSubdivider(vcfFileSegment).spliterator(), false)).collect(Collectors.toList());
        }
    }

}
