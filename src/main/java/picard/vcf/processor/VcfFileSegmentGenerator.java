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

import java.io.File;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Describes a mechanism for producing {@link VcfFileSegment}s from a VCF file.
 *
 * @author mccowan
 */
public abstract class VcfFileSegmentGenerator {
    final static Log LOG = Log.getInstance(VcfFileSegmentGenerator.class);

    public abstract Iterable<VcfFileSegment> forVcf(final File vcf);

    public static VcfFileSegmentGenerator byWholeContigSubdividingWithWidth(final long segmentWidth) {
        return WidthLimitingDecorator.wrapping(ByWholeContig.getInstance(), segmentWidth);
    }

    /**
     * Returns a decorated {@link VcfFileSegmentGenerator} that filters out {@link VcfFileSegment}s that have no overlap with the provided
     * {@link OverlapDetector}.
     */
    public static <T> VcfFileSegmentGenerator excludingNonOverlaps(final VcfFileSegmentGenerator strategy, final OverlapDetector<T> overlaps) {
        return  new VcfFileSegmentGenerator() {
            @Override
            public Iterable<VcfFileSegment> forVcf(final File vcf) {
                return StreamSupport.stream(strategy.forVcf(vcf).spliterator(), false).filter(segment -> {
                    final boolean keep = !overlaps.getOverlaps(new Interval(segment.contig(), segment.start(), segment.stop())).isEmpty();
                    if (!keep) {
                        LOG.debug(String.format("Ignoring segment because it does not overlap with detector, %s::%s:%s-%s",
                                segment.vcf().getName(), segment.contig(), segment.start(), segment.stop())
                        );
                    }
                    return keep;
                }).collect(Collectors.toList());
            }
        };
    }

    /**
     * A very simple {@link VcfFileSegmentGenerator} that breaks up the provided vcfs into contig-sized chunks.
     *
     * @author mccowan
     */
    static class ByWholeContig extends VcfFileSegmentGenerator {
        // Singleton!
        ByWholeContig() {
        }

        private static final ByWholeContig singleton = new ByWholeContig();

        public static ByWholeContig getInstance() {
            return singleton;
        }

        @Override
        public Iterable<VcfFileSegment> forVcf(final File vcf) {
            final List<SAMSequenceRecord> samSequenceRecords = readSequences(vcf);
            return samSequenceRecords.stream().map(samSequenceRecord -> VcfFileSegment.ofWholeSequence(samSequenceRecord, vcf)).collect(Collectors.toList());
        }

        private static List<SAMSequenceRecord> readSequences(final File vcf) {
            final VCFFileReader reader = new VCFFileReader(vcf);
            final VCFHeader header = reader.getFileHeader();
            final SAMSequenceDictionary dict = header.getSequenceDictionary();
            reader.close();
            return dict.getSequences();
        }
    }

    /**
     * Decorator to apply to other {@link VcfFileSegmentGenerator} to enforce that no segment is larger than the specified width.
     *
     * @author mccowan
     */
    static final class WidthLimitingDecorator extends VcfFileSegmentGenerator {
        final VcfFileSegmentGenerator underlyingStrategy;
        final long width;

        public static WidthLimitingDecorator wrapping(final VcfFileSegmentGenerator basis, final long maximumWidth) {
            return new WidthLimitingDecorator(basis, maximumWidth);
        }

        private WidthLimitingDecorator(final VcfFileSegmentGenerator underlyingStrategy, final long maximumWidth) {
            this.underlyingStrategy = underlyingStrategy;
            this.width = maximumWidth - 1;
        }

        /**
         * The thing that does the work; accepts a {@link VcfFileSegment} (produced by the parent {@link VcfFileSegmentGenerator}) and breaks
         * it down into subsegments.
         */
        private final class VcfFileSegmentSubdivider implements Iterable<VcfFileSegment> {
            final VcfFileSegment basis;

            private VcfFileSegmentSubdivider(final VcfFileSegment basis) {
                this.basis = basis;
            }

            @Override
            public Iterator<VcfFileSegment> iterator() {
                return new Iterator<VcfFileSegment>() {
                    int nextStart = basis.start();

                    @Override
                    public boolean hasNext() {
                        return nextStart <= basis.stop();
                    }

                    @Override
                    public VcfFileSegment next() {
                        final int start = nextStart;
                        final VcfFileSegment ret = new VcfFileSegment() {
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
                            public File vcf() {
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
        public Iterable<VcfFileSegment> forVcf(final File vcf) {
            // Turn the VCF into segments, and then apply our 
            return StreamSupport.stream(underlyingStrategy.forVcf(vcf).spliterator(), false).flatMap(vcfFileSegment -> StreamSupport.stream(new VcfFileSegmentSubdivider(vcfFileSegment).spliterator(), false)).collect(Collectors.toList());
        }
    }

}
