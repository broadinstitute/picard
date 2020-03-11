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

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * @author mccowan
 */
public class WidthLimitingDecoratorTest {

    class Segment extends VcfFileSegment {
        final int start, stop;

        Segment(final int start, final int stop) {
            this.start = start;
            this.stop = stop;
        }

        @Override
        public int start() {
            return start;
        }

        @Override
        public int stop() {
            return stop;
        }

        @Override
        public String contig() {
            return "A";
        }

        @Override
        public File vcf() {
            return new File("B");
        }

        @Override
        public String toString() {
            return "Segment{" +
                    "start=" + start +
                    ", stop=" + stop +
                    ", vcf=" + vcf() + 
                    ", contig=" + contig() +
                    '}';
        }
    }


    @Test
    public void testForVcf() throws Exception {

        final Segment entireThing = new Segment(1, 9942);
        final ImmutableList<Segment> expectedSubThings = ImmutableList.of(
                new Segment(1, 1000),
                new Segment(1001, 2000),
                new Segment(2001, 3000),
                new Segment(3001, 4000),
                new Segment(4001, 5000),
                new Segment(5001, 6000),
                new Segment(6001, 7000),
                new Segment(7001, 8000),
                new Segment(8001, 9000),
                new Segment(9001, 9942)
        );

        final VcfFileSegmentGenerator.WidthLimitingDecorator strategy = VcfFileSegmentGenerator.WidthLimitingDecorator.wrapping(new VcfFileSegmentGenerator() {
            @Override
            public Iterable<VcfFileSegment> forVcf(final File vcf) {
                return Collections.singleton((VcfFileSegment) entireThing);
            }
        }, 1000);

        final List<VcfFileSegment> observed = new ArrayList<VcfFileSegment>();
        Iterators.addAll(observed, strategy.forVcf(new File("B")).iterator());
        final Iterator<VcfFileSegment> observedIterator = observed.iterator();
        for (final VcfFileSegment e : expectedSubThings) {
            Assert.assertTrue(observedIterator.hasNext());
            final VcfFileSegment o = observedIterator.next();
            Assert.assertEquals(ComparisonChain.start()
                    .compare(o.contig(), e.contig())
                    .compare(o.start(), e.start())
                    .compare(o.stop(), e.stop())
                    .compare(o.vcf(), e.vcf())
                    .result(), 0, String.format(String.format("observed=%s@%s:%s-%s  expected=%s", o.vcf(), o.contig(), o.start(), 
                    o.stop(), e.toString())));
        }
    }
}
