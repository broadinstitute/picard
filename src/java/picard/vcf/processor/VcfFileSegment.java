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

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.SAMSequenceRecord;

import java.io.File;

/**
 * Describes a segment of a particular VCF file.
 * 
 * @author mccowan
 */
public abstract class VcfFileSegment {
    abstract public int start();
    abstract public int stop();
    abstract public String contig();
    abstract public File vcf();
    
    public Interval correspondingInterval() {
        return new Interval(contig(), start(), stop());
    }

    @Override
    public String toString() {
        return vcf().getName() + "::" + contig() + ":" + start() + "-" + stop();
    }

    static VcfFileSegment ofWholeSequence(final SAMSequenceRecord sequence, final File vcf) {
        return new SequenceSizedChunk(sequence, vcf);
    }
    
    static class SequenceSizedChunk extends VcfFileSegment {
        final SAMSequenceRecord sequence;
        final File vcf;

        private SequenceSizedChunk(final SAMSequenceRecord sequence, final File vcf) {
            this.sequence = sequence;
            this.vcf = vcf;
        }

        @Override
        public String toString() {
            return vcf().getAbsolutePath() + "::" + sequence.getSequenceName() + ":1-" + sequence.getSequenceLength();
        }

        @Override
        public int start() {
            return 1;
        }

        @Override
        public int stop() {
            return sequence.getSequenceLength();
        }

        @Override
        public String contig() {
            return sequence.getSequenceName();
        }

        @Override
        public File vcf() {
            return vcf;
        }
    }
}
