/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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

package picard.fingerprint;

import picard.PicardException;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReaderUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

import java.io.BufferedInputStream;
import java.io.File;
import java.util.Collection;

/**
 * Class that abstracts away the source of genotypes and provides abilities to read them
 * in from various sources into VariantContext objects
 *
 * @deprecated  Please use VCFFileReader instead of this class.
 */
public class GenotypeReader {
    /** Small class to encapsulate an iterator over variants, optionally with a sequence dictionary. */
    public abstract static class VariantIterator implements CloseableIterator<VariantContext> {
        private final SAMSequenceDictionary dictionary;
        private final Object header;

        protected VariantIterator(final SAMSequenceDictionary dictionary, final Object header) {
            this.dictionary = dictionary;
            this.header = header;
        }

        public SAMSequenceDictionary getSequenceDictionary() {
            return this.dictionary;
        }

        public Object getHeader() {
            return header;
        }
    }

    /**
     * Reads in a file that contains genotype data and returns an iterator over every entry
     * in the file as VariantContext objects.
     *
     * @deprecated  Please use VCFFileReader in Picard-public instead of this class.
     */
    public VariantIterator read(final File file) {
        if (isVcf(file)) return readVcf(file);
        else throw new PicardException("File doe not appear to be of a supported type: " + file);
    }

    /**
     * Reads in the file and returns an iterator over the set of entries in the file
     * that overlap with the interval supplied.
     *
     * @deprecated  Please use VCFFileReader in Picard-public instead of this class.
     */
    public VariantIterator read(final File file, final IntervalList intervals) {
        final VariantIterator i = read(file);
        final OverlapDetector<Interval> detector = new OverlapDetector<Interval>(0,0);
        detector.addAll(intervals.getIntervals(), intervals.getIntervals());

        // A little iterator that iterates over the full set of genotypes and returns just the ones
        // the client is interested in.
        return new VariantIterator(i.getSequenceDictionary(), i.getHeader()) {
            private VariantContext next = null;

            @Override public boolean hasNext() {
                if (next == null) {
                    while (i.hasNext()) {
                        final VariantContext ctx = i.next();
                        final Interval ctxInterval = new Interval(ctx.getChr(), ctx.getStart(), ctx.getEnd());
                        final Collection<Interval> hits = detector.getOverlaps(ctxInterval);
                        if (hits != null && !hits.isEmpty()) {
                            next = ctx;
                            break;
                        }
                    }
                }

                return next != null;
            }

            /** Returns the next VariantContext object if available. */
            @Override public VariantContext next() {
                if (!hasNext()) throw new IllegalStateException("next() called on exhausted iterator.");

                final VariantContext ctx = next;
                next = null;
                return ctx;
            }

            @Override public void remove() { throw new UnsupportedOperationException(); }

            @Override
            public void close() {
                i.close();
            }
        };
    }

    /** Tests whether a file is a VCF file or not. */
    boolean isVcf(final File f) {
        final String name = f.getName();
        return (name.endsWith(".vcf") || name.endsWith(".vcf.gz"));
    }

    /**
     * Opens a VCF file and returns an iterator over VariantContext objects.
     *
     * @deprecated  Please use VCFFileReader in Picard-public instead of this class.
     */
    VariantIterator readVcf(final File file) {
        final LineIterator reader = new LineIteratorImpl(LineReaderUtil.fromBufferedStream(new BufferedInputStream(IOUtil.openFileForReading(file))));
        final VCFCodec codec   = new VCFCodec();
        final Object header;
        header = codec.readActualHeader(reader);

        return new VariantIterator(null, header) {
            @Override public boolean hasNext() {
                return reader.hasNext();
            }

            @Override public VariantContext next() {
                return codec.decode(reader.next());
            }

            @Override public void remove() {
                throw new UnsupportedOperationException();
            }

            @Override
            public void close() {
                codec.close(reader);
            }
        };
    }
}
