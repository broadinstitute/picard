/*
* Copyright (c) 2013 The Broad Institute
*
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
*
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package net.sf.picard.vcf;

import net.sf.picard.PicardException;
import net.sf.samtools.util.CloserUtil;
import org.broad.tribble.FeatureCodecHeader;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.variant.bcf2.BCF2Codec;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFHeader;

import java.io.*;
import java.util.NoSuchElementException;

public class BcfIterator implements VariantContextIterator {
    private final BCF2Codec bcfCodec = new BCF2Codec();
    private final PositionalBufferedStream inputStream;
    private final FeatureCodecHeader codecHeader;

    public BcfIterator(final InputStream bcfStream) {
        inputStream = new PositionalBufferedStream(bcfStream);
        codecHeader = bcfCodec.readHeader(inputStream);
    }

    @Override
    public void close() {
        CloserUtil.close(inputStream);
    }

    @Override
    public boolean hasNext() {
        final boolean isDone;
        try {
            isDone = inputStream.isDone();
        } catch (IOException ioe) {
            throw new PicardException("Unable to determine if BcfIterator is exhausted", ioe);
        }
        return !isDone;
    }

    @Override
    public VariantContext next() {
        if (!this.hasNext()) {
            throw new NoSuchElementException("Called next() on an exhausted BcfIterator");
        }
        return bcfCodec.decode(inputStream);
    }

    public VCFHeader getHeader() {
        return (VCFHeader)codecHeader.getHeaderValue();
    }

    /**
     * Unsupported.
     */
    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }
}
