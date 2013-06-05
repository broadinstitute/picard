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

import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.util.RuntimeIOException;
import org.broad.tribble.readers.AsciiLineReader;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;
import java.io.IOException;
import java.io.InputStream;
import java.util.NoSuchElementException;

public class VcfIterator implements VariantContextIterator {
    private final VCFCodec vcfCodec = new VCFCodec();
    private final VCFHeader vcfHeader;
    private final AsciiLineReader reader;

    private String line = null;

    public VcfIterator(final InputStream vcfStream) {
        this.reader = new AsciiLineReader(vcfStream);
        final Object header = vcfCodec.readHeader(reader);
        if (!(header instanceof VCFHeader)) {
            throw new IllegalArgumentException("No VCF header found");
        }
        this.vcfHeader = (VCFHeader) header;
    }

    // TODO: Add a c'tor that reads intervals.

    @Override
    public void close() {
        CloserUtil.close(reader);
    }

    public VCFHeader getHeader() {
        return this.vcfHeader;
    }

    @Override
    public boolean hasNext() {
        try {
            if (line == null) line = reader.readLine();
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
        return line != null;
    }

    @Override
    public VariantContext next() {
        if ( ! this.hasNext()) throw new NoSuchElementException("Called next() on an exhausted VcfIterator");
        final String tmp = line;
        line = null;
        return vcfCodec.decode(tmp);
    }

    /**
     * Unsupported.
     */
    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }


}
