/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.variant.vcf;

import java.util.Map;

/**
 * A special class representing a contig VCF header line.  Nows the true contig order and sorts on that
 *
 * @author mdepristo
 */
public class VCFContigHeaderLine extends VCFSimpleHeaderLine {
    final Integer contigIndex;


    /**
     * create a VCF contig header line
     *
     * @param line      the header line
     * @param version   the vcf header version
     * @param key            the key for this header line
     */
    public VCFContigHeaderLine(final String line, final VCFHeaderVersion version, final String key, int contigIndex) {
        super(line, version, key, null);
        this.contigIndex = contigIndex;
    }

    public VCFContigHeaderLine(final Map<String, String> mapping, int contigIndex) {
        super(VCFHeader.CONTIG_KEY, mapping);
        this.contigIndex = contigIndex;
    }

    public Integer getContigIndex() {
        return contigIndex;
    }

    /**
     * IT IS CRITIAL THAT THIS BE OVERRIDDEN SO WE SORT THE CONTIGS IN THE CORRECT ORDER
     *
     * @param other
     * @return
     */
    @Override
    public int compareTo(final Object other) {
        if ( other instanceof VCFContigHeaderLine )
            return contigIndex.compareTo(((VCFContigHeaderLine) other).contigIndex);
        else {
            return super.compareTo(other);
        }
    }
}