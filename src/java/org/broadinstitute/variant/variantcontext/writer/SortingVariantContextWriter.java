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

package org.broadinstitute.variant.variantcontext.writer;

import org.broadinstitute.variant.variantcontext.VariantContext;

/**
 * this class writes VCF files, allowing records to be passed in unsorted (up to a certain genomic distance away)
 */
class SortingVariantContextWriter extends SortingVariantContextWriterBase {

    // the maximum START distance between records that we'll cache
    private int maxCachingStartDistance;

    /**
     * create a local-sorting VCF writer, given an inner VCF writer to write to
     *
     * @param innerWriter        the VCFWriter to write to
     * @param maxCachingStartDistance the maximum start distance between records that we'll cache
     * @param takeOwnershipOfInner Should this Writer close innerWriter when it's done with it
     */
    public SortingVariantContextWriter(VariantContextWriter innerWriter, int maxCachingStartDistance, boolean takeOwnershipOfInner) {
        super(innerWriter, takeOwnershipOfInner);
        this.maxCachingStartDistance = maxCachingStartDistance;
    }

    public SortingVariantContextWriter(VariantContextWriter innerWriter, int maxCachingStartDistance) {
        this(innerWriter, maxCachingStartDistance, false); // by default, don't own inner
    }

    protected void noteCurrentRecord(VariantContext vc) {
        super.noteCurrentRecord(vc); // first, check for errors

        // then, update mostUpstreamWritableLoc:
        int mostUpstreamWritableIndex = vc.getStart() - maxCachingStartDistance;
        this.mostUpstreamWritableLoc = Math.max(BEFORE_MOST_UPSTREAM_LOC, mostUpstreamWritableIndex);
    }
}