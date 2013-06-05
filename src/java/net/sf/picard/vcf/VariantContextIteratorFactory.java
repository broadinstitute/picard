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

import net.sf.picard.io.IoUtil;

import java.io.File;
import java.io.InputStream;

/**
 * Creates an iterator for a VCF/BCF based on the filename
 */
public class VariantContextIteratorFactory {
    private VariantContextIteratorFactory() {}

    public static VariantContextIterator create(final File location) {
        final InputStream inputStream = IoUtil.openFileForReading(location);
        // TODO: Both this and VariantContextWriterFactory base this on filename, in the future we may want to change this
        if (location.getName().toLowerCase().endsWith(".bcf")) {
            return new BcfIterator(inputStream);
        } else {
            return new VcfIterator(inputStream);
        }
    }
}

