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

import net.sf.samtools.Defaults;
import net.sf.samtools.SAMSequenceDictionary;

import java.io.*;
import java.util.EnumSet;

/**
 * Factory methods to create VariantContext writers
 *
 * @author depristo
 * @since 5/12
 */
public class VariantContextWriterFactory {

    public static final EnumSet<Options> DEFAULT_OPTIONS = EnumSet.of(Options.INDEX_ON_THE_FLY);
    public static final EnumSet<Options> NO_OPTIONS = EnumSet.noneOf(Options.class);

    private VariantContextWriterFactory() {}

    public static VariantContextWriter create(final File location, final SAMSequenceDictionary refDict) {
        return create(location, openOutputStream(location), refDict, DEFAULT_OPTIONS);
    }

    public static VariantContextWriter create(final File location, final SAMSequenceDictionary refDict, final EnumSet<Options> options) {
        return create(location, openOutputStream(location), refDict, options);
    }

    public static VariantContextWriter create(final File location,
                                              final OutputStream output,
                                              final SAMSequenceDictionary refDict) {
        return create(location, output, refDict, DEFAULT_OPTIONS);
    }

    public static VariantContextWriter create(final OutputStream output,
                                              final SAMSequenceDictionary refDict,
                                              final EnumSet<Options> options) {
        return create(null, output, refDict, options);
    }

    public static VariantContextWriter create(final File location,
                                              final OutputStream output,
                                              final SAMSequenceDictionary refDict,
                                              final EnumSet<Options> options) {
        final boolean enableBCF = isBCFOutput(location, options);

        if ( enableBCF )
            return new BCF2Writer(location, output, refDict,
                    options.contains(Options.INDEX_ON_THE_FLY),
                    options.contains(Options.DO_NOT_WRITE_GENOTYPES));
        else {
            return new VCFWriter(location, output, refDict,
                    options.contains(Options.INDEX_ON_THE_FLY),
                    options.contains(Options.DO_NOT_WRITE_GENOTYPES),
                    options.contains(Options.ALLOW_MISSING_FIELDS_IN_HEADER));
        }
    }

    /**
     * Should we output a BCF file based solely on the name of the file at location?
     *
     * @param location
     * @return
     */
    public static boolean isBCFOutput(final File location) {
        return isBCFOutput(location, EnumSet.noneOf(Options.class));
    }

    public static boolean isBCFOutput(final File location, final EnumSet<Options> options) {
        return options.contains(Options.FORCE_BCF) || (location != null && location.getName().contains(".bcf"));
    }

    public static VariantContextWriter sortOnTheFly(final VariantContextWriter innerWriter, int maxCachingStartDistance) {
        return sortOnTheFly(innerWriter, maxCachingStartDistance, false);
    }

    public static VariantContextWriter sortOnTheFly(final VariantContextWriter innerWriter, int maxCachingStartDistance, boolean takeOwnershipOfInner) {
        return new SortingVariantContextWriter(innerWriter, maxCachingStartDistance, takeOwnershipOfInner);
    }

    /**
     * Returns a output stream writing to location, or throws an exception if this fails
     * @param location
     * @return
     */
    protected static OutputStream openOutputStream(final File location) {
        try {
            return new BufferedOutputStream(new FileOutputStream(location), Defaults.BUFFER_SIZE);
        } catch (FileNotFoundException e) {
            throw new RuntimeException(location + ": Unable to create VCF writer", e);
        }
    }
}
