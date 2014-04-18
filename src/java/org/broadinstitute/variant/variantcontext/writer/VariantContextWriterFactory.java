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
import net.sf.samtools.util.BlockCompressedOutputStream;
import net.sf.samtools.util.IOUtil;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.index.IndexCreator;
import org.broad.tribble.index.tabix.TabixFormat;
import org.broad.tribble.index.tabix.TabixIndexCreator;

import java.io.*;
import java.util.EnumSet;

/**
 * Factory methods to create VariantContext writers
 *
 * @author depristo
 * @since 5/12
 *
 * @deprecated Replaced by {@link org.broadinstitute.variant.variantcontext.writer.VariantContextWriterBuilder}
 */
@Deprecated
public class VariantContextWriterFactory {

    public static final EnumSet<Options> DEFAULT_OPTIONS = EnumSet.of(Options.INDEX_ON_THE_FLY);
    public static final EnumSet<Options> NO_OPTIONS = EnumSet.noneOf(Options.class);

    static {
        if (Defaults.USE_ASYNC_IO) {
            DEFAULT_OPTIONS.add(Options.USE_ASYNC_IO);
        }
    }

    private VariantContextWriterFactory() {}

    public static VariantContextWriter create(final File location, final SAMSequenceDictionary refDict) {
        return create(location, openOutputStream(location), refDict, DEFAULT_OPTIONS);
    }

    public static VariantContextWriter create(final File location, final SAMSequenceDictionary refDict, final EnumSet<Options> options) {
        return create(location, openOutputStream(location), refDict, options);
    }

    /**
     * @param output If buffered writing is desired, caller must provide some kind of buffered OutputStream.
     */
    public static VariantContextWriter create(final File location,
                                              final OutputStream output,
                                              final SAMSequenceDictionary refDict) {
        return create(location, output, refDict, DEFAULT_OPTIONS);
    }

    /**
     * @param output If buffered writing is desired, caller must provide some kind of buffered OutputStream.
     */
    public static VariantContextWriter create(final OutputStream output,
                                              final SAMSequenceDictionary refDict,
                                              final EnumSet<Options> options) {
        return create(null, output, refDict, options);
    }

    /**
     * @param location Note that this parameter is used to producing intelligent log messages, and for naming the index,
     *                 but does not control where the file is written
     * @param output This is where the BCF is actually written. If buffered writing is desired, caller must provide
     *               some kind of buffered OutputStream.
     */
    public static VariantContextWriter createBcf2(final File location,
                                                  final OutputStream output,
                                                  final SAMSequenceDictionary refDict,
                                                  final EnumSet<Options> options) {
        return maybeWrapWithAsyncWriter(new BCF2Writer(location, output, refDict,
                options.contains(Options.INDEX_ON_THE_FLY),
                options.contains(Options.DO_NOT_WRITE_GENOTYPES)), options);
    }

    /**
     * @param location Note that this parameter is used to producing intelligent log messages, and for naming the index,
     *                 but does not control where the file is written
     * @param output This is where the BCF is actually written.  If buffered writing is desired, caller must provide
     *               some kind of buffered OutputStream.
     */
    public static VariantContextWriter createBcf2(final File location,
                                                  final OutputStream output,
                                                  final SAMSequenceDictionary refDict,
                                                  final IndexCreator indexCreator,
                                                  final EnumSet<Options> options) {
        return maybeWrapWithAsyncWriter(new BCF2Writer(location, output, refDict, indexCreator,
                options.contains(Options.INDEX_ON_THE_FLY),
                options.contains(Options.DO_NOT_WRITE_GENOTYPES)), options);
    }

    /**
     * @param location Note that this parameter is used to producing intelligent log messages, and for naming the index,
     *                 but does not control where the file is written
     * @param output This is where the VCF is actually written. If buffered writing is desired, caller must provide
     *               some kind of buffered OutputStream.
     */
    public static VariantContextWriter createVcf(final File location,
                                                 final OutputStream output,
                                                 final SAMSequenceDictionary refDict,
                                                 final EnumSet<Options> options) {
        return maybeWrapWithAsyncWriter(new VCFWriter(location, output, refDict,
                options.contains(Options.INDEX_ON_THE_FLY),
                options.contains(Options.DO_NOT_WRITE_GENOTYPES),
                options.contains(Options.ALLOW_MISSING_FIELDS_IN_HEADER)), options);
    }

    /**
     * @param location Note that this parameter is used to producing intelligent log messages, and for naming the index,
     *                 but does not control where the file is written
     * @param output This is where the VCF is actually written.  If buffered writing is desired, caller must provide
     *               some kind of buffered OutputStream.
     */
    public static VariantContextWriter createVcf(final File location,
                                                 final OutputStream output,
                                                 final SAMSequenceDictionary refDict,
                                                 final IndexCreator indexCreator,
                                                 final EnumSet<Options> options) {
        return maybeWrapWithAsyncWriter(new VCFWriter(location, output, refDict, indexCreator,
                options.contains(Options.INDEX_ON_THE_FLY),
                options.contains(Options.DO_NOT_WRITE_GENOTYPES),
                options.contains(Options.ALLOW_MISSING_FIELDS_IN_HEADER)), options);
    }

    /**
     * @param location Note that this parameter is used to producing intelligent log messages,
     *                 but does not control where the file is written
     * @param output This is where the VCF is actually written.  If buffered writing is desired, caller must provide
     *               some kind of buffered OutputStream.
     */
    public static VariantContextWriter createBlockCompressedVcf(final File location,
                                                                final OutputStream output,
                                                                final SAMSequenceDictionary refDict,
                                                                final EnumSet<Options> options) {
        final TabixIndexCreator indexCreator;
        if (options.contains(Options.INDEX_ON_THE_FLY)) {
            indexCreator = new TabixIndexCreator(refDict, TabixFormat.VCF);
        } else {
            indexCreator = null;
        }
        return maybeWrapWithAsyncWriter(new VCFWriter(location, BlockCompressedOutputStream.maybeBgzfWrapOutputStream(location, output),
                refDict, indexCreator,
                options.contains(Options.INDEX_ON_THE_FLY),
                options.contains(Options.DO_NOT_WRITE_GENOTYPES),
                options.contains(Options.ALLOW_MISSING_FIELDS_IN_HEADER)), options);
    }

    /**
     * @param location Note that this parameter is used to producing intelligent log messages,
     *                 but does not control where the file is written
     * @param output This is where the VCF is actually written. If buffered writing is desired, caller must provide
     *               some kind of buffered OutputStream.
     */
    public static VariantContextWriter createBlockCompressedVcf(final File location,
                                                                final OutputStream output,
                                                                final SAMSequenceDictionary refDict,
                                                                final IndexCreator indexCreator,
                                                                final EnumSet<Options> options) {
        return maybeWrapWithAsyncWriter(new VCFWriter(location, BlockCompressedOutputStream.maybeBgzfWrapOutputStream(location, output),
                refDict, indexCreator,
                options.contains(Options.INDEX_ON_THE_FLY),
                options.contains(Options.DO_NOT_WRITE_GENOTYPES),
                options.contains(Options.ALLOW_MISSING_FIELDS_IN_HEADER)), options);
    }

    public static VariantContextWriter create(final File location,
        final OutputStream output,
        final SAMSequenceDictionary refDict,
        final EnumSet<Options> options) {

        if (isBCFOutput(location, options)) {
            return createBcf2(location, output, refDict, options);
        } else if (isCompressedVcf(location)) {
            return createBlockCompressedVcf(location, output, refDict, options);
        } else {
            return createVcf(location, output, refDict, options);
        }
    }

    /**
     * @param output If buffered writing is desired, caller must provide some kind of buffered OutputStream.
     */
    public static VariantContextWriter create(final File location,
                                              final OutputStream output,
                                              final SAMSequenceDictionary refDict,
                                              final IndexCreator indexCreator,
                                              final EnumSet<Options> options) {

        if (isBCFOutput(location, options)) {
            return createBcf2(location, output, refDict, indexCreator, options);
        } else if (isCompressedVcf(location)) {
            return createBlockCompressedVcf(location, output, refDict, indexCreator, options);
        } else {
            return createVcf(location, output, refDict, indexCreator, options);
        }
    }

    private static VariantContextWriter maybeWrapWithAsyncWriter(final VariantContextWriter writer,
                                                                 final EnumSet<Options> options) {
        if (options.contains(Options.USE_ASYNC_IO)) {
            return new AsyncVariantContextWriter(writer, AsyncVariantContextWriter.DEFAULT_QUEUE_SIZE);
        }
        else return writer;
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

    public static boolean isCompressedVcf(final File location) {
        if (location == null)
            return false;

        return AbstractFeatureReader.hasBlockCompressedExtension(location);
    }

    public static VariantContextWriter sortOnTheFly(final VariantContextWriter innerWriter, final int maxCachingStartDistance) {
        return sortOnTheFly(innerWriter, maxCachingStartDistance, false);
    }

    public static VariantContextWriter sortOnTheFly(final VariantContextWriter innerWriter, final int maxCachingStartDistance, final boolean takeOwnershipOfInner) {
        return new SortingVariantContextWriter(innerWriter, maxCachingStartDistance, takeOwnershipOfInner);
    }

    /**
     * Returns a output stream writing to location, or throws an exception if this fails
     * @param location
     * @return
     */
    protected static OutputStream openOutputStream(final File location) {
        try {
            return IOUtil.maybeBufferOutputStream(new FileOutputStream(location));
        } catch (final FileNotFoundException e) {
            throw new RuntimeException(location + ": Unable to create VCF writer", e);
        }
    }
}
