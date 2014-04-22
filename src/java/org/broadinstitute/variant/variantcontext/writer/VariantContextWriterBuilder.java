/*
* Copyright (c) 2014 The Broad Institute
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
import net.sf.samtools.util.Md5CalculatingOutputStream;
import net.sf.samtools.util.RuntimeIOException;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.index.IndexCreator;
import org.broad.tribble.index.tabix.TabixFormat;
import org.broad.tribble.index.tabix.TabixIndexCreator;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.EnumSet;

/**
 * Created with IntelliJ IDEA.
 * User: thibault
 * Date: 3/7/14
 * Time: 2:07 PM
 *
 * Provides methods for creating VariantContextWriters using the Builder pattern.
 * Replaces VariantContextWriterFactory.
 *
 * The caller must choose an output file or an output stream for the VariantContextWriter to write to.
 * When a file is chosen, the output stream is created implicitly based on Defaults and options passed to the builder.
 * When a stream is chosen, it is passed unchanged to the VariantContextWriter.
 *
 * Example: Create a series of files with buffering and indexing on the fly.
 * Determine the appropriate file type based on filename.
 *
 *  VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
 *      .setReferenceDictionary(refDict)
 *      .setOption(Options.INDEX_ON_THE_FLY)
 *      .setBuffer(8192);
 *
 *  VariantContextWriter sample1_writer = builder
 *      .setOutputFile("sample1.vcf")
 *      .build();
 *  VariantContextWriter sample2_writer = builder
 *      .setOutputFile("sample2.bcf")
 *      .build();
 *  VariantContextWriter sample3_writer = builder
 *      .setOutputFile("sample3.vcf.bgzf")
 *      .build();
 *
 * Example: Explicitly turn off buffering and explicitly set the file type
 *
 *  VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
 *      .setReferenceDictionary(refDict)
 *      .setOption(Options.INDEX_ON_THE_FLY)
 *      .unsetBuffering();
 *
 *  VariantContextWriter sample1_writer = builder
 *      .setOutputFile("sample1.custom_extension")
 *      .setOutputFileType(OutputType.VCF)
 *      .build();
 *  VariantContextWriter sample2_writer = builder
 *      .setOutputFile("sample2.custom_extension")
 *      .setOutputFileType(OutputType.BLOCK_COMPRESSED_VCF)
 *      .build();

 */
public class VariantContextWriterBuilder {
    public static final EnumSet<Options> DEFAULT_OPTIONS = EnumSet.of(Options.INDEX_ON_THE_FLY);
    public static final EnumSet<Options> NO_OPTIONS = EnumSet.noneOf(Options.class);

    public enum OutputType {
        UNSPECIFIED,
        VCF,
        BCF,
        BLOCK_COMPRESSED_VCF,
        VCF_STREAM,
        BCF_STREAM
    }

    public static final EnumSet<OutputType> FILE_TYPES = EnumSet.of(OutputType.VCF, OutputType.BCF, OutputType.BLOCK_COMPRESSED_VCF);
    public static final EnumSet<OutputType> STREAM_TYPES = EnumSet.of(OutputType.VCF_STREAM, OutputType.BCF_STREAM);

    private SAMSequenceDictionary refDict = null;
    private OutputType outType = OutputType.UNSPECIFIED;
    private File outFile = null;
    private OutputStream outStream = null;
    private IndexCreator idxCreator = null;
    private int bufferSize = Defaults.BUFFER_SIZE;
    private boolean createMD5 = Defaults.CREATE_MD5;
    private EnumSet<Options> options = DEFAULT_OPTIONS.clone();

    /**
     * Default constructor.  Adds USE_ASYNC_IO to the Options if it is present in Defaults.
     */
    public VariantContextWriterBuilder() {
        if (Defaults.USE_ASYNC_IO)
            options.add(Options.USE_ASYNC_IO);
    }

    /**
     * Set the reference dictionary to be used by VariantContextWriters created by this builder
     *
     * @param refDict the reference dictionary
     * @return this VariantContextWriterBuilder
     */
    public VariantContextWriterBuilder setReferenceDictionary(final SAMSequenceDictionary refDict) {
        this.refDict = refDict;
        return this;
    }

    /**
     * Set the output file for the next VariantContextWriter created by this builder
     * Determines file type implicitly from the filename
     *
     * @param outFile the file the VariantContextWriter will write to
     * @return this VariantContextWriterBuilder
     */
    public VariantContextWriterBuilder setOutputFile(final File outFile) {
        this.outFile = outFile;
        this.outStream = null;
        determineOutputTypeFromFilename();
        return this;
    }

    /**
     * Set the output file for the next VariantContextWriter created by this builder
     * Determines file type implicitly from the filename
     *
     * @param outFile the file the VariantContextWriter will write to
     * @return this VariantContextWriterBuilder
     */
    public VariantContextWriterBuilder setOutputFile(final String outFile) {
        this.outFile = new File(outFile);
        this.outStream = null;
        determineOutputTypeFromFilename();
        return this;
    }

    /**
     * Set the output file type for the next VariantContextWriter created by this builder
     *
     * @param outType the type of file the VariantContextWriter will write to
     * @return this VariantContextWriterBuilder
     */
    public VariantContextWriterBuilder setOutputFileType(final OutputType outType) {
        if (!FILE_TYPES.contains(outType))
            throw new IllegalArgumentException("Must choose a file type, not other output types.");

        if (this.outFile == null || this.outStream != null)
            throw new IllegalArgumentException("Cannot set a file type if the output is not to a file.");

        this.outType = outType;
        return this;
    }

    /**
     * Set the output VCF stream for the next VariantContextWriter created by this builder
     * If buffered writing is desired, caller must provide some kind of buffered OutputStream.
     *
     * @param outStream the output stream to write to
     * @return this VariantContextWriterBuilder
     */
    public VariantContextWriterBuilder setOutputVCFStream(final OutputStream outStream) {
        this.outStream = outStream;
        this.outFile = null;
        this.outType = OutputType.VCF_STREAM;
        return this;
    }

    /**
     * Set the output BCF stream for the next VariantContextWriter created by this builder
     * If buffered writing is desired, caller must provide some kind of buffered OutputStream.
     *
     * @param outStream the output stream to write to
     * @return this VariantContextWriterBuilder
     */
    public VariantContextWriterBuilder setOutputBCFStream(final OutputStream outStream) {
        this.outStream = outStream;
        this.outFile = null;
        this.outType = OutputType.BCF_STREAM;
        return this;
    }

    /**
     * Set the output stream (VCF, by default) for the next VariantContextWriter created by this builder
     * If buffered writing is desired, caller must provide some kind of buffered OutputStream.
     *
     * @param outStream the output stream to write to
     * @return this VariantContextWriterBuilder
     */
    public VariantContextWriterBuilder setOutputStream(final OutputStream outStream) {
        return setOutputVCFStream(outStream);
    }

    /**
     * Set an IndexCreator for the next VariantContextWriter created by this builder
     *
     * @param idxCreator the IndexCreator to use
     * @return this VariantContextWriterBuilder
     */
    public VariantContextWriterBuilder setIndexCreator(final IndexCreator idxCreator) {
        this.idxCreator = idxCreator;
        return this;
    }

    /**
     * Do not pass an IndexCreator to the next VariantContextWriter created by this builder
     *
     * @return this VariantContextWriterBuilder
     */
    public VariantContextWriterBuilder clearIndexCreator() {
        this.idxCreator = null;
        return this;
    }

    /**
     * Set a buffer size for the file output stream passed to the next VariantContextWriter created by this builder
     * Set to 0 for no buffering
     * Does not affect OutputStreams passed directly to VariantContextWriterBuilder
     *
     * @param bufferSize the buffer size to use
     * @return this VariantContextWriterBuilder
     */
    public VariantContextWriterBuilder setBuffer(final int bufferSize) {
        this.bufferSize = bufferSize;
        return this;
    }

    /**
     * Do not use buffering in the next VariantContextWriter created by this builder
     * Does not affect OutputStreams passed directly to VariantContextWriterBuilder
     *
     * @return this VariantContextWriterBuilder
     */
    public VariantContextWriterBuilder unsetBuffering() {
        this.bufferSize = 0;
        return this;
    }

    /**
     * Choose whether to also create an MD5 digest file for the next VariantContextWriter created by this builder
     *
     * @param createMD5 boolean, true to create an MD5 digest
     * @return this VariantContextWriterBuilder
     */
    public VariantContextWriterBuilder setCreateMD5(final boolean createMD5) {
        this.createMD5 = createMD5;
        return this;
    }

    /**
     * Create an MD5 digest file for the next VariantContextWriter created by this builder
     *
     * @return this VariantContextWriterBuilder
     */
    public VariantContextWriterBuilder setCreateMD5() {
        return setCreateMD5(true);
    }

    /**
     * Don't create an MD5 digest file for the next VariantContextWriter created by this builder
     *
     * @return this VariantContextWriterBuilder
     */
    public VariantContextWriterBuilder unsetCreateMD5() {
        return setCreateMD5(false);
    }

    /**
     * Replace the set of Options for the VariantContextWriterBuilder with a new set
     *
     * @param options the complete set of options to use
     * @return this VariantContextWriterBuilder
     */
    public VariantContextWriterBuilder setOptions(final EnumSet<Options> options) {
        this.options = options;
        return this;
    }

    /**
     * Add one option to the set of Options for the VariantContextWriterBuilder, if it's not already present
     *
     * @param option the option to set
     * @return this VariantContextWriterBuilder
     */
    public VariantContextWriterBuilder setOption(final Options option) {
        this.options.add(option);
        return this;
    }

    /**
     * Remove one option from the set of Options for the VariantContextWriterBuilder, if it's present
     *
     * @param option the option to unset
     * @return this VariantContextWriterBuilder
     */
    public VariantContextWriterBuilder unsetOption(final Options option) {
        this.options.remove(option);
        return this;
    }

    /**
     * Remove all options from the set of Options for the VariantContextWriterBuilder
     *
     * @return this VariantContextWriterBuilder
     */
    public VariantContextWriterBuilder clearOptions() {
        this.options = NO_OPTIONS;
        return this;
    }

    /**
     * Validate and build the VariantContextWriter
     *
     * @return the VariantContextWriter as specified by previous method calls
     */
    public VariantContextWriter build() {
        VariantContextWriter writer = null;

        // don't allow FORCE_BCF to modify the outType state
        OutputType typeToBuild = this.outType;

        if (this.options.contains(Options.FORCE_BCF)) {
            if (FILE_TYPES.contains(this.outType))
                typeToBuild = OutputType.BCF;
            else if (STREAM_TYPES.contains(this.outType))
                typeToBuild = OutputType.BCF_STREAM;
        }

        OutputStream outStreamFromFile = this.outStream;
        if (FILE_TYPES.contains(this.outType)) {
            try {
                outStreamFromFile = IOUtil.maybeBufferOutputStream(new FileOutputStream(outFile), bufferSize);
            } catch (final FileNotFoundException e) {
                throw new RuntimeIOException("File not found: " + outFile, e);
            }

            if (createMD5)
                outStreamFromFile = new Md5CalculatingOutputStream(outStreamFromFile, new File(outFile.getAbsolutePath() + ".md5"));
        }

        switch (typeToBuild) {
            case UNSPECIFIED:
                throw new IllegalArgumentException("Must specify file or stream output type.");
            case VCF:
                if ((refDict == null) && (options.contains(Options.INDEX_ON_THE_FLY)))
                    throw new IllegalArgumentException("A reference dictionary is required for creating Tribble indices on the fly");

                writer = createVCFWriter(outFile, outStreamFromFile);
                break;
            case BLOCK_COMPRESSED_VCF:
                if (refDict == null)
                    idxCreator = new TabixIndexCreator(TabixFormat.VCF);
                else
                    idxCreator = new TabixIndexCreator(refDict, TabixFormat.VCF);

                writer = createVCFWriter(outFile, new BlockCompressedOutputStream(outStreamFromFile, outFile));
                break;
            case BCF:
                if ((refDict == null) && (options.contains(Options.INDEX_ON_THE_FLY)))
                    throw new IllegalArgumentException("A reference dictionary is required for creating Tribble indices on the fly");

                writer = createBCFWriter(outFile, outStreamFromFile);
                break;
            case VCF_STREAM:
                if (options.contains(Options.INDEX_ON_THE_FLY))
                    throw new IllegalArgumentException("VCF index creation not supported for stream output.");

                writer = createVCFWriter(null, outStream);
                break;
            case BCF_STREAM:
                if (options.contains(Options.INDEX_ON_THE_FLY))
                    throw new IllegalArgumentException("BCF index creation not supported for stream output.");

                writer = createBCFWriter(null, outStream);
                break;
        }

        if (this.options.contains(Options.USE_ASYNC_IO))
            writer = new AsyncVariantContextWriter(writer, AsyncVariantContextWriter.DEFAULT_QUEUE_SIZE);

        return writer;
     }

    private void determineOutputTypeFromFilename() {
        if (isBCF(this.outFile)) {
            this.outType = OutputType.BCF;
        } else if (isCompressedVCF(this.outFile)) {
            this.outType = OutputType.BLOCK_COMPRESSED_VCF;
        } else if (isVCF(this.outFile)) {
            this.outType = OutputType.VCF;
        }
        else {
            this.outType = OutputType.UNSPECIFIED;
        }
    }

    private boolean isVCF(final File outFile) {
        return outFile != null && outFile.getName().endsWith(".vcf");
    }

    private boolean isBCF(final File outFile) {
        return outFile != null && outFile.getName().endsWith(".bcf");
    }

    private boolean isCompressedVCF(final File outFile) {
        if (outFile == null)
            return false;

        return AbstractFeatureReader.hasBlockCompressedExtension(outFile);
    }

    private VariantContextWriter createVCFWriter(final File writerFile, final OutputStream writerStream) {
        if (idxCreator == null) {
            return new VCFWriter(writerFile, writerStream, refDict,
                    options.contains(Options.INDEX_ON_THE_FLY),
                    options.contains(Options.DO_NOT_WRITE_GENOTYPES),
                    options.contains(Options.ALLOW_MISSING_FIELDS_IN_HEADER));
        }
        else {
            return new VCFWriter(writerFile, writerStream, refDict, idxCreator,
                    options.contains(Options.INDEX_ON_THE_FLY),
                    options.contains(Options.DO_NOT_WRITE_GENOTYPES),
                    options.contains(Options.ALLOW_MISSING_FIELDS_IN_HEADER));
        }
    }

    private VariantContextWriter createBCFWriter(final File writerFile, final OutputStream writerStream) {
        if (idxCreator == null) {
            return new BCF2Writer(writerFile, writerStream, refDict,
                    options.contains(Options.INDEX_ON_THE_FLY),
                    options.contains(Options.DO_NOT_WRITE_GENOTYPES));
        }
        else {
            return new BCF2Writer(writerFile, writerStream, refDict, idxCreator,
                    options.contains(Options.INDEX_ON_THE_FLY),
                    options.contains(Options.DO_NOT_WRITE_GENOTYPES));
        }
    }
}
