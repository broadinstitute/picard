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

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.LocationAware;
import org.broad.tribble.index.*;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFHeader;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

/**
 * this class writes VCF files
 */
abstract class IndexingVariantContextWriter implements VariantContextWriter {
    private final String name;
    private final File location;
    private final SAMSequenceDictionary refDict;

    private OutputStream outputStream;
    private LocationAware locationSource = null;
    private IndexCreator indexer = null;

    private IndexingVariantContextWriter(final String name, final File location, final OutputStream output, final SAMSequenceDictionary refDict) {
        this.name = name;
        this.location = location;
        this.outputStream = output;
        this.refDict = refDict;
    }

    /**
     * Create a VariantContextWriter with an associated index using the default index creator
     *
     * @param name  the name of this writer (i.e. the file name or stream)
     * @param location  the path to the output file
     * @param output    the output stream to write to
     * @param refDict   the reference dictionary
     * @param enableOnTheFlyIndexing    is OTF indexing enabled?
     */
    @Requires({"name != null",
            "! ( location == null && output == null )",
            "! ( enableOnTheFlyIndexing && location == null )"})
    protected IndexingVariantContextWriter(final String name, final File location, final OutputStream output, final SAMSequenceDictionary refDict,
                                           final boolean enableOnTheFlyIndexing) {
        this(name, location, output, refDict);

        if ( enableOnTheFlyIndexing ) {
            initIndexingWriter(new DynamicIndexCreator(location, IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME));
        }
    }

    /**
     * Create a VariantContextWriter with an associated index using a custom index creator
     *
     * @param name  the name of this writer (i.e. the file name or stream)
     * @param location  the path to the output file
     * @param output    the output stream to write to
     * @param refDict   the reference dictionary
     * @param enableOnTheFlyIndexing    is OTF indexing enabled?
     * @param idxCreator    the custom index creator.  NOTE: must be initialized
     */
    @Requires({"name != null",
            "! ( location == null && output == null )",
            "! ( enableOnTheFlyIndexing && location == null )",
            "idxCreator != null"})
    protected IndexingVariantContextWriter(final String name, final File location, final OutputStream output, final SAMSequenceDictionary refDict,
                                           final boolean enableOnTheFlyIndexing, final IndexCreator idxCreator) {
        this(name, location, output, refDict);

        if ( enableOnTheFlyIndexing ) {
            // TODO: Handle non-Tribble IndexCreators
            initIndexingWriter(idxCreator);
        }
    }

    @Requires({"idxCreator != null"})
    private void initIndexingWriter(final IndexCreator idxCreator) {
        indexer = idxCreator;
        if (outputStream instanceof LocationAware) {
            locationSource = (LocationAware)outputStream;
        } else {
            final PositionalOutputStream positionalOutputStream = new PositionalOutputStream(outputStream);
            locationSource = positionalOutputStream;
            outputStream = positionalOutputStream;
        }
    }

    @Ensures("result != null")
    public OutputStream getOutputStream() {
        return outputStream;
    }

    @Ensures("result != null")
    public String getStreamName() {
        return name;
    }

    public abstract void writeHeader(VCFHeader header);

    /**
     * attempt to close the VCF file
     */
    public void close() {
        try {
            // close the underlying output stream
            outputStream.close();

            // close the index stream (keep it separate to help debugging efforts)
            if (indexer != null) {
                if (indexer instanceof TribbleIndexCreator) {
                    setIndexSequenceDictionary((TribbleIndexCreator)indexer, refDict);
                }
                final Index index = indexer.finalizeIndex(locationSource.getPosition());
                index.writeBasedOnFeatureFile(location);
            }


        } catch (final IOException e) {
            throw new RuntimeException("Unable to close index for " + getStreamName(), e);
        }
    }

    /**
     * @return the reference sequence dictionary used for the variant contexts being written
     */
    public SAMSequenceDictionary getRefDict() {
        return refDict;
    }

    /**
     * add a record to the file
     *
     * @param vc      the Variant Context object
     */
    public void add(final VariantContext vc) {
        // if we are doing on the fly indexing, add the record ***before*** we write any bytes
        if ( indexer != null )
            indexer.addFeature(vc, locationSource.getPosition());
    }

    /**
     * Returns a reasonable "name" for this writer, to display to the user if something goes wrong
     *
     * @param location
     * @param stream
     * @return
     */
    protected static final String writerName(final File location, final OutputStream stream) {
        return location == null ? stream.toString() : location.getAbsolutePath();
    }

    // a constant we use for marking sequence dictionary entries in the Tribble index property list
    private static final String SequenceDictionaryPropertyPredicate = "DICT:";

    private static void setIndexSequenceDictionary(final TribbleIndexCreator indexCreator, final SAMSequenceDictionary dict) {
        for (final SAMSequenceRecord seq : dict.getSequences()) {
            final String contig = SequenceDictionaryPropertyPredicate + seq.getSequenceName();
            final String length = String.valueOf(seq.getSequenceLength());
            indexCreator.addProperty(contig,length);
        }
    }
}

/**
 * Wraps output stream in a manner which keeps track of the position within the file and allowing writes
 * at arbitrary points
 */
final class PositionalOutputStream extends OutputStream implements LocationAware
{
    private final OutputStream out;
    private long position = 0;

    public PositionalOutputStream(final OutputStream out) {
        this.out = out;
    }

    public final void write(final byte[] bytes) throws IOException {
        write(bytes, 0, bytes.length);
    }

    public final void write(final byte[] bytes, final int startIndex, final int numBytes) throws IOException {
        position += numBytes;
        out.write(bytes, startIndex, numBytes);
    }

    public final void write(final int c)  throws IOException {
        position++;
        out.write(c);
    }

    public final long getPosition() { return position; }

    @Override
    public void close() throws IOException {
        super.close();
        out.close();
    }
}