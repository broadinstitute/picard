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

package org.broadinstitute.variant.bcf2;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.FeatureCodecHeader;
import org.broad.tribble.TribbleException;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.variant.utils.GeneralUtils;
import org.broadinstitute.variant.vcf.*;
import org.broadinstitute.variant.variantcontext.*;

import java.io.ByteArrayInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Decode BCF2 files
 */
public final class BCF2Codec implements FeatureCodec<VariantContext> {
    private final static int ALLOWED_MAJOR_VERSION = 2;
    private final static int MIN_MINOR_VERSION = 1;

    private BCFVersion bcfVersion = null;

    private VCFHeader header = null;

    /**
     * Maps offsets (encoded in BCF) into contig names (from header) for the CHROM field
     */
    private final ArrayList<String> contigNames = new ArrayList<String>();

    /**
     * Maps header string names (encoded in VCF) into strings found in the BCF header
     *
     * Initialized when processing the header
     */
    private ArrayList<String> dictionary;

    /**
     * Our decoder that reads low-level objects from the BCF2 records
     */
    private final BCF2Decoder decoder = new BCF2Decoder();

    /**
     * Provides some sanity checking on the header
     */
    private final static int MAX_HEADER_SIZE = 0x08000000;

    /**
     * Genotype field decoders that are initialized when the header is read
     */
    private BCF2GenotypeFieldDecoders gtFieldDecoders = null;

    /**
     * A cached array of GenotypeBuilders for efficient genotype decoding.
     *
     * Caching it allows us to avoid recreating this intermediate data
     * structure each time we decode genotypes
     */
    private GenotypeBuilder[] builders = null;

    // for error handling
    private int recordNo = 0;
    private int pos = 0;


    // ----------------------------------------------------------------------
    //
    // Feature codec interface functions
    //
    // ----------------------------------------------------------------------

    @Override
    public Feature decodeLoc( final PositionalBufferedStream inputStream ) {
        return decode(inputStream);
    }

    @Override
    public VariantContext decode( final PositionalBufferedStream inputStream ) {
        try {
            recordNo++;
            final VariantContextBuilder builder = new VariantContextBuilder();

            final int sitesBlockSize = decoder.readBlockSize(inputStream);
            final int genotypeBlockSize = decoder.readBlockSize(inputStream);

            decoder.readNextBlock(sitesBlockSize, inputStream);
            decodeSiteLoc(builder);
            final SitesInfoForDecoding info = decodeSitesExtendedInfo(builder);

            decoder.readNextBlock(genotypeBlockSize, inputStream);
            createLazyGenotypesDecoder(info, builder);
            return builder.fullyDecoded(true).make();
        } catch ( IOException e ) {
            throw new TribbleException("Failed to read BCF file", e);
        }
    }

    @Override
    public Class<VariantContext> getFeatureType() {
        return VariantContext.class;
    }

    @Override
    public FeatureCodecHeader readHeader( final PositionalBufferedStream inputStream ) {
        try {
            // note that this reads the magic as well, and so does double duty
            bcfVersion = BCFVersion.readBCFVersion(inputStream);
            if ( bcfVersion == null )
                error("Input stream does not contain a BCF encoded file; BCF magic header info not found");

            if ( bcfVersion.getMajorVersion() != ALLOWED_MAJOR_VERSION )
                error("BCF2Codec can only process BCF2 files, this file has major version " + bcfVersion.getMajorVersion());
            if ( bcfVersion.getMinorVersion() < MIN_MINOR_VERSION )
                error("BCF2Codec can only process BCF2 files with minor version >= " + MIN_MINOR_VERSION + " but this file has minor version " + bcfVersion.getMinorVersion());

            if ( GeneralUtils.DEBUG_MODE_ENABLED ) {
                System.err.println("Parsing data stream with BCF version " + bcfVersion);
            }

            final int headerSizeInBytes = BCF2Type.INT32.read(inputStream);

            if ( headerSizeInBytes <= 0 || headerSizeInBytes > MAX_HEADER_SIZE) // no bigger than 8 MB
                error("BCF2 header has invalid length: " + headerSizeInBytes + " must be >= 0 and < "+ MAX_HEADER_SIZE);

            final byte[] headerBytes = new byte[headerSizeInBytes];
            if ( inputStream.read(headerBytes) != headerSizeInBytes )
                error("Couldn't read all of the bytes specified in the header length = " + headerSizeInBytes);

            final PositionalBufferedStream bps = new PositionalBufferedStream(new ByteArrayInputStream(headerBytes));
            final AsciiLineReader headerReader = new AsciiLineReader(bps);
            final VCFCodec headerParser = new VCFCodec();
            this.header = (VCFHeader)headerParser.readHeader(headerReader);
            bps.close();
        } catch ( IOException e ) {
            throw new TribbleException("I/O error while reading BCF2 header");
        }

        // create the config offsets
        if ( ! header.getContigLines().isEmpty() ) {
            contigNames.clear();
            for ( final VCFContigHeaderLine contig : header.getContigLines()) {
                if ( contig.getID() == null || contig.getID().equals("") )
                    error("found a contig with an invalid ID " + contig);
                contigNames.add(contig.getID());
            }
        } else {
            error("Didn't find any contig lines in BCF2 file header");
        }

        // create the string dictionary
        dictionary = parseDictionary(header);

        // prepare the genotype field decoders
        gtFieldDecoders = new BCF2GenotypeFieldDecoders(header);

        // create and initialize the genotype builder array
        final int nSamples = header.getNGenotypeSamples();
        builders = new GenotypeBuilder[nSamples];
        for ( int i = 0; i < nSamples; i++ ) {
            builders[i] = new GenotypeBuilder(header.getGenotypeSamples().get(i));
        }

        // position right before next line (would be right before first real record byte at end of header)
        return new FeatureCodecHeader(header, inputStream.getPosition());
    }

    @Override
    public boolean canDecode( final String path ) {
        FileInputStream fis = null;
        try {
            fis = new FileInputStream(path);
            final BCFVersion version = BCFVersion.readBCFVersion(fis);
            return version != null && version.getMajorVersion() == ALLOWED_MAJOR_VERSION;
        } catch ( FileNotFoundException e ) {
            return false;
        } catch ( IOException e ) {
            return false;
        } finally {
            try {
                if ( fis != null ) fis.close();
            } catch ( IOException e ) {
                // do nothing
            }
        }
    }

    // --------------------------------------------------------------------------------
    //
    // implicit block
    //
    // The first four records of BCF are inline untype encoded data of:
    //
    // 4 byte integer chrom offset
    // 4 byte integer start
    // 4 byte integer ref length
    // 4 byte float qual
    //
    // --------------------------------------------------------------------------------

    /**
     * Decode the sites level data from this classes decoder
     *
     * @param builder
     * @return
     */
    @Requires({"builder != null"})
    private final void decodeSiteLoc(final VariantContextBuilder builder) throws IOException {
        final int contigOffset = decoder.decodeInt(BCF2Type.INT32);
        final String contig = lookupContigName(contigOffset);
        builder.chr(contig);

        this.pos = decoder.decodeInt(BCF2Type.INT32) + 1; // GATK is one based, BCF2 is zero-based
        final int refLength = decoder.decodeInt(BCF2Type.INT32);
        builder.start((long)pos);
        builder.stop((long)(pos + refLength - 1)); // minus one because GATK has closed intervals but BCF2 is open
    }

    /**
     * Decode the sites level data from this classes decoder
     *
     * @param builder
     * @return
     */
    @Requires({"builder != null", "decoder != null"})
    @Ensures({"result != null", "result.isValid()"})
    private final SitesInfoForDecoding decodeSitesExtendedInfo(final VariantContextBuilder builder) throws IOException {
        final Object qual = decoder.decodeSingleValue(BCF2Type.FLOAT);
        if ( qual != null ) {
            builder.log10PError(((Double)qual) / -10.0);
        }

        final int nAlleleInfo = decoder.decodeInt(BCF2Type.INT32);
        final int nFormatSamples = decoder.decodeInt(BCF2Type.INT32);
        final int nAlleles = nAlleleInfo >> 16;
        final int nInfo = nAlleleInfo & 0x0000FFFF;
        final int nFormatFields = nFormatSamples >> 24;
        final int nSamples = nFormatSamples & 0x00FFFFF;

        if ( header.getNGenotypeSamples() != nSamples )
            error("Reading BCF2 files with different numbers of samples per record " +
                    "is not currently supported.  Saw " + header.getNGenotypeSamples() +
                    " samples in header but have a record with " + nSamples + " samples");

        decodeID(builder);
        final List<Allele> alleles = decodeAlleles(builder, pos, nAlleles);
        decodeFilter(builder);
        decodeInfo(builder, nInfo);

        final SitesInfoForDecoding info = new SitesInfoForDecoding(nFormatFields, nSamples, alleles);
        if ( ! info.isValid() )
            error("Sites info is malformed: " + info);
        return info;
    }

    protected final static class SitesInfoForDecoding {
        final int nFormatFields;
        final int nSamples;
        final List<Allele> alleles;

        private SitesInfoForDecoding(final int nFormatFields, final int nSamples, final List<Allele> alleles) {
            this.nFormatFields = nFormatFields;
            this.nSamples = nSamples;
            this.alleles = alleles;
        }

        public boolean isValid() {
            return nFormatFields >= 0 &&
                    nSamples >= 0 &&
                    alleles != null && ! alleles.isEmpty() && alleles.get(0).isReference();
        }

        @Override
        public String toString() {
            return String.format("nFormatFields = %d, nSamples = %d, alleles = %s", nFormatFields, nSamples, alleles);
        }
    }

    /**
     * Decode the id field in this BCF2 file and store it in the builder
     * @param builder
     */
    private void decodeID( final VariantContextBuilder builder ) throws IOException {
        final String id = (String)decoder.decodeTypedValue();

        if ( id == null )
            builder.noID();
        else
            builder.id(id);
    }

    /**
     * Decode the alleles from this BCF2 file and put the results in builder
     * @param builder
     * @param pos
     * @param nAlleles
     * @return the alleles
     */
    @Requires("nAlleles > 0")
    private List<Allele> decodeAlleles( final VariantContextBuilder builder, final int pos, final int nAlleles ) throws IOException {
        // TODO -- probably need inline decoder for efficiency here (no sense in going bytes -> string -> vector -> bytes
        List<Allele> alleles = new ArrayList<Allele>(nAlleles);
        String ref = null;

        for ( int i = 0; i < nAlleles; i++ ) {
            final String alleleBases = (String)decoder.decodeTypedValue();

            final boolean isRef = i == 0;
            final Allele allele = Allele.create(alleleBases, isRef);
            if ( isRef ) ref = alleleBases;

            alleles.add(allele);
        }
        assert ref != null;

        builder.alleles(alleles);

        assert ref.length() > 0;

        return alleles;
    }

    /**
     * Decode the filter field of this BCF2 file and store the result in the builder
     * @param builder
     */
    private void decodeFilter( final VariantContextBuilder builder ) throws IOException {
        final Object value = decoder.decodeTypedValue();

        if ( value == null )
            builder.unfiltered();
        else {
            if ( value instanceof Integer ) {
                // fast path for single integer result
                final String filterString = getDictionaryString((Integer)value);
                if ( VCFConstants.PASSES_FILTERS_v4.equals(filterString))
                    builder.passFilters();
                else
                    builder.filter(filterString);
            } else {
                for ( final int offset : (List<Integer>)value )
                    builder.filter(getDictionaryString(offset));
            }
        }
    }

    /**
     * Loop over the info field key / value pairs in this BCF2 file and decode them into the builder
     *
     * @param builder
     * @param numInfoFields
     */
    @Requires("numInfoFields >= 0")
    private void decodeInfo( final VariantContextBuilder builder, final int numInfoFields ) throws IOException {
        if ( numInfoFields == 0 )
            // fast path, don't bother doing any work if there are no fields
            return;

        final Map<String, Object> infoFieldEntries = new HashMap<String, Object>(numInfoFields);
        for ( int i = 0; i < numInfoFields; i++ ) {
            final String key = getDictionaryString();
            Object value = decoder.decodeTypedValue();
            final VCFCompoundHeaderLine metaData = VariantContextUtils.getMetaDataForField(header, key);
            if ( metaData.getType() == VCFHeaderLineType.Flag ) value = true; // special case for flags
            infoFieldEntries.put(key, value);
        }

        builder.attributes(infoFieldEntries);
    }

    // --------------------------------------------------------------------------------
    //
    // Decoding Genotypes
    //
    // --------------------------------------------------------------------------------

    /**
     * Create the lazy loader for the genotypes data, and store it in the builder
     * so that the VC will be able to decode on demand the genotypes data
     *
     * @param siteInfo
     * @param builder
     */
    private void createLazyGenotypesDecoder( final SitesInfoForDecoding siteInfo,
                                             final VariantContextBuilder builder ) {
        if (siteInfo.nSamples > 0) {
            final LazyGenotypesContext.LazyParser lazyParser =
                    new BCF2LazyGenotypesDecoder(this, siteInfo.alleles, siteInfo.nSamples, siteInfo.nFormatFields, builders);

            final LazyData lazyData = new LazyData(header, siteInfo.nFormatFields, decoder.getRecordBytes());
            final LazyGenotypesContext lazy = new LazyGenotypesContext(lazyParser, lazyData, header.getNGenotypeSamples());

            // did we resort the sample names?  If so, we need to load the genotype data
            if ( !header.samplesWereAlreadySorted() )
                lazy.decode();

            builder.genotypesNoValidation(lazy);
        }
    }

    public static class LazyData {
        final public VCFHeader header;
        final public int nGenotypeFields;
        final public byte[] bytes;

        @Requires({"nGenotypeFields > 0", "bytes != null"})
        public LazyData(final VCFHeader header, final int nGenotypeFields, final byte[] bytes) {
            this.header = header;
            this.nGenotypeFields = nGenotypeFields;
            this.bytes = bytes;
        }
    }

    @Ensures("result != null")
    private final String getDictionaryString() throws IOException {
        return getDictionaryString((Integer) decoder.decodeTypedValue());
    }

    @Requires("offset < dictionary.size()")
    @Ensures("result != null")
    protected final String getDictionaryString(final int offset) {
        return dictionary.get(offset);
    }

    /**
     * Translate the config offset as encoded in the BCF file into the actual string
     * name of the contig from the dictionary
     *
     * @param contigOffset
     * @return
     */
    @Requires({"contigOffset >= 0", "contigOffset < contigNames.size()"})
    @Ensures("result != null")
    private final String lookupContigName( final int contigOffset ) {
        return contigNames.get(contigOffset);
    }

    @Requires("header != null")
    @Ensures({"result != null", "! result.isEmpty()"})
    private final ArrayList<String> parseDictionary(final VCFHeader header) {
        final ArrayList<String> dict = BCF2Utils.makeDictionary(header);

        // if we got here we never found a dictionary, or there are no elements in the dictionary
        if ( dict.isEmpty() )
            error("Dictionary header element was absent or empty");

        return dict;
    }

    /**
     * @return the VCFHeader we found in this BCF2 file
     */
    protected VCFHeader getHeader() {
        return header;
    }

    @Requires("field != null")
    @Ensures("result != null")
    protected BCF2GenotypeFieldDecoders.Decoder getGenotypeFieldDecoder(final String field) {
        return gtFieldDecoders.getDecoder(field);
    }

    private void error(final String message) throws RuntimeException {
        throw new TribbleException(String.format("%s, at record %d with position %d:", message, recordNo, pos));
    }
}
