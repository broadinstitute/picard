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
import org.broadinstitute.variant.bcf2.BCF2Codec;
import org.broadinstitute.variant.bcf2.BCF2Type;
import org.broadinstitute.variant.bcf2.BCF2Utils;
import org.broadinstitute.variant.bcf2.BCFVersion;
import org.broadinstitute.variant.utils.GeneralUtils;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.variant.vcf.VCFContigHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.vcf.VCFUtils;

import java.io.*;
import java.util.*;

/**
 * VariantContextWriter that emits BCF2 binary encoding
 *
 * Overall structure of this writer is complex for efficiency reasons
 *
 * -- The BCF2Writer manages the low-level BCF2 encoder, the mappings
 * from contigs and strings to offsets, the VCF header, and holds the
 * lower-level encoders that map from VC and Genotype fields to their
 * specific encoders.  This class also writes out the standard BCF2 fields
 * like POS, contig, the size of info and genotype data, QUAL, etc.  It
 * has loops over the INFO and GENOTYPES to encode each individual datum
 * with the generic field encoders, but the actual encoding work is
 * done with by the FieldWriters classes themselves
 *
 * -- BCF2FieldWriter are specialized classes for writing out SITE and
 * genotype information for specific SITE/GENOTYPE fields (like AC for
 * sites and GQ for genotypes).  These are objects in themselves because
 * the manage all of the complexity of relating the types in the VCF header
 * with the proper encoding in BCF as well as the type representing this
 * in java.  Relating all three of these pieces of information together
 * is the main complexity challenge in the encoder.  The piece of code
 * that determines which FieldWriters to associate with each SITE and
 * GENOTYPE field is the BCF2FieldWriterManager.  These FieldWriters
 * are specialized for specific combinations of encoders (see below)
 * and contexts (genotypes) for efficiency, so they smartly manage
 * the writing of PLs (encoded as int[]) directly into the lowest
 * level BCFEncoder.
 *
 * -- At the third level is the BCF2FieldEncoder, relatively simple
 * pieces of code that handle the task of determining the right
 * BCF2 type for specific field values, as well as reporting back
 * information such as the number of elements used to encode it
 * (simple for atomic values like Integer but complex for PLs
 * or lists of strings)
 *
 * -- At the lowest level is the BCF2Encoder itself.  This provides
 * just the limited encoding methods specified by the BCF2 specification.  This encoder
 * doesn't do anything but make it possible to conveniently write out valid low-level
 * BCF2 constructs.
 *
 * @author Mark DePristo
 * @since 06/12
 */
class BCF2Writer extends IndexingVariantContextWriter {
    public static final int MAJOR_VERSION = 2;
    public static final int MINOR_VERSION = 1;

    final private static boolean ALLOW_MISSING_CONTIG_LINES = false;

    private final OutputStream outputStream;      // Note: do not flush until completely done writing, to avoid issues with eventual BGZF support
    private VCFHeader header;
    private final Map<String, Integer> contigDictionary = new HashMap<String, Integer>();
    private final Map<String, Integer> stringDictionaryMap = new LinkedHashMap<String, Integer>();
    private final boolean doNotWriteGenotypes;
    private String[] sampleNames = null;

    private final BCF2Encoder encoder = new BCF2Encoder(); // initialized after the header arrives
    final BCF2FieldWriterManager fieldManager = new BCF2FieldWriterManager();

    /**
     * cached results for whether we can write out raw genotypes data.
     */
    private VCFHeader lastVCFHeaderOfUnparsedGenotypes = null;
    private boolean canPassOnUnparsedGenotypeDataForLastVCFHeader = false;


    public BCF2Writer(final File location, final OutputStream output, final SAMSequenceDictionary refDict, final boolean enableOnTheFlyIndexing, final boolean doNotWriteGenotypes) {
        super(writerName(location, output), location, output, refDict, enableOnTheFlyIndexing);
        this.outputStream = getOutputStream();
        this.doNotWriteGenotypes = doNotWriteGenotypes;
    }

    // --------------------------------------------------------------------------------
    //
    // Interface functions
    //
    // --------------------------------------------------------------------------------

    @Override
    public void writeHeader(VCFHeader header) {
        // make sure the header is sorted correctly
        header = new VCFHeader(header.getMetaDataInSortedOrder(), header.getGenotypeSamples());

        // create the config offsets map
        if ( header.getContigLines().isEmpty() ) {
            if ( ALLOW_MISSING_CONTIG_LINES ) {
                if ( GeneralUtils.DEBUG_MODE_ENABLED ) {
                    System.err.println("No contig dictionary found in header, falling back to reference sequence dictionary");
                }
                createContigDictionary(VCFUtils.makeContigHeaderLines(getRefDict(), null));
            } else {
                throw new IllegalStateException("Cannot write BCF2 file with missing contig lines");
            }
        } else {
            createContigDictionary(header.getContigLines());
        }

        // set up the map from dictionary string values -> offset
        final ArrayList<String> dict = BCF2Utils.makeDictionary(header);
        for ( int i = 0; i < dict.size(); i++ ) {
            stringDictionaryMap.put(dict.get(i), i);
        }

        sampleNames = header.getGenotypeSamples().toArray(new String[header.getNGenotypeSamples()]);

        // setup the field encodings
        fieldManager.setup(header, encoder, stringDictionaryMap);

        try {
            // write out the header into a byte stream, get its length, and write everything to the file
            final ByteArrayOutputStream capture = new ByteArrayOutputStream();
            final OutputStreamWriter writer = new OutputStreamWriter(capture);
            this.header = VCFWriter.writeHeader(header, writer, doNotWriteGenotypes, VCFWriter.getVersionLine(), "BCF2 stream");
            writer.append('\0'); // the header is null terminated by a byte
            writer.close();

            final byte[] headerBytes = capture.toByteArray();
            new BCFVersion(MAJOR_VERSION, MINOR_VERSION).write(outputStream);
            BCF2Type.INT32.write(headerBytes.length, outputStream);
            outputStream.write(headerBytes);
        } catch (IOException e) {
            throw new RuntimeException("BCF2 stream: Got IOException while trying to write BCF2 header", e);
        }
    }

    @Override
    public void add( VariantContext vc ) {
        if ( doNotWriteGenotypes )
            vc = new VariantContextBuilder(vc).noGenotypes().make();
        vc = vc.fullyDecode(header, false);

        super.add(vc); // allow on the fly indexing

        try {
            final byte[] infoBlock = buildSitesData(vc);
            final byte[] genotypesBlock = buildSamplesData(vc);

            // write the two blocks to disk
            writeBlock(infoBlock, genotypesBlock);
        }
        catch ( IOException e ) {
            throw new RuntimeException("Error writing record to BCF2 file: " + vc.toString(), e);
        }
    }

    @Override
    public void close() {
        try {
            outputStream.flush();
            outputStream.close();
        }
        catch ( IOException e ) {
            throw new RuntimeException("Failed to close BCF2 file");
        }
        super.close();
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
    private byte[] buildSitesData( VariantContext vc ) throws IOException {
        final int contigIndex = contigDictionary.get(vc.getChr());
        if ( contigIndex == -1 )
            throw new IllegalStateException(String.format("Contig %s not found in sequence dictionary from reference", vc.getChr()));

        // note use of encodeRawValue to not insert the typing byte
        encoder.encodeRawValue(contigIndex, BCF2Type.INT32);

        // pos.  GATK is 1 based, BCF2 is 0 based
        encoder.encodeRawValue(vc.getStart() - 1, BCF2Type.INT32);

        // ref length.  GATK is closed, but BCF2 is open so the ref length is GATK end - GATK start + 1
        // for example, a SNP is in GATK at 1:10-10, which has ref length 10 - 10 + 1 = 1
        encoder.encodeRawValue(vc.getEnd() - vc.getStart() + 1, BCF2Type.INT32);

        // qual
        if ( vc.hasLog10PError() )
            encoder.encodeRawFloat((float) vc.getPhredScaledQual());
        else
            encoder.encodeRawMissingValue(BCF2Type.FLOAT);

        // info fields
        final int nAlleles = vc.getNAlleles();
        final int nInfo = vc.getAttributes().size();
        final int nGenotypeFormatFields = getNGenotypeFormatFields(vc);
        final int nSamples = header.getNGenotypeSamples();

        encoder.encodeRawInt((nAlleles << 16) | (nInfo & 0x0000FFFF), BCF2Type.INT32);
        encoder.encodeRawInt((nGenotypeFormatFields << 24) | (nSamples & 0x00FFFFF), BCF2Type.INT32);

        buildID(vc);
        buildAlleles(vc);
        buildFilter(vc);
        buildInfo(vc);

        return encoder.getRecordBytes();
    }


    /**
     * Can we safely write on the raw (undecoded) genotypes of an input VC?
     *
     * The cache depends on the undecoded lazy data header == lastVCFHeaderOfUnparsedGenotypes, in
     * which case we return the previous result.  If it's not cached, we use the BCF2Util to
     * compare the VC header with our header (expensive) and cache it.
     *
     * @param lazyData
     * @return
     */
    private boolean canSafelyWriteRawGenotypesBytes(final BCF2Codec.LazyData lazyData) {
        if ( lazyData.header != lastVCFHeaderOfUnparsedGenotypes ) {
            // result is already cached
            canPassOnUnparsedGenotypeDataForLastVCFHeader = BCF2Utils.headerLinesAreOrderedConsistently(this.header,lazyData.header);
            lastVCFHeaderOfUnparsedGenotypes = lazyData.header;
        }

        return canPassOnUnparsedGenotypeDataForLastVCFHeader;
    }

    private BCF2Codec.LazyData getLazyData(final VariantContext vc) {
        if ( vc.getGenotypes().isLazyWithData() ) {
            final LazyGenotypesContext lgc = (LazyGenotypesContext)vc.getGenotypes();

            if ( lgc.getUnparsedGenotypeData() instanceof BCF2Codec.LazyData &&
                    canSafelyWriteRawGenotypesBytes((BCF2Codec.LazyData) lgc.getUnparsedGenotypeData())) {
                return (BCF2Codec.LazyData)lgc.getUnparsedGenotypeData();
            } else {
                lgc.decode(); // WARNING -- required to avoid keeping around bad lazy data for too long
            }
        }

        return null;
    }

    /**
     * Try to get the nGenotypeFields as efficiently as possible.
     *
     * If this is a lazy BCF2 object just grab the field count from there,
     * otherwise do the whole counting by types test in the actual data
     *
     * @param vc
     * @return
     */
    private int getNGenotypeFormatFields(final VariantContext vc) {
        final BCF2Codec.LazyData lazyData = getLazyData(vc);
        return lazyData != null ? lazyData.nGenotypeFields : vc.calcVCFGenotypeKeys(header).size();
    }

    private void buildID( VariantContext vc ) throws IOException {
        encoder.encodeTypedString(vc.getID());
    }

    private void buildAlleles( VariantContext vc ) throws IOException {
        for ( Allele allele : vc.getAlleles() ) {
            final byte[] s = allele.getDisplayBases();
            if ( s == null )
                throw new IllegalStateException("BUG: BCF2Writer encountered null padded allele" + allele);
            encoder.encodeTypedString(s);
        }
    }

    private void buildFilter( VariantContext vc ) throws IOException {
        if ( vc.isFiltered() ) {
            encodeStringsByRef(vc.getFilters());
        } else if ( vc.filtersWereApplied() ) {
            encodeStringsByRef(Collections.singleton(VCFConstants.PASSES_FILTERS_v4));
        } else {
            encoder.encodeTypedMissing(BCF2Type.INT8);
        }
    }

    private void buildInfo( VariantContext vc ) throws IOException {
        for ( Map.Entry<String, Object> infoFieldEntry : vc.getAttributes().entrySet() ) {
            final String field = infoFieldEntry.getKey();
            final BCF2FieldWriter.SiteWriter writer = fieldManager.getSiteFieldWriter(field);
            if ( writer == null ) errorUnexpectedFieldToWrite(vc, field, "INFO");
            writer.start(encoder, vc);
            writer.site(encoder, vc);
            writer.done(encoder, vc);
        }
    }

    private byte[] buildSamplesData(final VariantContext vc) throws IOException {
        final BCF2Codec.LazyData lazyData = getLazyData(vc);  // has critical side effects
        if ( lazyData != null ) {
            // we never decoded any data from this BCF file, so just pass it back
            return lazyData.bytes;
        }

        // we have to do work to convert the VC into a BCF2 byte stream
        final List<String> genotypeFields = vc.calcVCFGenotypeKeys(header);
        for ( final String field : genotypeFields ) {
            final BCF2FieldWriter.GenotypesWriter writer = fieldManager.getGenotypeFieldWriter(field);
            if ( writer == null ) errorUnexpectedFieldToWrite(vc, field, "FORMAT");

            assert writer != null;

            writer.start(encoder, vc);
            for ( final String name : sampleNames ) {
                Genotype g = vc.getGenotype(name);
                if ( g == null ) g = GenotypeBuilder.createMissing(name, writer.nValuesPerGenotype);
                writer.addGenotype(encoder, vc, g);
            }
            writer.done(encoder, vc);
        }
        return encoder.getRecordBytes();
    }

    /**
     * Throws a meaningful error message when a field (INFO or FORMAT) is found when writing out a file
     * but there's no header line for it.
     *
     * @param vc
     * @param field
     * @param fieldType
     */
    private void errorUnexpectedFieldToWrite(final VariantContext vc, final String field, final String fieldType) {
        throw new IllegalStateException("Found field " + field + " in the " + fieldType + " fields of VariantContext at " +
                vc.getChr() + ":" + vc.getStart() + " from " + vc.getSource() + " but this hasn't been defined in the VCFHeader");
    }

    // --------------------------------------------------------------------------------
    //
    // Low-level block encoding
    //
    // --------------------------------------------------------------------------------

    /**
     * Write the data in the encoder to the outputstream as a length encoded
     * block of data.  After this call the encoder stream will be ready to
     * start a new data block
     *
     * @throws IOException
     */
    @Requires({"infoBlock.length > 0", "genotypesBlock.length >= 0"})
    private void writeBlock(final byte[] infoBlock, final byte[] genotypesBlock) throws IOException {
        BCF2Type.INT32.write(infoBlock.length, outputStream);
        BCF2Type.INT32.write(genotypesBlock.length, outputStream);
        outputStream.write(infoBlock);
        outputStream.write(genotypesBlock);
    }

    @Requires("! strings.isEmpty()")
    @Ensures("result.isIntegerType()")
    private BCF2Type encodeStringsByRef(final Collection<String> strings) throws IOException {
        final List<Integer> offsets = new ArrayList<Integer>(strings.size());

        // iterate over strings until we find one that needs 16 bits, and break
        for ( final String string : strings ) {
            final Integer got = stringDictionaryMap.get(string);
            if ( got == null ) throw new IllegalStateException("Format error: could not find string " + string + " in header as required by BCF");
            final int offset = got;
            offsets.add(offset);
        }

        final BCF2Type type = BCF2Utils.determineIntegerType(offsets);
        encoder.encodeTyped(offsets, type);
        return type;
    }

    /**
     * Create the contigDictionary from the contigLines extracted from the VCF header
     *
     * @param contigLines
     */
    @Requires("contigDictionary.isEmpty()")
    private void createContigDictionary(final Collection<VCFContigHeaderLine> contigLines) {
        int offset = 0;
        for ( VCFContigHeaderLine contig : contigLines )
            contigDictionary.put(contig.getID(), offset++);
    }
}
