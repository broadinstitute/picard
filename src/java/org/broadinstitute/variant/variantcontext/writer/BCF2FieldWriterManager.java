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
import org.broadinstitute.variant.utils.GeneralUtils;
import org.broadinstitute.variant.vcf.*;

import java.util.HashMap;
import java.util.Map;

/**
 * See #BCFWriter for documentation on this classes role in encoding BCF2 files
 *
 * @author Mark DePristo
 * @since 06/12
 */
public class BCF2FieldWriterManager {
    final Map<String, BCF2FieldWriter.SiteWriter> siteWriters = new HashMap<String, BCF2FieldWriter.SiteWriter>();
    final Map<String, BCF2FieldWriter.GenotypesWriter> genotypesWriters = new HashMap<String, BCF2FieldWriter.GenotypesWriter>();
    final IntGenotypeFieldAccessors intGenotypeFieldAccessors = new IntGenotypeFieldAccessors();

    public BCF2FieldWriterManager() { }

    /**
     * Setup the FieldWriters appropriate to each INFO and FORMAT in the VCF header
     *
     * Must be called before any of the getter methods will work
     *
     * @param header a VCFHeader containing description for every INFO and FORMAT field we'll attempt to write out to BCF
     * @param encoder the encoder we are going to use to write out the BCF2 data
     * @param stringDictionary a map from VCFHeader strings to their offsets for encoding
     */
    public void setup(final VCFHeader header, final BCF2Encoder encoder, final Map<String, Integer> stringDictionary) {
        for (final VCFInfoHeaderLine line : header.getInfoHeaderLines()) {
            final String field = line.getID();
            final BCF2FieldWriter.SiteWriter writer = createInfoWriter(header, line, encoder, stringDictionary);
            add(siteWriters, field, writer);
        }

        for (final VCFFormatHeaderLine line : header.getFormatHeaderLines()) {
            final String field = line.getID();
            final BCF2FieldWriter.GenotypesWriter writer = createGenotypesWriter(header, line, encoder, stringDictionary);
            add(genotypesWriters, field, writer);
        }
    }

    @Requires({"field != null", "writer != null"})
    @Ensures("map.containsKey(field)")
    private final <T> void add(final Map<String, T> map, final String field, final T writer) {
        if ( map.containsKey(field) )
            throw new IllegalStateException("BUG: field " + field + " already seen in VCFHeader while building BCF2 field encoders");
        map.put(field, writer);
    }

    // -----------------------------------------------------------------
    //
    // Master routine to look at the header, a specific line, and
    // build an appropriate SiteWriter for that header element
    //
    // -----------------------------------------------------------------

    private BCF2FieldWriter.SiteWriter createInfoWriter(final VCFHeader header,
                                                        final VCFInfoHeaderLine line,
                                                        final BCF2Encoder encoder,
                                                        final Map<String, Integer> dict) {
        return new BCF2FieldWriter.GenericSiteWriter(header, createFieldEncoder(line, encoder, dict, false));
    }

    private BCF2FieldEncoder createFieldEncoder(final VCFCompoundHeaderLine line,
                                                final BCF2Encoder encoder,
                                                final Map<String, Integer> dict,
                                                final boolean createGenotypesEncoders ) {

        if ( createGenotypesEncoders && intGenotypeFieldAccessors.getAccessor(line.getID()) != null ) {
            if ( GeneralUtils.DEBUG_MODE_ENABLED && line.getType() != VCFHeaderLineType.Integer )
                System.err.println("Warning: field " + line.getID() + " expected to encode an integer but saw " + line.getType() + " for record " + line);
            return new BCF2FieldEncoder.IntArray(line, dict);
        } else if ( createGenotypesEncoders && line.getID().equals(VCFConstants.GENOTYPE_KEY) ) {
            return new BCF2FieldEncoder.GenericInts(line, dict);
        } else {
            switch ( line.getType() ) {
                case Character:
                case String:
                    return new BCF2FieldEncoder.StringOrCharacter(line, dict);
                case Flag:
                    return new BCF2FieldEncoder.Flag(line, dict);
                case Float:
                    return new BCF2FieldEncoder.Float(line, dict);
                case Integer:
                    if ( line.isFixedCount() && line.getCount() == 1 )
                        return new BCF2FieldEncoder.AtomicInt(line, dict);
                    else
                        return new BCF2FieldEncoder.GenericInts(line, dict);
                default:
                    throw new IllegalArgumentException("Unexpected type for field " + line.getID());
            }
        }
    }

    // -----------------------------------------------------------------
    //
    // Master routine to look at the header, a specific line, and
    // build an appropriate Genotypes for that header element
    //
    // -----------------------------------------------------------------

    private BCF2FieldWriter.GenotypesWriter createGenotypesWriter(final VCFHeader header,
                                                                  final VCFFormatHeaderLine line,
                                                                  final BCF2Encoder encoder,
                                                                  final Map<String, Integer> dict) {
        final String field = line.getID();
        final BCF2FieldEncoder fieldEncoder = createFieldEncoder(line, encoder, dict, true);

        if ( field.equals(VCFConstants.GENOTYPE_KEY) ) {
            return new BCF2FieldWriter.GTWriter(header, fieldEncoder);
        } else if ( line.getID().equals(VCFConstants.GENOTYPE_FILTER_KEY) ) {
            return new BCF2FieldWriter.FTGenotypesWriter(header, fieldEncoder);
        } else if ( intGenotypeFieldAccessors.getAccessor(field) != null ) {
            return new BCF2FieldWriter.IGFGenotypesWriter(header, fieldEncoder, intGenotypeFieldAccessors.getAccessor(field));
        } else if ( line.getType() == VCFHeaderLineType.Integer ) {
            return new BCF2FieldWriter.IntegerTypeGenotypesWriter(header, fieldEncoder);
        } else {
            return new BCF2FieldWriter.StaticallyTypeGenotypesWriter(header, fieldEncoder);
        }
    }

    // -----------------------------------------------------------------
    //
    // Accessors to get site / genotype writers
    //
    // -----------------------------------------------------------------

    /**
     * Get a site writer specialized to encode values for site info field
     * @param field key found in the VCF header INFO records
     * @return non-null writer if one can be found, or null if none exists for field
     */
    public BCF2FieldWriter.SiteWriter getSiteFieldWriter(final String field) {
        return getWriter(field, siteWriters);
    }

    /**
     * Get a genotypes writer specialized to encode values for genotypes field
     * @param field key found in the VCF header FORMAT records
     * @return non-null writer if one can be found, or null if none exists for field
     */
    public BCF2FieldWriter.GenotypesWriter getGenotypeFieldWriter(final String field) {
        return getWriter(field, genotypesWriters);
    }

    @Requires({"map != null", "key != null"})
    public <T> T getWriter(final String key, final Map<String, T> map) {
        return map.get(key);
    }
}
