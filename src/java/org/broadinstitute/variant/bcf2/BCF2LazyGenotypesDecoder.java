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

import com.google.java.contract.Requires;
import org.broad.tribble.TribbleException;
import org.broadinstitute.variant.variantcontext.*;

import java.io.IOException;
import java.util.*;

/**
 * Lazy version of genotypes decoder for BCF2 genotypes
 *
 * @author Mark DePristo
 * @since 5/12
 */
public class BCF2LazyGenotypesDecoder implements LazyGenotypesContext.LazyParser {
    // the essential information for us to use to decode the genotypes data
    // initialized when this lazy decoder is created, as we know all of this from the BCF2Codec
    // and its stored here again for code cleanliness
    private final BCF2Codec codec;
    private final List<Allele> siteAlleles;
    private final int nSamples;
    private final int nFields;
    private final GenotypeBuilder[] builders;

    @Requires("codec.getHeader().getNGenotypeSamples() == builders.length")
    BCF2LazyGenotypesDecoder(final BCF2Codec codec, final List<Allele> alleles, final int nSamples,
                             final int nFields, final GenotypeBuilder[] builders) {
        this.codec = codec;
        this.siteAlleles = alleles;
        this.nSamples = nSamples;
        this.nFields = nFields;
        this.builders = builders;
    }

    @Override
    public LazyGenotypesContext.LazyData parse(final Object data) {
        try {

            // load our byte[] data into the decoder
            final BCF2Decoder decoder = new BCF2Decoder(((BCF2Codec.LazyData)data).bytes);

            for ( int i = 0; i < nSamples; i++ )
                builders[i].reset(true);

            for ( int i = 0; i < nFields; i++ ) {
                // get the field name
                final int offset = (Integer) decoder.decodeTypedValue();
                final String field = codec.getDictionaryString(offset);

                // the type of each element
                final byte typeDescriptor = decoder.readTypeDescriptor();
                final int numElements = decoder.decodeNumberOfElements(typeDescriptor);
                final BCF2GenotypeFieldDecoders.Decoder fieldDecoder = codec.getGenotypeFieldDecoder(field);
                try {
                    fieldDecoder.decode(siteAlleles, field, decoder, typeDescriptor, numElements, builders);
                } catch ( ClassCastException e ) {
                    throw new TribbleException("BUG: expected encoding of field " + field
                            + " inconsistent with the value observed in the decoded value");
                }
            }

            final ArrayList<Genotype> genotypes = new ArrayList<Genotype>(nSamples);
            for ( final GenotypeBuilder gb : builders )
                genotypes.add(gb.make());

            return new LazyGenotypesContext.LazyData(genotypes, codec.getHeader().getSampleNamesInOrder(), codec.getHeader().getSampleNameToOffset());
        } catch ( IOException e ) {
            throw new TribbleException("Unexpected IOException parsing already read genotypes data block", e);
        }
    }
}
