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
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.GenotypeBuilder;

import java.io.IOException;
import java.util.*;

/**
 * An efficient scheme for building and obtaining specialized
 * genotype field decoders.  Used by the BCFCodec to parse
 * with little overhead the fields from BCF2 encoded genotype
 * records
 *
 * @author Mark DePristo
 * @since 6/12
 */
public class BCF2GenotypeFieldDecoders {
    private final static boolean ENABLE_FASTPATH_GT = true;
    private final static int MIN_SAMPLES_FOR_FASTPATH_GENOTYPES = 0; // TODO -- update to reasonable number

    // initialized once per writer to allow parallel writers to work
    private final HashMap<String, Decoder> genotypeFieldDecoder = new HashMap<String, Decoder>();
    private final Decoder defaultDecoder = new GenericDecoder();

    public BCF2GenotypeFieldDecoders(final VCFHeader header) {
        // TODO -- fill in appropriate decoders for each FORMAT field in the header

        genotypeFieldDecoder.put(VCFConstants.GENOTYPE_KEY, new GTDecoder());
        // currently the generic decoder handles FILTER values properly, in so far as we don't tolerate multiple filter field values per genotype
        genotypeFieldDecoder.put(VCFConstants.GENOTYPE_FILTER_KEY, new FTDecoder());
        genotypeFieldDecoder.put(VCFConstants.DEPTH_KEY, new DPDecoder());
        genotypeFieldDecoder.put(VCFConstants.GENOTYPE_ALLELE_DEPTHS, new ADDecoder());
        genotypeFieldDecoder.put(VCFConstants.GENOTYPE_PL_KEY, new PLDecoder());
        genotypeFieldDecoder.put(VCFConstants.GENOTYPE_QUALITY_KEY, new GQDecoder());
    }

    // -----------------------------------------------------------------
    //
    // Genotype field decoder
    //
    // -----------------------------------------------------------------

    /**
     * Return decoder appropriate for field, or the generic decoder if no
     * specialized one is bound
     * @param field the GT field to decode
     * @return a non-null decoder
     */
    @Requires("field != null")
    @Ensures("result != null")
    public Decoder getDecoder(final String field) {
        final Decoder d = genotypeFieldDecoder.get(field);
        return d == null ? defaultDecoder : d;
    }

    /**
     * Decoder a field (implicit from creation) encoded as
     * typeDescriptor in the decoder object in the GenotypeBuilders
     * one for each sample in order.
     *
     * The way this works is that this decode method
     * iterates over the builders, decoding a genotype field
     * in BCF2 for each sample from decoder.
     *
     * This system allows us to easily use specialized
     * decoders for specific genotype field values. For example,
     * we use a special decoder to directly read the BCF2 data for
     * the PL field into a int[] rather than the generic List of Integer
     */
    public interface Decoder {
        @Requires({"siteAlleles != null", "! siteAlleles.isEmpty()",
                "field != null", "decoder != null", "gbs != null", "gbs.length != 0"})
        public void decode(final List<Allele> siteAlleles,
                           final String field,
                           final BCF2Decoder decoder,
                           final byte typeDescriptor,
                           final int numElements,
                           final GenotypeBuilder[] gbs) throws IOException;
    }

    private class GTDecoder implements Decoder {
        @Override
        public void decode(final List<Allele> siteAlleles, final String field, final BCF2Decoder decoder, final byte typeDescriptor, final int numElements, final GenotypeBuilder[] gbs) throws IOException {
            if ( ENABLE_FASTPATH_GT && siteAlleles.size() == 2 && numElements == 2 && gbs.length >= MIN_SAMPLES_FOR_FASTPATH_GENOTYPES )
                fastBiallelicDiploidDecode(siteAlleles, decoder, typeDescriptor, gbs);
            else {
                generalDecode(siteAlleles, numElements, decoder, typeDescriptor, gbs);
            }
        }

        /**
         * fast path for many samples with diploid genotypes
         *
         * The way this would work is simple.  Create a List<Allele> diploidGenotypes[] object
         * After decoding the offset, if that sample is diploid compute the
         * offset into the alleles vector which is simply offset = allele0 * nAlleles + allele1
         * if there's a value at diploidGenotypes[offset], use it, otherwise create the genotype
         * cache it and use that
         *
         * Some notes.  If there are nAlleles at the site, there are implicitly actually
         * n + 1 options including
         */
        @Requires("siteAlleles.size() == 2")
        @SuppressWarnings({"unchecked"})
        private final void fastBiallelicDiploidDecode(final List<Allele> siteAlleles,
                                                      final BCF2Decoder decoder,
                                                      final byte typeDescriptor,
                                                      final GenotypeBuilder[] gbs) throws IOException {
            final BCF2Type type = BCF2Utils.decodeType(typeDescriptor);

            final int nPossibleGenotypes = 3 * 3;
            final Object allGenotypes[] = new Object[nPossibleGenotypes];

            for ( final GenotypeBuilder gb : gbs ) {
                final int a1 = decoder.decodeInt(type);
                final int a2 = decoder.decodeInt(type);

                if ( a1 == type.getMissingBytes() ) {
                    assert a2 == type.getMissingBytes();
                    // no called sample GT = .
                    gb.alleles(null);
                } else if ( a2 == type.getMissingBytes() ) {
                    gb.alleles(Arrays.asList(getAlleleFromEncoded(siteAlleles, a1)));
                } else {
                    // downshift to remove phase
                    final int offset = (a1 >> 1) * 3 + (a2 >> 1);
                    assert offset < allGenotypes.length;

                    // TODO -- how can I get rid of this cast?
                    List<Allele> gt = (List<Allele>)allGenotypes[offset];
                    if ( gt == null ) {
                        final Allele allele1 = getAlleleFromEncoded(siteAlleles, a1);
                        final Allele allele2 = getAlleleFromEncoded(siteAlleles, a2);
                        gt = Arrays.asList(allele1, allele2);
                        allGenotypes[offset] = gt;
                    }

                    gb.alleles(gt);
                }

                final boolean phased = (a1 & 0x01) == 1;
                gb.phased(phased);
            }
        }

        private final void generalDecode(final List<Allele> siteAlleles,
                                         final int ploidy,
                                         final BCF2Decoder decoder,
                                         final byte typeDescriptor,
                                         final GenotypeBuilder[] gbs) throws IOException {
            final BCF2Type type = BCF2Utils.decodeType(typeDescriptor);

            // a single cache for the encoded genotypes, since we don't actually need this vector
            final int[] tmp = new int[ploidy];

            for ( final GenotypeBuilder gb : gbs ) {
                final int[] encoded = decoder.decodeIntArray(ploidy, type, tmp);
                if ( encoded == null )
                    // no called sample GT = .
                    gb.alleles(null);
                else {
                    assert encoded.length > 0;

                    // we have at least some alleles to decode
                    final List<Allele> gt = new ArrayList<Allele>(encoded.length);

                    // note that the auto-pruning of fields magically handles different
                    // ploidy per sample at a site
                    for ( final int encode : encoded )
                        gt.add(getAlleleFromEncoded(siteAlleles, encode));

                    gb.alleles(gt);
                    final boolean phased = (encoded[0] & 0x01) == 1;
                    gb.phased(phased);
                }
            }
        }

        @Requires({"siteAlleles != null && ! siteAlleles.isEmpty()", "encode >= 0"})
        @Ensures("result != null")
        private final Allele getAlleleFromEncoded(final List<Allele> siteAlleles, final int encode) {
            final int offset = encode >> 1;
            return offset == 0 ? Allele.NO_CALL : siteAlleles.get(offset - 1);
        }
    }

    private class DPDecoder implements Decoder {
        @Override
        public void decode(final List<Allele> siteAlleles, final String field, final BCF2Decoder decoder, final byte typeDescriptor, final int numElements, final GenotypeBuilder[] gbs) throws IOException {
            for ( final GenotypeBuilder gb : gbs ) {
                // the -1 is for missing
                gb.DP(decoder.decodeInt(typeDescriptor, -1));
            }
        }
    }

    private class GQDecoder implements Decoder {
        @Override
        public void decode(final List<Allele> siteAlleles, final String field, final BCF2Decoder decoder, final byte typeDescriptor, final int numElements, final GenotypeBuilder[] gbs) throws IOException {
            for ( final GenotypeBuilder gb : gbs ) {
                // the -1 is for missing
                gb.GQ(decoder.decodeInt(typeDescriptor, -1));
            }
        }
    }

    private class ADDecoder implements Decoder {
        @Override
        public void decode(final List<Allele> siteAlleles, final String field, final BCF2Decoder decoder, final byte typeDescriptor, final int numElements, final GenotypeBuilder[] gbs) throws IOException {
            for ( final GenotypeBuilder gb : gbs ) {
                gb.AD(decoder.decodeIntArray(typeDescriptor, numElements));
            }
        }
    }

    private class PLDecoder implements Decoder {
        @Override
        public void decode(final List<Allele> siteAlleles, final String field, final BCF2Decoder decoder, final byte typeDescriptor, final int numElements, final GenotypeBuilder[] gbs) throws IOException {
            for ( final GenotypeBuilder gb : gbs ) {
                gb.PL(decoder.decodeIntArray(typeDescriptor, numElements));
            }
        }
    }

    private class GenericDecoder implements Decoder {
        @Override
        public void decode(final List<Allele> siteAlleles, final String field, final BCF2Decoder decoder, final byte typeDescriptor, final int numElements, final GenotypeBuilder[] gbs) throws IOException {
            for ( final GenotypeBuilder gb : gbs ) {
                Object value = decoder.decodeTypedValue(typeDescriptor, numElements);
                if ( value != null ) { // don't add missing values
                    if ( value instanceof List && ((List)value).size() == 1) {
                        // todo -- I really hate this, and it suggests that the code isn't completely right
                        // the reason it's here is that it's possible to prune down a vector to a singleton
                        // value and there we have the contract that the value comes back as an atomic value
                        // not a vector of size 1
                        value = ((List)value).get(0);
                    }
                    gb.attribute(field, value);
                }
            }
        }
    }

    private class FTDecoder implements Decoder {
        @Override
        public void decode(final List<Allele> siteAlleles, final String field, final BCF2Decoder decoder, final byte typeDescriptor, final int numElements, final GenotypeBuilder[] gbs) throws IOException {
            for ( final GenotypeBuilder gb : gbs ) {
                Object value = decoder.decodeTypedValue(typeDescriptor, numElements);
                assert value == null || value instanceof String;
                gb.filter((String)value);
            }
        }
    }
}
