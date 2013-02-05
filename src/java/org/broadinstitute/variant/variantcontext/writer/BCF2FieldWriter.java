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
import org.broadinstitute.variant.bcf2.BCF2Type;
import org.broadinstitute.variant.bcf2.BCF2Utils;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * See #BCFWriter for documentation on this classes role in encoding BCF2 files
 *
 * @author Mark DePristo
 * @since 06/12
 */
public abstract class BCF2FieldWriter {
    private final VCFHeader header;
    private final BCF2FieldEncoder fieldEncoder;

    @Requires({"header != null", "fieldEncoder != null"})
    protected BCF2FieldWriter(final VCFHeader header, final BCF2FieldEncoder fieldEncoder) {
        this.header = header;
        this.fieldEncoder = fieldEncoder;
    }

    @Ensures("result != null")
    protected VCFHeader getHeader() { return header; }
    @Ensures("result != null")
    protected BCF2FieldEncoder getFieldEncoder() {
        return fieldEncoder;
    }
    @Ensures("result != null")
    protected String getField() { return getFieldEncoder().getField(); }

    @Requires("vc != null")
    public void start(final BCF2Encoder encoder, final VariantContext vc) throws IOException {
        fieldEncoder.writeFieldKey(encoder);
    }

    public void done(final BCF2Encoder encoder, final VariantContext vc) throws IOException { } // TODO -- overload done so that we null out values and test for correctness

    @Override
    public String toString() {
        return "BCF2FieldWriter " + getClass().getSimpleName() + " with encoder " + getFieldEncoder();
    }

    // --------------------------------------------------------------------------------
    //
    // Sites writers
    //
    // --------------------------------------------------------------------------------

    public static abstract class SiteWriter extends BCF2FieldWriter {
        protected SiteWriter(final VCFHeader header, final BCF2FieldEncoder fieldEncoder) {
            super(header, fieldEncoder);
        }

        public abstract void site(final BCF2Encoder encoder, final VariantContext vc) throws IOException;
    }

    public static class GenericSiteWriter extends SiteWriter {
        public GenericSiteWriter(final VCFHeader header, final BCF2FieldEncoder fieldEncoder) {
            super(header, fieldEncoder);
        }

        @Override
        public void site(final BCF2Encoder encoder, final VariantContext vc) throws IOException {
            final Object rawValue = vc.getAttribute(getField(), null);
            final BCF2Type type = getFieldEncoder().getType(rawValue);
            if ( rawValue == null ) {
                // the value is missing, just write in null
                encoder.encodeType(0, type);
            } else {
                final int valueCount = getFieldEncoder().numElements(vc, rawValue);
                encoder.encodeType(valueCount, type);
                getFieldEncoder().encodeValue(encoder, rawValue, type, valueCount);
            }
        }
    }

    // --------------------------------------------------------------------------------
    //
    // Genotypes writers
    //
    // --------------------------------------------------------------------------------

    public static abstract class GenotypesWriter extends BCF2FieldWriter {
        int nValuesPerGenotype = -1;
        BCF2Type encodingType = null;

        protected GenotypesWriter(final VCFHeader header, final BCF2FieldEncoder fieldEncoder) {
            super(header, fieldEncoder);

            if ( fieldEncoder.hasConstantNumElements() ) {
                nValuesPerGenotype = getFieldEncoder().numElements();
            }
        }

        @Override
        @Requires({"encodingType != null",
                "nValuesPerGenotype >= 0 || ! getFieldEncoder().hasConstantNumElements()"})
        @Ensures("nValuesPerGenotype >= 0")
        public void start(final BCF2Encoder encoder, final VariantContext vc) throws IOException {
            // writes the key information
            super.start(encoder, vc);

            // only update if we need to
            if ( ! getFieldEncoder().hasConstantNumElements() ) {
                if ( getFieldEncoder().hasContextDeterminedNumElements() )
                    // we are cheap -- just depends on genotype of allele counts
                    nValuesPerGenotype = getFieldEncoder().numElements(vc);
                else
                    // we have to go fishing through the values themselves (expensive)
                    nValuesPerGenotype = computeMaxSizeOfGenotypeFieldFromValues(vc);
            }

            encoder.encodeType(nValuesPerGenotype, encodingType);
        }

        @Requires({"encodingType != null", "nValuesPerGenotype >= 0"})
        public void addGenotype(final BCF2Encoder encoder, final VariantContext vc, final Genotype g) throws IOException {
            final Object fieldValue = g.getExtendedAttribute(getField(), null);
            getFieldEncoder().encodeValue(encoder, fieldValue, encodingType, nValuesPerGenotype);
        }

        @Ensures({"result >= 0"})
        protected int numElements(final VariantContext vc, final Genotype g) {
            return getFieldEncoder().numElements(vc, g.getExtendedAttribute(getField()));
        }

        @Ensures({"result >= 0"})
        private final int computeMaxSizeOfGenotypeFieldFromValues(final VariantContext vc) {
            int size = -1;

            for ( final Genotype g : vc.getGenotypes() ) {
                size = Math.max(size, numElements(vc, g));
            }

            return size;
        }
    }

    public static class StaticallyTypeGenotypesWriter extends GenotypesWriter {
        public StaticallyTypeGenotypesWriter(final VCFHeader header, final BCF2FieldEncoder fieldEncoder) {
            super(header, fieldEncoder);
            encodingType = getFieldEncoder().getStaticType();
        }
    }

    public static class IntegerTypeGenotypesWriter extends GenotypesWriter {
        public IntegerTypeGenotypesWriter(final VCFHeader header, final BCF2FieldEncoder fieldEncoder) {
            super(header, fieldEncoder);
        }

        @Override
        public void start(final BCF2Encoder encoder, final VariantContext vc) throws IOException {
            // the only value that is dynamic are integers
            final List<Integer> values = new ArrayList<Integer>(vc.getNSamples());
            for ( final Genotype g : vc.getGenotypes() ) {
                for ( final Object i : BCF2Utils.toList(g.getExtendedAttribute(getField(), null)) ) {
                    if ( i != null ) values.add((Integer)i); // we know they are all integers
                }
            }

            encodingType = BCF2Utils.determineIntegerType(values);
            super.start(encoder, vc);
        }
    }

    public static class IGFGenotypesWriter extends GenotypesWriter {
        final IntGenotypeFieldAccessors.Accessor ige;

        public IGFGenotypesWriter(final VCFHeader header, final BCF2FieldEncoder fieldEncoder, final IntGenotypeFieldAccessors.Accessor ige) {
            super(header, fieldEncoder);
            this.ige = ige;

            if ( ! (fieldEncoder instanceof BCF2FieldEncoder.IntArray) )
                throw new IllegalArgumentException("BUG: IntGenotypesWriter requires IntArray encoder for field " + getField());
        }

        @Override
        public void start(final BCF2Encoder encoder, final VariantContext vc) throws IOException {
            // TODO
            // TODO this piece of code consumes like 10% of the runtime alone because fo the vc.getGenotypes() iteration
            // TODO
            encodingType = BCF2Type.INT8;
            for ( final Genotype g : vc.getGenotypes() ) {
                final int[] pls = ige.getValues(g);
                final BCF2Type plsType = getFieldEncoder().getType(pls);
                encodingType = BCF2Utils.maxIntegerType(encodingType, plsType);
                if ( encodingType == BCF2Type.INT32 )
                    break; // stop early
            }

            super.start(encoder, vc);
        }

        @Override
        public void addGenotype(final BCF2Encoder encoder, final VariantContext vc, final Genotype g) throws IOException {
            getFieldEncoder().encodeValue(encoder, ige.getValues(g), encodingType, nValuesPerGenotype);
        }

        @Override
        protected int numElements(final VariantContext vc, final Genotype g) {
            return ige.getSize(g);
        }
    }

    public static class FTGenotypesWriter extends StaticallyTypeGenotypesWriter {
        public FTGenotypesWriter(final VCFHeader header, final BCF2FieldEncoder fieldEncoder) {
            super(header, fieldEncoder);
        }

        public void addGenotype(final BCF2Encoder encoder, final VariantContext vc, final Genotype g) throws IOException {
            final String fieldValue = g.getFilters();
            getFieldEncoder().encodeValue(encoder, fieldValue, encodingType, nValuesPerGenotype);
        }

        @Override
        protected int numElements(final VariantContext vc, final Genotype g) {
            return getFieldEncoder().numElements(vc, g.getFilters());
        }
    }

    public static class GTWriter extends GenotypesWriter {
        final Map<Allele, Integer> alleleMapForTriPlus = new HashMap<Allele, Integer>(5);
        Allele ref, alt1;

        public GTWriter(final VCFHeader header, final BCF2FieldEncoder fieldEncoder) {
            super(header, fieldEncoder);
        }

        @Override
        public void start(final BCF2Encoder encoder, final VariantContext vc) throws IOException {
            if ( vc.getNAlleles() > BCF2Utils.MAX_ALLELES_IN_GENOTYPES )
                throw new IllegalStateException("Current BCF2 encoder cannot handle sites " +
                        "with > " + BCF2Utils.MAX_ALLELES_IN_GENOTYPES + " alleles, but you have "
                        + vc.getNAlleles() + " at " + vc.getChr() + ":" + vc.getStart());

            encodingType = BCF2Type.INT8;
            buildAlleleMap(vc);
            nValuesPerGenotype = vc.getMaxPloidy(2);

            super.start(encoder, vc);
        }

        @Override
        public void addGenotype(final BCF2Encoder encoder, final VariantContext vc, final Genotype g) throws IOException {
            final int samplePloidy = g.getPloidy();
            for ( int i = 0; i < nValuesPerGenotype; i++ ) {
                if ( i < samplePloidy ) {
                    // we encode the actual allele
                    final Allele a = g.getAllele(i);
                    final int offset = getAlleleOffset(a);
                    final int encoded = ((offset+1) << 1) | (g.isPhased() ? 0x01 : 0x00);
                    encoder.encodeRawBytes(encoded, encodingType);
                } else {
                    // we need to pad with missing as we have ploidy < max for this sample
                    encoder.encodeRawBytes(encodingType.getMissingBytes(), encodingType);
                }
            }
        }

        /**
         * Fast path code to determine the offset.
         *
         * Inline tests for == against ref (most common, first test)
         * == alt1 (second most common, second test)
         * == NO_CALL (third)
         * and finally in the map from allele => offset for all alt 2+ alleles
         *
         * @param a the allele whose offset we wish to determine
         * @return the offset (from 0) of the allele in the list of variant context alleles (-1 means NO_CALL)
         */
        @Requires("a != null")
        private final int getAlleleOffset(final Allele a) {
            if ( a == ref ) return 0;
            else if ( a == alt1 ) return 1;
            else if ( a == Allele.NO_CALL ) return -1;
            else {
                final Integer o = alleleMapForTriPlus.get(a);
                if ( o == null ) throw new IllegalStateException("BUG: Couldn't find allele offset for allele " + a);
                return o;
            }
        }

        private final void buildAlleleMap(final VariantContext vc) {
            // these are fast path options to determine the offsets for
            final int nAlleles = vc.getNAlleles();
            ref = vc.getReference();
            alt1 = nAlleles > 1 ? vc.getAlternateAllele(0) : null;

            if ( nAlleles > 2 ) {
                // for multi-allelics we need to clear the map, and add additional looks
                alleleMapForTriPlus.clear();
                final List<Allele> alleles = vc.getAlleles();
                for ( int i = 2; i < alleles.size(); i++ ) {
                    alleleMapForTriPlus.put(alleles.get(i), i);
                }
            }
        }
    }
}

