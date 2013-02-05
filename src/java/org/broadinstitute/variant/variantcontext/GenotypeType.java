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

package org.broadinstitute.variant.variantcontext;

/**
 * Summary types for Genotype objects
 *
 * @author Your Name
 * @since Date created
 */
public enum GenotypeType {
    /** The sample is no-called (all alleles are NO_CALL */
    NO_CALL,
    /** The sample is homozygous reference */
    HOM_REF,
    /** The sample is heterozygous, with at least one ref and at least one one alt in any order */
    HET,
    /** All alleles are non-reference */
    HOM_VAR,
    /** There is no allele data availble for this sample (alleles.isEmpty) */
    UNAVAILABLE,
    /** Some chromosomes are NO_CALL and others are called */
    MIXED  // no-call and call in the same genotype
}
