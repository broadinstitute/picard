package picard.vcf.concordance;

/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/**
 * Enum class that lists the possible variant call states for a call (a called variant, as opposed to the truth variant)
 *
 * @author George Grant
 * */
public enum CalledVariantCallState {
    HomRef,
    Het,                // call and truth alleles agree.  call is A*C and A*C
    Het1Mismatch,       // called variant is a het, but differs from truth variant (het) in one allele (i.e. Call is A*C and Truth is A*T)
    Het2Mismatch,       // called variant is a het, but differs from truth variant (het) in one allele (i.e. Call is GC and Truth is A*T)
    HomVar,             // call and truth agree.  Call is CC, truth is CC
    HomVarMismatch,     // called variant is HomVar, but differs from truth variant (HomVar) in alleles (i.e. Call is AA and Truth is TT)
    NoCall,
    Filtered
}

