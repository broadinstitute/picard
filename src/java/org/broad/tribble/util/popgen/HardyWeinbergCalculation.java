/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broad.tribble.util.popgen;

/**
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2004 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/

/**
* This class calculates a HardyWeinberg p-value given three values representing
* the observed frequences of homozygous and heterozygous genotypes in the
* test population.
*
* @author Bob Handsaker
*/
public final class HardyWeinbergCalculation {
   /**
    * This class is not instantiable.
    */
   private HardyWeinbergCalculation() {
   }

   /**
    * Calculates exact two-sided hardy-weinberg p-value. Parameters
    * are number of genotypes, number of rare alleles observed and
    * number of heterozygotes observed.
    *
    * (c) 2003 Jan Wigginton, Goncalo Abecasis (goncalo@umich.edu)
    */
   public static double hwCalculate(int obsAA, int obsAB, int obsBB) {
       int diplotypes = obsAA + obsAB + obsBB;
       int rare = (obsAA * 2) + obsAB;
       int hets = obsAB;

       //make sure "rare" allele is really the rare allele
       if (rare > diplotypes) {
           rare = (2 * diplotypes) - rare;
       }

       //make sure numbers aren't screwy
       if (hets > rare) {
           return -1;
       }

       double[] tailProbs = new double[rare + 1];

       for (int z = 0; z < tailProbs.length; z++) {
           tailProbs[z] = 0;
       }

       //start at midpoint
       int mid = (rare * ((2 * diplotypes) - rare)) / (2 * diplotypes);

       //check to ensure that midpoint and rare alleles have same parity
       if (((rare & 1) ^ (mid & 1)) != 0) {
           mid++;
       }

       int het = mid;
       int hom_r = (rare - mid) / 2;
       int hom_c = diplotypes - het - hom_r;

       //Calculate probability for each possible observed heterozygote
       //count up to a scaling constant, to avoid underflow and overflow
       tailProbs[mid] = 1.0;

       double sum = tailProbs[mid];

       for (het = mid; het > 1; het -= 2) {
           tailProbs[het - 2] = (tailProbs[het] * het * (het - 1.0)) / (4.0 * (hom_r + 1.0) * (hom_c +
               1.0));
           sum += tailProbs[het - 2];

           //2 fewer hets for next iteration -> add one rare and one common homozygote
           hom_r++;
           hom_c++;
       }

       het = mid;
       hom_r = (rare - mid) / 2;
       hom_c = diplotypes - het - hom_r;

       for (het = mid; het <= (rare - 2); het += 2) {
           tailProbs[het + 2] = (tailProbs[het] * 4.0 * hom_r * hom_c) / ((het + 2.0) * (het +
               1.0));
           sum += tailProbs[het + 2];

           //2 more hets for next iteration -> subtract one rare and one common homozygote
           hom_r--;
           hom_c--;
       }

       for (int z = 0; z < tailProbs.length; z++) {
           tailProbs[z] /= sum;
       }

       double top = tailProbs[hets];

       for (int i = hets + 1; i <= rare; i++) {
           top += tailProbs[i];
       }

       double otherSide = tailProbs[hets];

       for (int i = hets - 1; i >= 0; i--) {
           otherSide += tailProbs[i];
       }

       if ((top > 0.5) && (otherSide > 0.5)) {
           return 1.0;
       }

       if (top < otherSide) {
           return top * 2;
       }

       return otherSide * 2;
   }
}