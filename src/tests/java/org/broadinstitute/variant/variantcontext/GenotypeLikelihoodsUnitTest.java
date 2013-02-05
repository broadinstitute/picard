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


// the imports for unit testing.


import org.broad.tribble.TribbleException;
import org.broadinstitute.variant.VariantBaseTest;
import org.broadinstitute.variant.utils.GeneralUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.EnumMap;
import java.util.List;


/**
 * Basic unit test for Genotype likelihoods objects
 */
public class GenotypeLikelihoodsUnitTest extends VariantBaseTest {
    double [] v = new double[]{-10.5, -1.25, -5.11};
    final static String vGLString = "-10.50,-1.25,-5.11";
    final static String vPLString = "93,0,39";
    double[] triAllelic = new double[]{-4.2,-2.0,-3.0,-1.6,0.0,-4.0}; //AA,AB,AC,BB,BC,CC

    @Test
    public void testFromVector2() {
        GenotypeLikelihoods gl = GenotypeLikelihoods.fromLog10Likelihoods(v);
        assertDoubleArraysAreEqual(gl.getAsVector(), v);
        Assert.assertEquals(gl.getAsString(), vPLString);
    }

    @Test
    public void testFromString1() {
        GenotypeLikelihoods gl = GenotypeLikelihoods.fromPLField(vPLString);
        assertDoubleArraysAreEqual(gl.getAsVector(), new double[]{-9.3, 0, -3.9});
        Assert.assertEquals(gl.getAsString(), vPLString);
    }

    @Test
    public void testFromString2() {
        GenotypeLikelihoods gl = GenotypeLikelihoods.fromGLField(vGLString);
        assertDoubleArraysAreEqual(gl.getAsVector(), v);
        Assert.assertEquals(gl.getAsString(), vPLString);
    }

    @Test (expectedExceptions = TribbleException.class)
    public void testErrorBadFormat() {
        GenotypeLikelihoods gl = GenotypeLikelihoods.fromPLField("adf,b,c");
        gl.getAsVector();
    }

    @Test
    public void testGetAsMap(){
        GenotypeLikelihoods gl = GenotypeLikelihoods.fromLog10Likelihoods(v);
        //Log scale
        EnumMap<GenotypeType,Double> glMap = gl.getAsMap(false);
        Assert.assertEquals(v[GenotypeType.HOM_REF.ordinal()-1],glMap.get(GenotypeType.HOM_REF));
        Assert.assertEquals(v[GenotypeType.HET.ordinal()-1],glMap.get(GenotypeType.HET));
        Assert.assertEquals(v[GenotypeType.HOM_VAR.ordinal()-1],glMap.get(GenotypeType.HOM_VAR));

        //Linear scale
        glMap = gl.getAsMap(true);
        double [] vl = GeneralUtils.normalizeFromLog10(v);
        Assert.assertEquals(vl[GenotypeType.HOM_REF.ordinal()-1],glMap.get(GenotypeType.HOM_REF));
        Assert.assertEquals(vl[GenotypeType.HET.ordinal()-1],glMap.get(GenotypeType.HET));
        Assert.assertEquals(vl[GenotypeType.HOM_VAR.ordinal()-1],glMap.get(GenotypeType.HOM_VAR));

        //Test missing likelihoods
        gl = GenotypeLikelihoods.fromPLField(".");
        glMap = gl.getAsMap(false);
        Assert.assertNull(glMap);

    }

    @Test
    public void testCalculateNumLikelihoods() {    
        
        for (int nAlleles=2; nAlleles<=5; nAlleles++)
            // simplest case: diploid
            Assert.assertEquals(GenotypeLikelihoods.numLikelihoods(nAlleles, 2), nAlleles*(nAlleles+1)/2);

        // some special cases: ploidy = 20, #alleles = 4
        Assert.assertEquals(GenotypeLikelihoods.numLikelihoods(4, 20), 1771);
    }
    
    @Test
    public void testGetLog10GQ(){
        GenotypeLikelihoods gl = GenotypeLikelihoods.fromPLField(vPLString);

        //GQ for the best guess genotype
        Assert.assertEquals(gl.getLog10GQ(GenotypeType.HET),-3.9);

        double[] test = GeneralUtils.normalizeFromLog10(gl.getAsVector());

        //GQ for the other genotypes
        Assert.assertEquals(gl.getLog10GQ(GenotypeType.HOM_REF), Math.log10(1.0 - test[GenotypeType.HOM_REF.ordinal()-1]));
        Assert.assertEquals(gl.getLog10GQ(GenotypeType.HOM_VAR), Math.log10(1.0 - test[GenotypeType.HOM_VAR.ordinal()-1]));

       //Test missing likelihoods
        gl = GenotypeLikelihoods.fromPLField(".");
        Assert.assertEquals(gl.getLog10GQ(GenotypeType.HOM_REF),Double.NEGATIVE_INFINITY);
        Assert.assertEquals(gl.getLog10GQ(GenotypeType.HET),Double.NEGATIVE_INFINITY);
        Assert.assertEquals(gl.getLog10GQ(GenotypeType.HOM_VAR),Double.NEGATIVE_INFINITY);

    }

    @Test
    public void testgetQualFromLikelihoods() {
        double[] likelihoods = new double[]{-1, 0, -2};
        // qual values we expect for each possible "best" genotype
        double[] expectedQuals = new double[]{-0.04100161, -1, -0.003930294};

        for ( int i = 0; i < likelihoods.length; i++ ) {
            Assert.assertEquals(GenotypeLikelihoods.getGQLog10FromLikelihoods(i, likelihoods), expectedQuals[i], 1e-6,
                    "GQ value for genotype " + i + " was not calculated correctly");
        }
    }

    // this test is completely broken, the method is wrong.
    public void testGetQualFromLikelihoodsMultiAllelicBroken() {
        GenotypeLikelihoods gl = GenotypeLikelihoods.fromLog10Likelihoods(triAllelic);
        double actualGQ = gl.getLog10GQ(GenotypeType.HET);
        double expectedGQ = 1.6;
        Assert.assertEquals(actualGQ,expectedGQ);
    }

    public void testGetQualFromLikelihoodsMultiAllelic() {
        GenotypeLikelihoods gl = GenotypeLikelihoods.fromLog10Likelihoods(triAllelic);
        Allele ref = Allele.create((byte)'A',true);
        Allele alt1 = Allele.create((byte)'C');
        Allele alt2 = Allele.create((byte)'T');
        List<Allele> allAlleles = Arrays.asList(ref,alt1,alt2);
        List<Allele> gtAlleles = Arrays.asList(alt1,alt2);
        GenotypeBuilder gtBuilder = new GenotypeBuilder();
        gtBuilder.alleles(gtAlleles);
        double actualGQ = gl.getLog10GQ(gtBuilder.make(),allAlleles);
        double expectedGQ = 1.6;
        Assert.assertEquals(actualGQ,expectedGQ);
    }

    private void assertDoubleArraysAreEqual(double[] v1, double[] v2) {
        Assert.assertEquals(v1.length, v2.length);
        for ( int i = 0; i < v1.length; i++ ) {
            Assert.assertEquals(v1[i], v2[i], 1e-6);
        }
    }

    @Test
    public void testCalculatePLindex(){
        int counter = 0;
        for ( int i = 0; i <= 3; i++ ) {
            for ( int j = i; j <= 3; j++ ) {
                Assert.assertEquals(GenotypeLikelihoods.calculatePLindex(i, j), GenotypeLikelihoods.PLindexConversion[counter++], "PL index of alleles " + i + "," + j + " was not calculated correctly");
            }
        }
    }

    @Test
    public void testGetAllelePair(){
        allelePairTest(0, 0, 0);
        allelePairTest(1, 0, 1);
        allelePairTest(2, 1, 1);
        allelePairTest(3, 0, 2);
        allelePairTest(4, 1, 2);
        allelePairTest(5, 2, 2);
        allelePairTest(6, 0, 3);
        allelePairTest(7, 1, 3);
        allelePairTest(8, 2, 3);
        allelePairTest(9, 3, 3);
    }
        
    private void allelePairTest(int PLindex, int allele1, int allele2) {
        Assert.assertEquals(GenotypeLikelihoods.getAllelePair(PLindex).alleleIndex1, allele1, "allele index " + allele1 + " from PL index " + PLindex + " was not calculated correctly");
        Assert.assertEquals(GenotypeLikelihoods.getAllelePair(PLindex).alleleIndex2, allele2, "allele index " + allele2 + " from PL index " + PLindex + " was not calculated correctly");
    }
}