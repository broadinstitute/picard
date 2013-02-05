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


import org.broadinstitute.variant.VariantBaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;


public class GenotypeUnitTest extends VariantBaseTest {
    Allele A, Aref, T;

    @BeforeSuite
    public void before() {
        A = Allele.create("A");
        Aref = Allele.create("A", true);
        T = Allele.create("T");
    }

    private static final GenotypeBuilder makeGB() {
        return new GenotypeBuilder("misc");
    }

    @Test
    public void testFilters() {
        Assert.assertFalse(makeGB().make().isFiltered(), "by default Genotypes must be PASS");
        Assert.assertNull(makeGB().make().getFilters(), "by default Genotypes must be PASS => getFilters() == null");
        Assert.assertFalse(makeGB().filter(null).make().isFiltered(), "setting filter == null => Genotypes must be PASS");
        Assert.assertNull(makeGB().filter(null).make().getFilters(), "Genotypes PASS => getFilters == null");
        Assert.assertFalse(makeGB().filter("PASS").make().isFiltered(), "setting filter == PASS => Genotypes must be PASS");
        Assert.assertNull(makeGB().filter("PASS").make().getFilters(), "Genotypes PASS => getFilters == null");
        Assert.assertTrue(makeGB().filter("x").make().isFiltered(), "setting filter != null => Genotypes must be PASS");
        Assert.assertEquals(makeGB().filter("x").make().getFilters(), "x", "Should get back the expected filter string");
        Assert.assertEquals(makeGB().filters("x", "y").make().getFilters(), "x;y", "Multiple filter field values should be joined with ;");
        Assert.assertEquals(makeGB().filters("x", "y", "z").make().getFilters(), "x;y;z", "Multiple filter field values should be joined with ;");
        Assert.assertTrue(makeGB().filters("x", "y", "z").make().isFiltered(), "Multiple filter values should be filtered");
        Assert.assertEquals(makeGB().filter("x;y;z").make().getFilters(), "x;y;z", "Multiple filter field values should be joined with ;");
    }

//    public Genotype(String sampleName, List<Allele> alleles, double negLog10PError, Set<String> filters, Map<String, ?> attributes, boolean isPhased) {
//    public Genotype(String sampleName, List<Allele> alleles, double negLog10PError, Set<String> filters, Map<String, ?> attributes, boolean isPhased, double[] log10Likelihoods) {
//    public Genotype(String sampleName, List<Allele> alleles, double negLog10PError, double[] log10Likelihoods)
//    public Genotype(String sampleName, List<Allele> alleles, double negLog10PError)
//    public Genotype(String sampleName, List<Allele> alleles)
//    public List<Allele> getAlleles()
//    public List<Allele> getAlleles(Allele allele)
//    public Allele getAllele(int i)
//    public boolean isPhased()
//    public int getPloidy()
//    public Type getType()
//    public boolean isHom()
//    public boolean isHomRef()
//    public boolean isHomVar()
//    public boolean isHet()
//    public boolean isNoCall()
//    public boolean isCalled()
//    public boolean isAvailable()
//    public boolean hasLikelihoods()
//    public GenotypeLikelihoods getLikelihoods()
//    public boolean sameGenotype(Genotype other)
//    public boolean sameGenotype(Genotype other, boolean ignorePhase)
//    public String getSampleName()
//    public boolean hasLog10PError()
//    public double getLog10PError()
//    public double getPhredScaledQual()
//    public boolean hasExtendedAttribute(String key)
//    public Object getExtendedAttribute(String key)
//    public Object getExtendedAttribute(String key, Object defaultValue)
//    public String getAttributeAsString(String key, String defaultValue)
//    public int getAttributeAsInt(String key, int defaultValue)
//    public double getAttributeAsDouble(String key, double  defaultValue)
//    public boolean getAttributeAsBoolean(String key, boolean  defaultValue)
}
