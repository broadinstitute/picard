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

import org.broadinstitute.variant.VariantBaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Map;


/**
 * 
 * @author aaron 
 * 
 * Class VariantJEXLContextUnitTest
 *
 * Test out parts of the VariantJEXLContext
 */
public class VariantJEXLContextUnitTest extends VariantBaseTest {

    private static String expression = "QUAL > 500.0";
    private static VariantContextUtils.JexlVCMatchExp exp;

    Allele A, Aref, T, Tref;

    Allele ATC, ATCref;
    // A [ref] / T at 10

    // - / ATC [ref] from 20-23

    @BeforeClass
    public void beforeClass() {
        try {
            exp = new VariantContextUtils.JexlVCMatchExp("name", VariantContextUtils.engine.createExpression(expression));
        } catch (Exception e) {
            Assert.fail("Unable to create expression" + e.getMessage());
        }
    }

    @BeforeMethod
    public void before() {
        A = Allele.create("A");
        Aref = Allele.create("A", true);
        T = Allele.create("T");
        Tref = Allele.create("T", true);

        ATC = Allele.create("ATC");
        ATCref = Allele.create("ATC", true);
    }


    @Test
    public void testGetValue() {
        Map<VariantContextUtils.JexlVCMatchExp, Boolean> map = getVarContext();

        // make sure the context has a value
        Assert.assertTrue(!map.isEmpty());
        Assert.assertEquals(map.size(), 1);

        // eval our known expression
        Assert.assertTrue(!map.get(exp));
    }

    @Test(expectedExceptions=UnsupportedOperationException.class)
    public void testContainsValue() {
        Map<VariantContextUtils.JexlVCMatchExp, Boolean> map = getVarContext();

        map.containsValue(exp);
    }

    @Test(expectedExceptions=UnsupportedOperationException.class)
    public void testRemove() {
        Map<VariantContextUtils.JexlVCMatchExp, Boolean> map = getVarContext();

        map.remove(exp);
    }

    @Test(expectedExceptions=UnsupportedOperationException.class)
    public void testEntrySet() {
        Map<VariantContextUtils.JexlVCMatchExp, Boolean> map = getVarContext();

        map.entrySet();
    }

    @Test(expectedExceptions=UnsupportedOperationException.class)
    public void testClear() {
        Map<VariantContextUtils.JexlVCMatchExp, Boolean> map = getVarContext();

        map.clear();
    }

    /**
     * helper method
     * @return a VariantJEXLContext
     */
    private JEXLMap getVarContext() {
        List<Allele> alleles = Arrays.asList(Aref, T);

        VariantContext vc = new VariantContextBuilder("test", "chr1", 10, 10, alleles).make();
        return new JEXLMap(Arrays.asList(exp),vc);
    }
}
