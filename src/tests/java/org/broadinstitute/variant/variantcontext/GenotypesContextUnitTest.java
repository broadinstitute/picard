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


import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.variant.VariantBaseTest;
import org.broadinstitute.variant.utils.GeneralUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;


public class GenotypesContextUnitTest extends VariantBaseTest {
    Allele Aref, C, T;
    Genotype AA, AT, TT, AC, CT, CC, MISSING;
    List<Genotype> allGenotypes;

    @BeforeSuite
    public void before() {
        C = Allele.create("C");
        Aref = Allele.create("A", true);
        T = Allele.create("T");
        AA = GenotypeBuilder.create("AA", Arrays.asList(Aref, Aref));
        AT = GenotypeBuilder.create("AT", Arrays.asList(Aref, T));
        TT = GenotypeBuilder.create("TT", Arrays.asList(T, T));
        AC = GenotypeBuilder.create("AC", Arrays.asList(Aref, C));
        CT = GenotypeBuilder.create("CT", Arrays.asList(C, T));
        CC = GenotypeBuilder.create("CC", Arrays.asList(C, C));
        MISSING = GenotypeBuilder.create("MISSING", Arrays.asList(C, C));

        allGenotypes = Arrays.asList(AA, AT, TT, AC, CT, CC);
    }

    // --------------------------------------------------------------------------------
    //
    // Provider
    //
    // --------------------------------------------------------------------------------

    private interface ContextMaker {
        public GenotypesContext make(List<Genotype> initialSamples);
    }

    private ContextMaker baseMaker = new ContextMaker() {
        @Override
        public GenotypesContext make(final List<Genotype> initialSamples) {
            return GenotypesContext.copy(initialSamples);
        }

        @Override
        public String toString() {
            return "GenotypesContext";
        }
    };

    private final class lazyMaker implements LazyGenotypesContext.LazyParser, ContextMaker {
        @Override
        public LazyGenotypesContext.LazyData parse(final Object data) {
            GenotypesContext gc = GenotypesContext.copy((List<Genotype>)data);
            gc.ensureSampleNameMap();
            gc.ensureSampleOrdering();
            return new LazyGenotypesContext.LazyData(gc.notToBeDirectlyAccessedGenotypes, gc.sampleNamesInOrder, gc.sampleNameToOffset);
        }

        @Override
        public GenotypesContext make(final List<Genotype> initialSamples) {
            return new LazyGenotypesContext(this, initialSamples, initialSamples.size());
        }

        @Override
        public String toString() {
            return "LazyGenotypesContext";
        }
    }

    private Collection<ContextMaker> allMakers = Arrays.asList(baseMaker, new lazyMaker());

    private class GenotypesContextProvider {
        String name;
        ContextMaker maker;
        final List<Genotype> initialSamples;

        private GenotypesContextProvider(ContextMaker maker, List<Genotype> initialSamples) {
            this.name = String.format("%s with %d samples", maker.toString(), initialSamples.size());
            this.maker = maker;
            this.initialSamples = initialSamples;
        }

        public GenotypesContext makeContext() {
            return maker.make(initialSamples);
        }
    }

    @DataProvider(name = "GenotypesContextProvider")
    public Object[][] MakeSampleNamesTest() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( ContextMaker maker : allMakers ) {
            for ( int i = 0; i < allGenotypes.size(); i++ ) {
                List<Genotype> samples = allGenotypes.subList(0, i);
                // sorted
                tests.add(new Object[]{new GenotypesContextProvider(maker, samples)});
                // unsorted
                tests.add(new Object[]{new GenotypesContextProvider(maker, GeneralUtils.reverse(samples))});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    private final static void testIterable(Iterable<Genotype> genotypeIterable, Set<String> expectedNames) {
        int count = 0;
        for ( final Genotype g : genotypeIterable ) {
            Assert.assertTrue(expectedNames.contains(g.getSampleName()));
            count++;
        }
        Assert.assertEquals(count, expectedNames.size(), "Iterable returned unexpected number of genotypes");
    }

    @Test(dataProvider = "GenotypesContextProvider")
    public void testInitialSamplesAreAsExpected(GenotypesContextProvider cfg) {
        testGenotypesContextContainsExpectedSamples(cfg.makeContext(), cfg.initialSamples);
    }

    private final void testGenotypesContextContainsExpectedSamples(GenotypesContext gc, List<Genotype> expectedSamples) {
        Assert.assertEquals(gc.isEmpty(), expectedSamples.isEmpty());
        Assert.assertEquals(gc.size(), expectedSamples.size());

        // get(index) is doing the right thing
        for ( int i = 0; i < expectedSamples.size(); i++ ) {
            Assert.assertEquals(gc.get(i), expectedSamples.get(i));
        }
        Assert.assertFalse(gc.containsSample(MISSING.getSampleName()));

        // we can fetch samples by name
        final Set<String> genotypeNames = VariantContextUtils.genotypeNames(expectedSamples);
        for ( final String name : genotypeNames ) {
            Assert.assertTrue(gc.containsSample(name));
        }
        Assert.assertFalse(gc.containsSample(MISSING.getSampleName()));

        // all of the iterators are working
        testIterable(gc.iterateInSampleNameOrder(), genotypeNames);
        testIterable(gc, genotypeNames);
        testIterable(gc.iterateInSampleNameOrder(genotypeNames), genotypeNames);
        if ( ! genotypeNames.isEmpty() ) {
            Set<String> first = Collections.singleton(genotypeNames.iterator().next());
            testIterable(gc.iterateInSampleNameOrder(first), first);
        }

        // misc. utils are working as expected
        assertEqualsSet(gc.getSampleNames(), genotypeNames, "gc sample names vs. expected sample names");
        Assert.assertTrue(ParsingUtils.isSorted(gc.getSampleNamesOrderedByName()));
        Assert.assertTrue(ParsingUtils.isSorted(gc.iterateInSampleNameOrder()));
        Assert.assertTrue(gc.containsSamples(genotypeNames));

        final Set<String> withMissing = new HashSet<String>(Arrays.asList(MISSING.getSampleName()));
        withMissing.addAll(genotypeNames);
        Assert.assertFalse(gc.containsSamples(withMissing));
    }

    @Test(enabled = true, dataProvider = "GenotypesContextProvider")
    public void testImmutable(GenotypesContextProvider cfg) {
        GenotypesContext gc = cfg.makeContext();
        Assert.assertEquals(gc.isMutable(), true);
        gc.immutable();
        Assert.assertEquals(gc.isMutable(), false);
    }

    @Test(enabled = true, dataProvider = "GenotypesContextProvider", expectedExceptions = Throwable.class )
    public void testImmutableCall1(GenotypesContextProvider cfg) {
        GenotypesContext gc = cfg.makeContext();
        gc.immutable();
        gc.add(MISSING);
    }

    @Test(enabled = true, dataProvider = "GenotypesContextProvider")
    public void testClear(GenotypesContextProvider cfg) {
        GenotypesContext gc = cfg.makeContext();
        gc.clear();
        testGenotypesContextContainsExpectedSamples(gc, Collections.<Genotype>emptyList());
    }

    private static final List<Genotype> with(List<Genotype> genotypes, Genotype ... add) {
        List<Genotype> l = new ArrayList<Genotype>(genotypes);
        l.addAll(Arrays.asList(add));
        return l;
    }

    private static final List<Genotype> without(List<Genotype> genotypes, Genotype ... remove) {
        List<Genotype> l = new ArrayList<Genotype>(genotypes);
        l.removeAll(Arrays.asList(remove));
        return l;
    }

    @Test(enabled = true, dataProvider = "GenotypesContextProvider")
    public void testAdds(GenotypesContextProvider cfg) {
        Genotype add1 = GenotypeBuilder.create("add1", Arrays.asList(Aref, Aref));
        Genotype add2 = GenotypeBuilder.create("add2", Arrays.asList(Aref, Aref));

        GenotypesContext gc = cfg.makeContext();
        gc.add(add1);
        testGenotypesContextContainsExpectedSamples(gc, with(cfg.initialSamples, add1));

        gc = cfg.makeContext();
        gc.add(add1);
        gc.add(add2);
        testGenotypesContextContainsExpectedSamples(gc, with(cfg.initialSamples, add1, add2));

        gc = cfg.makeContext();
        gc.addAll(Arrays.asList(add1, add2));
        testGenotypesContextContainsExpectedSamples(gc, with(cfg.initialSamples, add1, add2));
    }

    @Test(enabled = true, dataProvider = "GenotypesContextProvider")
    public void testRemoves(GenotypesContextProvider cfg) {
        Genotype rm1 = AA;
        Genotype rm2 = AC;

        GenotypesContext gc = cfg.makeContext();
        if (gc.size() > 1) {
            Genotype rm = gc.get(0);
            gc.remove(rm);
            testGenotypesContextContainsExpectedSamples(gc, without(cfg.initialSamples, rm));
        }

        gc = cfg.makeContext();
        gc.remove(rm1);
        testGenotypesContextContainsExpectedSamples(gc, without(cfg.initialSamples, rm1));

        gc = cfg.makeContext();
        gc.remove(rm1);
        gc.remove(rm2);
        testGenotypesContextContainsExpectedSamples(gc, without(cfg.initialSamples, rm1, rm2));

        gc = cfg.makeContext();
        gc.removeAll(Arrays.asList(rm1, rm2));
        testGenotypesContextContainsExpectedSamples(gc, without(cfg.initialSamples, rm1, rm2));

        gc = cfg.makeContext();
        HashSet<Genotype> expected = new HashSet<Genotype>();
        if ( gc.contains(rm1) ) expected.add(rm1);
        if ( gc.contains(rm2) ) expected.add(rm2);
        gc.retainAll(Arrays.asList(rm1, rm2));

        // ensure that the two lists are the same
        assertEqualsSet(new HashSet<Genotype>(gc.getGenotypes()), expected, "gc genotypes vs. expected");
        // because the list order can change, we use the gc's list itself
        testGenotypesContextContainsExpectedSamples(gc, gc.getGenotypes());
    }

    @Test(enabled = true, dataProvider = "GenotypesContextProvider")
    public void testSet(GenotypesContextProvider cfg) {
        Genotype set = GenotypeBuilder.create("replace", Arrays.asList(Aref, Aref));
        int n = cfg.makeContext().size();
        for ( int i = 0; i < n; i++ ) {
            GenotypesContext gc = cfg.makeContext();
            Genotype setted = gc.set(i, set);
            Assert.assertNotNull(setted);
            ArrayList<Genotype> l = new ArrayList<Genotype>(cfg.initialSamples);
            l.set(i, set);
            testGenotypesContextContainsExpectedSamples(gc, l);
        }
    }

    @Test(enabled = true, dataProvider = "GenotypesContextProvider")
    public void testReplace(GenotypesContextProvider cfg) {
        int n = cfg.makeContext().size();
        for ( int i = 0; i < n; i++ ) {
            GenotypesContext gc = cfg.makeContext();
            Genotype toReplace = gc.get(i);
            Genotype replacement = GenotypeBuilder.create(toReplace.getSampleName(), Arrays.asList(Aref, Aref));
            gc.replace(replacement);
            ArrayList<Genotype> l = new ArrayList<Genotype>(cfg.initialSamples);
            l.set(i, replacement);
            Assert.assertEquals(replacement, gc.get(i));
            testGenotypesContextContainsExpectedSamples(gc, l);
        }
    }

    // subset to samples tested in VariantContextUnitTest
}
