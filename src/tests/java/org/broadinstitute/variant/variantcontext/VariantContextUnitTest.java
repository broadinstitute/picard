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
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.util.*;


public class VariantContextUnitTest extends VariantBaseTest {
    Allele A, Aref, C, T, Tref;
    Allele del, delRef, ATC, ATCref;

    // A [ref] / T at 10
    String snpLoc = "chr1";
    int snpLocStart = 10;
    int snpLocStop = 10;

    // - / ATC [ref] from 20-22
    String delLoc = "chr1";
    int delLocStart = 20;
    int delLocStop = 22;

    // - [ref] / ATC from 20-20
    String insLoc = "chr1";
    int insLocStart = 20;
    int insLocStop = 20;

    VariantContextBuilder basicBuilder, snpBuilder, insBuilder;

    @BeforeSuite
    public void before() {
        del = Allele.create("A");
        delRef = Allele.create("A", true);

        A = Allele.create("A");
        C = Allele.create("C");
        Aref = Allele.create("A", true);
        T = Allele.create("T");
        Tref = Allele.create("T", true);

        ATC = Allele.create("ATC");
        ATCref = Allele.create("ATC", true);
    }

    @BeforeMethod
    public void beforeTest() {
        basicBuilder = new VariantContextBuilder("test", snpLoc,snpLocStart, snpLocStop, Arrays.asList(Aref, T));
        snpBuilder = new VariantContextBuilder("test", snpLoc,snpLocStart, snpLocStop, Arrays.asList(Aref, T));
        insBuilder = new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, Arrays.asList(delRef, ATC));
    }

    @Test
    public void testDetermineTypes() {
        Allele ACref = Allele.create("AC", true);
        Allele AC = Allele.create("AC");
        Allele AT = Allele.create("AT");
        Allele C = Allele.create("C");
        Allele CAT = Allele.create("CAT");
        Allele TAref = Allele.create("TA", true);
        Allele TA = Allele.create("TA");
        Allele TC = Allele.create("TC");
        Allele symbolic = Allele.create("<FOO>");

        // test REF
        List<Allele> alleles = Arrays.asList(Tref);
        VariantContext vc = snpBuilder.alleles(alleles).stop(snpLocStop).make();
        Assert.assertEquals(vc.getType(), VariantContext.Type.NO_VARIATION);

        // test SNPs
        alleles = Arrays.asList(Tref, A);
        vc = snpBuilder.alleles(alleles).stop(snpLocStop).make();
        Assert.assertEquals(vc.getType(), VariantContext.Type.SNP);

        alleles = Arrays.asList(Tref, A, C);
        vc = snpBuilder.alleles(alleles).stop(snpLocStop).make();
        Assert.assertEquals(vc.getType(), VariantContext.Type.SNP);

        // test MNPs
        alleles = Arrays.asList(ACref, TA);
        vc = snpBuilder.alleles(alleles).stop(snpLocStop+1).make();
        Assert.assertEquals(vc.getType(), VariantContext.Type.MNP);

        alleles = Arrays.asList(ATCref, CAT, Allele.create("GGG"));
        vc = basicBuilder.alleles(alleles).stop(snpLocStop+2).make();
        Assert.assertEquals(vc.getType(), VariantContext.Type.MNP);

        // test INDELs
        alleles = Arrays.asList(Aref, ATC);
        vc = basicBuilder.alleles(alleles).stop(snpLocStop).make();
        Assert.assertEquals(vc.getType(), VariantContext.Type.INDEL);

        alleles = Arrays.asList(ATCref, A);
        vc = basicBuilder.alleles(alleles).stop(snpLocStop+2).make();
        Assert.assertEquals(vc.getType(), VariantContext.Type.INDEL);

        alleles = Arrays.asList(Tref, TA, TC);
        vc = basicBuilder.alleles(alleles).stop(snpLocStop).make();
        Assert.assertEquals(vc.getType(), VariantContext.Type.INDEL);

        alleles = Arrays.asList(ATCref, A, AC);
        vc = basicBuilder.alleles(alleles).stop(snpLocStop+2).make();
        Assert.assertEquals(vc.getType(), VariantContext.Type.INDEL);

        alleles = Arrays.asList(ATCref, A, Allele.create("ATCTC"));
        vc = basicBuilder.alleles(alleles).stop(snpLocStop+2).make();
        Assert.assertEquals(vc.getType(), VariantContext.Type.INDEL);

        // test MIXED
        alleles = Arrays.asList(TAref, T, TC);
        vc = basicBuilder.alleles(alleles).stop(snpLocStop+1).make();
        Assert.assertEquals(vc.getType(), VariantContext.Type.MIXED);

        alleles = Arrays.asList(TAref, T, AC);
        vc = basicBuilder.alleles(alleles).stop(snpLocStop+1).make();
        Assert.assertEquals(vc.getType(), VariantContext.Type.MIXED);

        alleles = Arrays.asList(ACref, ATC, AT);
        vc = basicBuilder.alleles(alleles).stop(snpLocStop+1).make();
        Assert.assertEquals(vc.getType(), VariantContext.Type.MIXED);

        alleles = Arrays.asList(Aref, T, symbolic);
        vc = basicBuilder.alleles(alleles).stop(snpLocStop).make();
        Assert.assertEquals(vc.getType(), VariantContext.Type.MIXED);

        // test SYMBOLIC
        alleles = Arrays.asList(Tref, symbolic);
        vc = basicBuilder.alleles(alleles).stop(snpLocStop).make();
        Assert.assertEquals(vc.getType(), VariantContext.Type.SYMBOLIC);
    }

    @Test
    public void testMultipleSNPAlleleOrdering() {
        final List<Allele> allelesNaturalOrder = Arrays.asList(Aref, C, T);
        final List<Allele> allelesUnnaturalOrder = Arrays.asList(Aref, T, C);
        VariantContext naturalVC = snpBuilder.alleles(allelesNaturalOrder).make();
        VariantContext unnaturalVC = snpBuilder.alleles(allelesUnnaturalOrder).make();
        Assert.assertEquals(new ArrayList<Allele>(naturalVC.getAlleles()), allelesNaturalOrder);
        Assert.assertEquals(new ArrayList<Allele>(unnaturalVC.getAlleles()), allelesUnnaturalOrder);
    }

    @Test
    public void testCreatingSNPVariantContext() {

        List<Allele> alleles = Arrays.asList(Aref, T);
        VariantContext vc = snpBuilder.alleles(alleles).make();

        Assert.assertEquals(vc.getChr(), snpLoc);
        Assert.assertEquals(vc.getStart(), snpLocStart);
        Assert.assertEquals(vc.getEnd(), snpLocStop);
        Assert.assertEquals(vc.getType(), VariantContext.Type.SNP);
        Assert.assertTrue(vc.isSNP());
        Assert.assertFalse(vc.isIndel());
        Assert.assertFalse(vc.isSimpleInsertion());
        Assert.assertFalse(vc.isSimpleDeletion());
        Assert.assertFalse(vc.isMixed());
        Assert.assertTrue(vc.isBiallelic());
        Assert.assertEquals(vc.getNAlleles(), 2);

        Assert.assertEquals(vc.getReference(), Aref);
        Assert.assertEquals(vc.getAlleles().size(), 2);
        Assert.assertEquals(vc.getAlternateAlleles().size(), 1);
        Assert.assertEquals(vc.getAlternateAllele(0), T);

        Assert.assertFalse(vc.hasGenotypes());

        Assert.assertEquals(vc.getSampleNames().size(), 0);
    }

    @Test
    public void testCreatingRefVariantContext() {
        List<Allele> alleles = Arrays.asList(Aref);
        VariantContext vc = snpBuilder.alleles(alleles).make();

        Assert.assertEquals(vc.getChr(), snpLoc);
        Assert.assertEquals(vc.getStart(), snpLocStart);
        Assert.assertEquals(vc.getEnd(), snpLocStop);
        Assert.assertEquals(VariantContext.Type.NO_VARIATION, vc.getType());
        Assert.assertFalse(vc.isSNP());
        Assert.assertFalse(vc.isIndel());
        Assert.assertFalse(vc.isSimpleInsertion());
        Assert.assertFalse(vc.isSimpleDeletion());
        Assert.assertFalse(vc.isMixed());
        Assert.assertFalse(vc.isBiallelic());
        Assert.assertEquals(vc.getNAlleles(), 1);

        Assert.assertEquals(vc.getReference(), Aref);
        Assert.assertEquals(vc.getAlleles().size(), 1);
        Assert.assertEquals(vc.getAlternateAlleles().size(), 0);
        //Assert.assertEquals(vc.getAlternateAllele(0), T);

        Assert.assertFalse(vc.hasGenotypes());
        Assert.assertEquals(vc.getSampleNames().size(), 0);
    }

    @Test
    public void testCreatingDeletionVariantContext() {
        List<Allele> alleles = Arrays.asList(ATCref, del);
        VariantContext vc = new VariantContextBuilder("test", delLoc, delLocStart, delLocStop, alleles).make();

        Assert.assertEquals(vc.getChr(), delLoc);
        Assert.assertEquals(vc.getStart(), delLocStart);
        Assert.assertEquals(vc.getEnd(), delLocStop);
        Assert.assertEquals(vc.getType(), VariantContext.Type.INDEL);
        Assert.assertFalse(vc.isSNP());
        Assert.assertTrue(vc.isIndel());
        Assert.assertFalse(vc.isSimpleInsertion());
        Assert.assertTrue(vc.isSimpleDeletion());
        Assert.assertFalse(vc.isMixed());
        Assert.assertTrue(vc.isBiallelic());
        Assert.assertEquals(vc.getNAlleles(), 2);

        Assert.assertEquals(vc.getReference(), ATCref);
        Assert.assertEquals(vc.getAlleles().size(), 2);
        Assert.assertEquals(vc.getAlternateAlleles().size(), 1);
        Assert.assertEquals(vc.getAlternateAllele(0), del);

        Assert.assertFalse(vc.hasGenotypes());

        Assert.assertEquals(vc.getSampleNames().size(), 0);
    }

    @Test
    public void testMatchingAlleles() {
        List<Allele> alleles = Arrays.asList(ATCref, del);
        VariantContext vc = new VariantContextBuilder("test", delLoc, delLocStart, delLocStop, alleles).make();
        VariantContext vc2 = new VariantContextBuilder("test2", delLoc, delLocStart+12, delLocStop+12, alleles).make();

        Assert.assertTrue(vc.hasSameAllelesAs(vc2));
        Assert.assertTrue(vc.hasSameAlternateAllelesAs(vc2));
    }

    @Test
    public void testCreatingInsertionVariantContext() {
        List<Allele> alleles = Arrays.asList(delRef, ATC);
        VariantContext vc = insBuilder.alleles(alleles).make();

        Assert.assertEquals(vc.getChr(), insLoc);
        Assert.assertEquals(vc.getStart(), insLocStart);
        Assert.assertEquals(vc.getEnd(), insLocStop);
        Assert.assertEquals(vc.getType(), VariantContext.Type.INDEL);
        Assert.assertFalse(vc.isSNP());
        Assert.assertTrue(vc.isIndel());
        Assert.assertTrue(vc.isSimpleInsertion());
        Assert.assertFalse(vc.isSimpleDeletion());
        Assert.assertFalse(vc.isMixed());
        Assert.assertTrue(vc.isBiallelic());
        Assert.assertEquals(vc.getNAlleles(), 2);

        Assert.assertEquals(vc.getReference(), delRef);
        Assert.assertEquals(vc.getAlleles().size(), 2);
        Assert.assertEquals(vc.getAlternateAlleles().size(), 1);
        Assert.assertEquals(vc.getAlternateAllele(0), ATC);
        Assert.assertFalse(vc.hasGenotypes());

        Assert.assertEquals(vc.getSampleNames().size(), 0);
    }

    @Test
    public void testCreatingPartiallyCalledGenotype() {
        List<Allele> alleles = Arrays.asList(Aref, C);
        Genotype g = GenotypeBuilder.create("foo", Arrays.asList(C, Allele.NO_CALL));
        VariantContext vc = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, alleles).genotypes(g).make();

        Assert.assertTrue(vc.isSNP());
        Assert.assertEquals(vc.getNAlleles(), 2);
        Assert.assertTrue(vc.hasGenotypes());
        Assert.assertFalse(vc.isMonomorphicInSamples());
        Assert.assertTrue(vc.isPolymorphicInSamples());
        Assert.assertEquals(vc.getGenotype("foo"), g);
        Assert.assertEquals(vc.getCalledChrCount(), 1); // we only have 1 called chromosomes, we exclude the NO_CALL one isn't called
        Assert.assertEquals(vc.getCalledChrCount(Aref), 0);
        Assert.assertEquals(vc.getCalledChrCount(C), 1);
        Assert.assertFalse(vc.getGenotype("foo").isHet());
        Assert.assertFalse(vc.getGenotype("foo").isHom());
        Assert.assertFalse(vc.getGenotype("foo").isNoCall());
        Assert.assertFalse(vc.getGenotype("foo").isHom());
        Assert.assertTrue(vc.getGenotype("foo").isMixed());
        Assert.assertEquals(vc.getGenotype("foo").getType(), GenotypeType.MIXED);
    }

    @Test (expectedExceptions = Exception.class)
    public void testBadConstructorArgs1() {
        new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, Arrays.asList(delRef, ATCref)).make();
    }

    @Test (expectedExceptions = Exception.class)
    public void testBadConstructorArgs2() {
        new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, Arrays.asList(delRef, del)).make();
    }

    @Test (expectedExceptions = Exception.class)
    public void testBadConstructorArgs3() {
        new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, Arrays.asList(del)).make();
    }

    @Test (expectedExceptions = Throwable.class)
    public void testBadConstructorArgs4() {
        new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, Collections.<Allele>emptyList()).make();
    }

    @Test (expectedExceptions = Exception.class)
    public void testBadConstructorArgsDuplicateAlleles1() {
        new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, Arrays.asList(Aref, T, T)).make();
    }

    @Test (expectedExceptions = Exception.class)
    public void testBadConstructorArgsDuplicateAlleles2() {
        new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, Arrays.asList(Aref, A)).make();
    }

    @Test (expectedExceptions = Throwable.class)
    public void testBadLoc1() {
        List<Allele> alleles = Arrays.asList(Aref, T, del);
        new VariantContextBuilder("test", delLoc, delLocStart, delLocStop, alleles).make();
    }

    @Test (expectedExceptions = Throwable.class)
    public void testBadID1() {
        new VariantContextBuilder("test", delLoc, delLocStart, delLocStop, Arrays.asList(Aref, T)).id(null).make();
    }

    @Test (expectedExceptions = Exception.class)
    public void testBadID2() {
        new VariantContextBuilder("test", delLoc, delLocStart, delLocStop, Arrays.asList(Aref, T)).id("").make();
    }

    @Test (expectedExceptions = Throwable.class)
    public void testBadPError() {
        new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, Arrays.asList(delRef, ATCref)).log10PError(0.5).make();
    }

    @Test
    public void testAccessingSimpleSNPGenotypes() {
        List<Allele> alleles = Arrays.asList(Aref, T);

        Genotype g1 = GenotypeBuilder.create("AA", Arrays.asList(Aref, Aref));
        Genotype g2 = GenotypeBuilder.create("AT", Arrays.asList(Aref, T));
        Genotype g3 = GenotypeBuilder.create("TT", Arrays.asList(T, T));

        VariantContext vc = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, alleles)
                .genotypes(g1, g2, g3).make();

        Assert.assertTrue(vc.hasGenotypes());
        Assert.assertFalse(vc.isMonomorphicInSamples());
        Assert.assertTrue(vc.isPolymorphicInSamples());
        Assert.assertEquals(vc.getSampleNames().size(), 3);

        Assert.assertEquals(vc.getGenotypes().size(), 3);
        Assert.assertEquals(vc.getGenotypes().get("AA"), g1);
        Assert.assertEquals(vc.getGenotype("AA"), g1);
        Assert.assertEquals(vc.getGenotypes().get("AT"), g2);
        Assert.assertEquals(vc.getGenotype("AT"), g2);
        Assert.assertEquals(vc.getGenotypes().get("TT"), g3);
        Assert.assertEquals(vc.getGenotype("TT"), g3);

        Assert.assertTrue(vc.hasGenotype("AA"));
        Assert.assertTrue(vc.hasGenotype("AT"));
        Assert.assertTrue(vc.hasGenotype("TT"));
        Assert.assertFalse(vc.hasGenotype("foo"));
        Assert.assertFalse(vc.hasGenotype("TTT"));
        Assert.assertFalse(vc.hasGenotype("at"));
        Assert.assertFalse(vc.hasGenotype("tt"));

        Assert.assertEquals(vc.getCalledChrCount(), 6);
        Assert.assertEquals(vc.getCalledChrCount(Aref), 3);
        Assert.assertEquals(vc.getCalledChrCount(T), 3);
    }

    @Test
    public void testAccessingCompleteGenotypes() {
        List<Allele> alleles = Arrays.asList(Aref, T, ATC);

        Genotype g1 = GenotypeBuilder.create("AA", Arrays.asList(Aref, Aref));
        Genotype g2 = GenotypeBuilder.create("AT", Arrays.asList(Aref, T));
        Genotype g3 = GenotypeBuilder.create("TT", Arrays.asList(T, T));
        Genotype g4 = GenotypeBuilder.create("Td", Arrays.asList(T, ATC));
        Genotype g5 = GenotypeBuilder.create("dd", Arrays.asList(ATC, ATC));
        Genotype g6 = GenotypeBuilder.create("..", Arrays.asList(Allele.NO_CALL, Allele.NO_CALL));

        VariantContext vc = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, alleles)
                .genotypes(g1, g2, g3, g4, g5, g6).make();

        Assert.assertTrue(vc.hasGenotypes());
        Assert.assertFalse(vc.isMonomorphicInSamples());
        Assert.assertTrue(vc.isPolymorphicInSamples());
        Assert.assertEquals(vc.getGenotypes().size(), 6);

        Assert.assertEquals(3, vc.getGenotypes(Arrays.asList("AA", "Td", "dd")).size());

        Assert.assertEquals(10, vc.getCalledChrCount());
        Assert.assertEquals(3, vc.getCalledChrCount(Aref));
        Assert.assertEquals(4, vc.getCalledChrCount(T));
        Assert.assertEquals(3, vc.getCalledChrCount(ATC));
        Assert.assertEquals(2, vc.getCalledChrCount(Allele.NO_CALL));
    }

    @Test
    public void testAccessingRefGenotypes() {
        List<Allele> alleles1 = Arrays.asList(Aref, T);
        List<Allele> alleles2 = Arrays.asList(Aref);
        List<Allele> alleles3 = Arrays.asList(Aref, T);
        for ( List<Allele> alleles : Arrays.asList(alleles1, alleles2, alleles3)) {
            Genotype g1 = GenotypeBuilder.create("AA1", Arrays.asList(Aref, Aref));
            Genotype g2 = GenotypeBuilder.create("AA2", Arrays.asList(Aref, Aref));
            Genotype g3 = GenotypeBuilder.create("..", Arrays.asList(Allele.NO_CALL, Allele.NO_CALL));
            VariantContext vc = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, alleles)
                    .genotypes(g1, g2, g3).make();

            Assert.assertTrue(vc.hasGenotypes());
            Assert.assertTrue(vc.isMonomorphicInSamples());
            Assert.assertFalse(vc.isPolymorphicInSamples());
            Assert.assertEquals(vc.getGenotypes().size(), 3);

            Assert.assertEquals(4, vc.getCalledChrCount());
            Assert.assertEquals(4, vc.getCalledChrCount(Aref));
            Assert.assertEquals(0, vc.getCalledChrCount(T));
            Assert.assertEquals(2, vc.getCalledChrCount(Allele.NO_CALL));
        }
    }

    @Test
    public void testFilters() {
        List<Allele> alleles = Arrays.asList(Aref, T);
        Genotype g1 = GenotypeBuilder.create("AA", Arrays.asList(Aref, Aref));
        Genotype g2 = GenotypeBuilder.create("AT", Arrays.asList(Aref, T));

        VariantContext vc = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, alleles).genotypes(g1, g2).make();

        Assert.assertTrue(vc.isNotFiltered());
        Assert.assertFalse(vc.isFiltered());
        Assert.assertEquals(0, vc.getFilters().size());
        Assert.assertFalse(vc.filtersWereApplied());
        Assert.assertNull(vc.getFiltersMaybeNull());

        vc = new VariantContextBuilder(vc).filters("BAD_SNP_BAD!").make();

        Assert.assertFalse(vc.isNotFiltered());
        Assert.assertTrue(vc.isFiltered());
        Assert.assertEquals(1, vc.getFilters().size());
        Assert.assertTrue(vc.filtersWereApplied());
        Assert.assertNotNull(vc.getFiltersMaybeNull());

        Set<String> filters = new HashSet<String>(Arrays.asList("BAD_SNP_BAD!", "REALLY_BAD_SNP", "CHRIST_THIS_IS_TERRIBLE"));
        vc = new VariantContextBuilder(vc).filters(filters).make();

        Assert.assertFalse(vc.isNotFiltered());
        Assert.assertTrue(vc.isFiltered());
        Assert.assertEquals(3, vc.getFilters().size());
        Assert.assertTrue(vc.filtersWereApplied());
        Assert.assertNotNull(vc.getFiltersMaybeNull());
    }

    @Test
    public void testGetGenotypeCounts() {
        List<Allele> alleles = Arrays.asList(Aref, T);
        Genotype g1 = GenotypeBuilder.create("AA", Arrays.asList(Aref, Aref));
        Genotype g2 = GenotypeBuilder.create("AT", Arrays.asList(Aref, T));
        Genotype g3 = GenotypeBuilder.create("TT", Arrays.asList(T, T));
        Genotype g4 = GenotypeBuilder.create("A.", Arrays.asList(Aref, Allele.NO_CALL));
        Genotype g5 = GenotypeBuilder.create("..", Arrays.asList(Allele.NO_CALL, Allele.NO_CALL));

        // we need to create a new VariantContext each time
        VariantContext vc = new VariantContextBuilder("foo", snpLoc, snpLocStart, snpLocStop, alleles).genotypes(g1,g2,g3,g4,g5).make();
        Assert.assertEquals(1, vc.getHetCount());
        vc = new VariantContextBuilder("foo", snpLoc, snpLocStart, snpLocStop, alleles).genotypes(g1,g2,g3,g4,g5).make();
        Assert.assertEquals(1, vc.getHomRefCount());
        vc = new VariantContextBuilder("foo", snpLoc, snpLocStart, snpLocStop, alleles).genotypes(g1,g2,g3,g4,g5).make();
        Assert.assertEquals(1, vc.getHomVarCount());
        vc = new VariantContextBuilder("foo", snpLoc, snpLocStart, snpLocStop, alleles).genotypes(g1,g2,g3,g4,g5).make();
        Assert.assertEquals(1, vc.getMixedCount());
        vc = new VariantContextBuilder("foo", snpLoc, snpLocStart, snpLocStop, alleles).genotypes(g1,g2,g3,g4,g5).make();
        Assert.assertEquals(1, vc.getNoCallCount());
    }

        @Test
    public void testVCFfromGenotypes() {
        List<Allele> alleles = Arrays.asList(Aref, T);
        Genotype g1 = GenotypeBuilder.create("AA", Arrays.asList(Aref, Aref));
        Genotype g2 = GenotypeBuilder.create("AT", Arrays.asList(Aref, T));
        Genotype g3 = GenotypeBuilder.create("TT", Arrays.asList(T, T));
        Genotype g4 = GenotypeBuilder.create("..", Arrays.asList(Allele.NO_CALL, Allele.NO_CALL));
        VariantContext vc = new VariantContextBuilder("genotypes", snpLoc, snpLocStart, snpLocStop, alleles).genotypes(g1,g2,g3,g4).make();

        VariantContext vc12 = vc.subContextFromSamples(new HashSet<String>(Arrays.asList(g1.getSampleName(), g2.getSampleName())), true);
        VariantContext vc1 = vc.subContextFromSamples(new HashSet<String>(Arrays.asList(g1.getSampleName())), true);
        VariantContext vc23 = vc.subContextFromSamples(new HashSet<String>(Arrays.asList(g2.getSampleName(), g3.getSampleName())), true);
        VariantContext vc4 = vc.subContextFromSamples(new HashSet<String>(Arrays.asList(g4.getSampleName())), true);
        VariantContext vc14 = vc.subContextFromSamples(new HashSet<String>(Arrays.asList(g1.getSampleName(), g4.getSampleName())), true);

        Assert.assertTrue(vc12.isPolymorphicInSamples());
        Assert.assertTrue(vc23.isPolymorphicInSamples());
        Assert.assertTrue(vc1.isMonomorphicInSamples());
        Assert.assertTrue(vc4.isMonomorphicInSamples());
        Assert.assertTrue(vc14.isMonomorphicInSamples());

        Assert.assertTrue(vc12.isSNP());
        Assert.assertTrue(vc12.isVariant());
        Assert.assertTrue(vc12.isBiallelic());

        Assert.assertFalse(vc1.isSNP());
        Assert.assertFalse(vc1.isVariant());
        Assert.assertFalse(vc1.isBiallelic());

        Assert.assertTrue(vc23.isSNP());
        Assert.assertTrue(vc23.isVariant());
        Assert.assertTrue(vc23.isBiallelic());

        Assert.assertFalse(vc4.isSNP());
        Assert.assertFalse(vc4.isVariant());
        Assert.assertFalse(vc4.isBiallelic());

        Assert.assertFalse(vc14.isSNP());
        Assert.assertFalse(vc14.isVariant());
        Assert.assertFalse(vc14.isBiallelic());

        Assert.assertEquals(3, vc12.getCalledChrCount(Aref));
        Assert.assertEquals(1, vc23.getCalledChrCount(Aref));
        Assert.assertEquals(2, vc1.getCalledChrCount(Aref));
        Assert.assertEquals(0, vc4.getCalledChrCount(Aref));
        Assert.assertEquals(2, vc14.getCalledChrCount(Aref));
    }

    public void testGetGenotypeMethods() {
        Genotype g1 = GenotypeBuilder.create("AA", Arrays.asList(Aref, Aref));
        Genotype g2 = GenotypeBuilder.create("AT", Arrays.asList(Aref, T));
        Genotype g3 = GenotypeBuilder.create("TT", Arrays.asList(T, T));
        GenotypesContext gc = GenotypesContext.create(g1, g2, g3);
        VariantContext vc = new VariantContextBuilder("genotypes", snpLoc, snpLocStart, snpLocStop, Arrays.asList(Aref, T)).genotypes(gc).make();

        Assert.assertEquals(vc.getGenotype("AA"), g1);
        Assert.assertEquals(vc.getGenotype("AT"), g2);
        Assert.assertEquals(vc.getGenotype("TT"), g3);
        Assert.assertEquals(vc.getGenotype("CC"), null);

        Assert.assertEquals(vc.getGenotypes(), gc);
        Assert.assertEquals(vc.getGenotypes(Arrays.asList("AA", "AT")), Arrays.asList(g1, g2));
        Assert.assertEquals(vc.getGenotypes(Arrays.asList("AA", "TT")), Arrays.asList(g1, g3));
        Assert.assertEquals(vc.getGenotypes(Arrays.asList("AA", "AT", "TT")), Arrays.asList(g1, g2, g3));
        Assert.assertEquals(vc.getGenotypes(Arrays.asList("AA", "AT", "CC")), Arrays.asList(g1, g2));

        Assert.assertEquals(vc.getGenotype(0), g1);
        Assert.assertEquals(vc.getGenotype(1), g2);
        Assert.assertEquals(vc.getGenotype(2), g3);
    }

    // --------------------------------------------------------------------------------
    //
    // Test allele merging
    //
    // --------------------------------------------------------------------------------

    private class GetAllelesTest {
        List<Allele> alleles;
        String name;

        private GetAllelesTest(String name, Allele... arg) {
            this.name = name;
            this.alleles = Arrays.asList(arg);
        }

        public String toString() {
            return String.format("%s input=%s", name, alleles);
        }
    }

    @DataProvider(name = "getAlleles")
    public Object[][] mergeAllelesData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{new GetAllelesTest("A*",   Aref)});
        tests.add(new Object[]{new GetAllelesTest("A*/C", Aref, C)});
        tests.add(new Object[]{new GetAllelesTest("A*/C/T", Aref, C, T)});
        tests.add(new Object[]{new GetAllelesTest("A*/T/C", Aref, T, C)});
        tests.add(new Object[]{new GetAllelesTest("A*/C/T/ATC", Aref, C, T, ATC)});
        tests.add(new Object[]{new GetAllelesTest("A*/T/C/ATC", Aref, T, C, ATC)});
        tests.add(new Object[]{new GetAllelesTest("A*/ATC/T/C", Aref, ATC, T, C)});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "getAlleles")
    public void testMergeAlleles(GetAllelesTest cfg) {
        final List<Allele> altAlleles = cfg.alleles.subList(1, cfg.alleles.size());
        final VariantContext vc = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, cfg.alleles).make();

        Assert.assertEquals(vc.getAlleles(), cfg.alleles, "VC alleles not the same as input alleles");
        Assert.assertEquals(vc.getNAlleles(), cfg.alleles.size(), "VC getNAlleles not the same as input alleles size");
        Assert.assertEquals(vc.getAlternateAlleles(), altAlleles, "VC alt alleles not the same as input alt alleles");


        for ( int i = 0; i < cfg.alleles.size(); i++ ) {
            final Allele inputAllele = cfg.alleles.get(i);

            Assert.assertTrue(vc.hasAllele(inputAllele));
            if ( inputAllele.isReference() ) {
                final Allele nonRefVersion = Allele.create(inputAllele.getBases(), false);
                Assert.assertTrue(vc.hasAllele(nonRefVersion, true));
                Assert.assertFalse(vc.hasAllele(nonRefVersion, false));
            }

            Assert.assertEquals(inputAllele, vc.getAllele(inputAllele.getBaseString()));
            Assert.assertEquals(inputAllele, vc.getAllele(inputAllele.getBases()));

            if ( i > 0 ) { // it's an alt allele
                Assert.assertEquals(inputAllele, vc.getAlternateAllele(i-1));
            }
        }

        final Allele missingAllele = Allele.create("AACCGGTT"); // does not exist
        Assert.assertNull(vc.getAllele(missingAllele.getBases()));
        Assert.assertFalse(vc.hasAllele(missingAllele));
        Assert.assertFalse(vc.hasAllele(missingAllele, true));
    }

    private class SitesAndGenotypesVC {
        VariantContext vc, copy;
        String name;

        private SitesAndGenotypesVC(String name, VariantContext original) {
            this.name = name;
            this.vc = original;
            this.copy = new VariantContextBuilder(original).make();
        }

        public String toString() {
            return String.format("%s input=%s", name, vc);
        }
    }

    @DataProvider(name = "SitesAndGenotypesVC")
    public Object[][] MakeSitesAndGenotypesVCs() {
        Genotype g1 = GenotypeBuilder.create("AA", Arrays.asList(Aref, Aref));
        Genotype g2 = GenotypeBuilder.create("AT", Arrays.asList(Aref, T));
        Genotype g3 = GenotypeBuilder.create("TT", Arrays.asList(T, T));

        VariantContext sites = new VariantContextBuilder("sites", snpLoc, snpLocStart, snpLocStop, Arrays.asList(Aref, T)).make();
        VariantContext genotypes = new VariantContextBuilder(sites).source("genotypes").genotypes(g1, g2, g3).make();

        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{new SitesAndGenotypesVC("sites", sites)});
        tests.add(new Object[]{new SitesAndGenotypesVC("genotypes", genotypes)});

        return tests.toArray(new Object[][]{});
    }

    // --------------------------------------------------------------------------------
    //
    // Test modifying routines
    //
    // --------------------------------------------------------------------------------
    @Test(dataProvider = "SitesAndGenotypesVC")
    public void runModifyVCTests(SitesAndGenotypesVC cfg) {
        VariantContext modified = new VariantContextBuilder(cfg.vc).loc("chr2", 123, 123).make();
        Assert.assertEquals(modified.getChr(), "chr2");
        Assert.assertEquals(modified.getStart(), 123);
        Assert.assertEquals(modified.getEnd(), 123);

        modified = new VariantContextBuilder(cfg.vc).id("newID").make();
        Assert.assertEquals(modified.getID(), "newID");

        Set<String> newFilters = Collections.singleton("newFilter");
        modified = new VariantContextBuilder(cfg.vc).filters(newFilters).make();
        Assert.assertEquals(modified.getFilters(), newFilters);

        // test the behavior when the builder's attribute object is null
        modified = new VariantContextBuilder(modified).attributes(null).make();
        Assert.assertTrue(modified.getAttributes().isEmpty());
        modified = new VariantContextBuilder(modified).attributes(null).rmAttribute("AC").make();
        Assert.assertTrue(modified.getAttributes().isEmpty());
        modified = new VariantContextBuilder(modified).attributes(null).attribute("AC", 1).make();
        Assert.assertEquals(modified.getAttribute("AC"), 1);

        // test the behavior when the builder's attribute object is not initialized
        modified = new VariantContextBuilder(modified.getSource(), modified.getChr(), modified.getStart(), modified.getEnd(), modified.getAlleles()).attribute("AC", 1).make();

        // test normal attribute modification
        modified = new VariantContextBuilder(cfg.vc).attribute("AC", 1).make();
        Assert.assertEquals(modified.getAttribute("AC"), 1);
        modified = new VariantContextBuilder(modified).attribute("AC", 2).make();
        Assert.assertEquals(modified.getAttribute("AC"), 2);

        Genotype g1 = GenotypeBuilder.create("AA2", Arrays.asList(Aref, Aref));
        Genotype g2 = GenotypeBuilder.create("AT2", Arrays.asList(Aref, T));
        Genotype g3 = GenotypeBuilder.create("TT2", Arrays.asList(T, T));
        GenotypesContext gc = GenotypesContext.create(g1,g2,g3);
        modified = new VariantContextBuilder(cfg.vc).genotypes(gc).make();
        Assert.assertEquals(modified.getGenotypes(), gc);
        modified = new VariantContextBuilder(cfg.vc).noGenotypes().make();
        Assert.assertTrue(modified.getGenotypes().isEmpty());

        // test that original hasn't changed
        Assert.assertEquals(cfg.vc.getChr(), cfg.copy.getChr());
        Assert.assertEquals(cfg.vc.getStart(), cfg.copy.getStart());
        Assert.assertEquals(cfg.vc.getEnd(), cfg.copy.getEnd());
        Assert.assertEquals(cfg.vc.getAlleles(), cfg.copy.getAlleles());
        Assert.assertEquals(cfg.vc.getAttributes(), cfg.copy.getAttributes());
        Assert.assertEquals(cfg.vc.getID(), cfg.copy.getID());
        Assert.assertEquals(cfg.vc.getGenotypes(), cfg.copy.getGenotypes());
        Assert.assertEquals(cfg.vc.getLog10PError(), cfg.copy.getLog10PError());
        Assert.assertEquals(cfg.vc.getFilters(), cfg.copy.getFilters());
    }

    // --------------------------------------------------------------------------------
    //
    // Test subcontext
    //
    // --------------------------------------------------------------------------------
    private class SubContextTest {
        Set<String> samples;
        boolean updateAlleles;

        private SubContextTest(Collection<String> samples, boolean updateAlleles) {
            this.samples = new HashSet<String>(samples);
            this.updateAlleles = updateAlleles;
        }

        public String toString() {
            return String.format("%s samples=%s updateAlleles=%b", "SubContextTest", samples, updateAlleles);
        }
    }

    @DataProvider(name = "SubContextTest")
    public Object[][] MakeSubContextTest() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( boolean updateAlleles : Arrays.asList(true, false)) {
            tests.add(new Object[]{new SubContextTest(Collections.<String>emptySet(), updateAlleles)});
            tests.add(new Object[]{new SubContextTest(Collections.singleton("MISSING"), updateAlleles)});
            tests.add(new Object[]{new SubContextTest(Collections.singleton("AA"), updateAlleles)});
            tests.add(new Object[]{new SubContextTest(Collections.singleton("AT"), updateAlleles)});
            tests.add(new Object[]{new SubContextTest(Collections.singleton("TT"), updateAlleles)});
            tests.add(new Object[]{new SubContextTest(Arrays.asList("AA", "AT"), updateAlleles)});
            tests.add(new Object[]{new SubContextTest(Arrays.asList("AA", "AT", "TT"), updateAlleles)});
            tests.add(new Object[]{new SubContextTest(Arrays.asList("AA", "AT", "MISSING"), updateAlleles)});
            tests.add(new Object[]{new SubContextTest(Arrays.asList("AA", "AT", "TT", "MISSING"), updateAlleles)});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "SubContextTest")
    public void runSubContextTest(SubContextTest cfg) {
        Genotype g1 = GenotypeBuilder.create("AA", Arrays.asList(Aref, Aref));
        Genotype g2 = GenotypeBuilder.create("AT", Arrays.asList(Aref, T));
        Genotype g3 = GenotypeBuilder.create("TT", Arrays.asList(T, T));

        GenotypesContext gc = GenotypesContext.create(g1, g2, g3);
        VariantContext vc = new VariantContextBuilder("genotypes", snpLoc, snpLocStart, snpLocStop, Arrays.asList(Aref, T)).genotypes(gc).make();
        VariantContext sub = vc.subContextFromSamples(cfg.samples, cfg.updateAlleles);

        // unchanged attributes should be the same
        Assert.assertEquals(sub.getChr(), vc.getChr());
        Assert.assertEquals(sub.getStart(), vc.getStart());
        Assert.assertEquals(sub.getEnd(), vc.getEnd());
        Assert.assertEquals(sub.getLog10PError(), vc.getLog10PError());
        Assert.assertEquals(sub.getFilters(), vc.getFilters());
        Assert.assertEquals(sub.getID(), vc.getID());
        Assert.assertEquals(sub.getAttributes(), vc.getAttributes());

        Set<Genotype> expectedGenotypes = new HashSet<Genotype>();
        if ( cfg.samples.contains(g1.getSampleName()) ) expectedGenotypes.add(g1);
        if ( cfg.samples.contains(g2.getSampleName()) ) expectedGenotypes.add(g2);
        if ( cfg.samples.contains(g3.getSampleName()) ) expectedGenotypes.add(g3);
        GenotypesContext expectedGC = GenotypesContext.copy(expectedGenotypes);

        // these values depend on the results of sub
        if ( cfg.updateAlleles ) {
            // do the work to see what alleles should be here, and which not
            Set<Allele> alleles = new HashSet<Allele>();
            for ( final Genotype g : expectedGC ) alleles.addAll(g.getAlleles());
            if ( ! alleles.contains(Aref) ) alleles.add(Aref); // always have the reference
            Assert.assertEquals(new HashSet<Allele>(sub.getAlleles()), alleles);
        } else {
            // not updating alleles -- should be the same
            Assert.assertEquals(sub.getAlleles(), vc.getAlleles());
        }

        // same sample names => success
        Assert.assertEquals(sub.getGenotypes().getSampleNames(), expectedGC.getSampleNames());
    }

    // --------------------------------------------------------------------------------
    //
    // Test sample name functions
    //
    // --------------------------------------------------------------------------------
    private class SampleNamesTest {
        List<String> sampleNames;
        List<String> sampleNamesInOrder;

        private SampleNamesTest(List<String> sampleNames, List<String> sampleNamesInOrder) {
            this.sampleNamesInOrder = sampleNamesInOrder;
            this.sampleNames = sampleNames;
        }

        public String toString() {
            return String.format("%s samples=%s order=%s", "SampleNamesTest", sampleNames, sampleNamesInOrder);
        }
    }

    @DataProvider(name = "SampleNamesTest")
    public Object[][] MakeSampleNamesTest() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{new SampleNamesTest(Arrays.asList("1"), Arrays.asList("1"))});
        tests.add(new Object[]{new SampleNamesTest(Arrays.asList("2", "1"), Arrays.asList("1", "2"))});
        tests.add(new Object[]{new SampleNamesTest(Arrays.asList("1", "2"), Arrays.asList("1", "2"))});
        tests.add(new Object[]{new SampleNamesTest(Arrays.asList("1", "2", "3"), Arrays.asList("1", "2", "3"))});
        tests.add(new Object[]{new SampleNamesTest(Arrays.asList("2", "1", "3"), Arrays.asList("1", "2", "3"))});
        tests.add(new Object[]{new SampleNamesTest(Arrays.asList("2", "3", "1"), Arrays.asList("1", "2", "3"))});
        tests.add(new Object[]{new SampleNamesTest(Arrays.asList("3", "1", "2"), Arrays.asList("1", "2", "3"))});
        tests.add(new Object[]{new SampleNamesTest(Arrays.asList("3", "2", "1"), Arrays.asList("1", "2", "3"))});
        tests.add(new Object[]{new SampleNamesTest(Arrays.asList("NA2", "NA1"), Arrays.asList("NA1", "NA2"))});

        return tests.toArray(new Object[][]{});
    }

    private final static void assertGenotypesAreInOrder(Iterable<Genotype> gIt, List<String> names) {
        int i = 0;
        for ( final Genotype g : gIt ) {
            Assert.assertEquals(g.getSampleName(), names.get(i), "Unexpected genotype ordering");
            i++;
        }
    }


    @Test(dataProvider = "SampleNamesTest")
    public void runSampleNamesTest(SampleNamesTest cfg) {
        GenotypesContext gc = GenotypesContext.create(cfg.sampleNames.size());
        for ( final String name : cfg.sampleNames ) {
            gc.add(GenotypeBuilder.create(name, Arrays.asList(Aref, T)));
        }

        VariantContext vc = new VariantContextBuilder("genotypes", snpLoc, snpLocStart, snpLocStop, Arrays.asList(Aref, T)).genotypes(gc).make();

        // same sample names => success
        Assert.assertTrue(vc.getSampleNames().equals(new HashSet<String>(cfg.sampleNames)), "vc.getSampleNames() = " + vc.getSampleNames());
        Assert.assertEquals(vc.getSampleNamesOrderedByName(), cfg.sampleNamesInOrder, "vc.getSampleNamesOrderedByName() = " + vc.getSampleNamesOrderedByName());

        assertGenotypesAreInOrder(vc.getGenotypesOrderedByName(), cfg.sampleNamesInOrder);
        assertGenotypesAreInOrder(vc.getGenotypesOrderedBy(cfg.sampleNames), cfg.sampleNames);
    }

    @Test
    public void testGenotypeCounting() {
        Genotype noCall = GenotypeBuilder.create("nocall", Arrays.asList(Allele.NO_CALL));
        Genotype mixed  = GenotypeBuilder.create("mixed", Arrays.asList(Aref, Allele.NO_CALL));
        Genotype homRef = GenotypeBuilder.create("homRef", Arrays.asList(Aref, Aref));
        Genotype het    = GenotypeBuilder.create("het", Arrays.asList(Aref, T));
        Genotype homVar = GenotypeBuilder.create("homVar", Arrays.asList(T, T));

        List<Genotype> allGenotypes = Arrays.asList(noCall, mixed, homRef, het, homVar);
        final int nCycles = allGenotypes.size() * 10;

        for ( int i = 0; i < nCycles; i++ ) {
            int nNoCall = 0, nNoCallAlleles = 0, nA = 0, nT = 0, nMixed = 0, nHomRef = 0, nHet = 0, nHomVar = 0;
            int nSamples = 0;
            GenotypesContext gc = GenotypesContext.create();
            for ( int j = 0; j < i; j++ ) {
                nSamples++;
                Genotype g = allGenotypes.get(j % allGenotypes.size());
                final String name = String.format("%s_%d%d", g.getSampleName(), i, j);
                gc.add(GenotypeBuilder.create(name, g.getAlleles()));
                switch ( g.getType() ) {
                    case NO_CALL: nNoCall++; nNoCallAlleles++; break;
                    case HOM_REF: nA += 2; nHomRef++; break;
                    case HET: nA++; nT++; nHet++; break;
                    case HOM_VAR: nT += 2; nHomVar++; break;
                    case MIXED: nA++; nNoCallAlleles++; nMixed++; break;
                    default: throw new RuntimeException("Unexpected genotype type " + g.getType());
                }

            }

            VariantContext vc = new VariantContextBuilder("genotypes", snpLoc, snpLocStart, snpLocStop, Arrays.asList(Aref, T)).genotypes(gc).make();
            Assert.assertEquals(vc.getNSamples(), nSamples);
            if ( nSamples > 0 ) {
                Assert.assertEquals(vc.isPolymorphicInSamples(), nT > 0);
                Assert.assertEquals(vc.isMonomorphicInSamples(), nT == 0);
            }
            Assert.assertEquals(vc.getCalledChrCount(), nA + nT);

            Assert.assertEquals(vc.getCalledChrCount(Allele.NO_CALL), nNoCallAlleles);
            Assert.assertEquals(vc.getCalledChrCount(Aref), nA);
            Assert.assertEquals(vc.getCalledChrCount(T), nT);

            Assert.assertEquals(vc.getNoCallCount(), nNoCall);
            Assert.assertEquals(vc.getHomRefCount(), nHomRef);
            Assert.assertEquals(vc.getHetCount(), nHet);
            Assert.assertEquals(vc.getHomVarCount(), nHomVar);
            Assert.assertEquals(vc.getMixedCount(), nMixed);
        }
    }
}
