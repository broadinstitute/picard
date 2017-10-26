package picard.util;

import htsjdk.variant.variantcontext.*;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.vcf.VcfTestUtils;

import java.util.*;

public class AlleleSubsettingUtilsTest {

    private static final Allele Aref = Allele.create("A", true);
    private static final Allele C = Allele.create("C");
    private static final Allele G = Allele.create("G");

    final List<Allele> AA = Arrays.asList(Aref, Aref);
    final List<Allele> AC = Arrays.asList(Aref, C);
    final List<Allele> CC = Arrays.asList(C, C);
    final List<Allele> AG = Arrays.asList(Aref, G);
    final List<Allele> GG = Arrays.asList(G, G);
    final List<Allele> ACG = Arrays.asList(Aref, C, G);

    @Test(dataProvider = "updatePLsAndADData")
    public void testUpdatePLsAndADData(final VariantContext originalVC,
                                       final VariantContext selectedVC,
                                       final List<Genotype> expectedGenotypes) {
        // initialize cache of allele anyploid indices
        for (final Genotype genotype : originalVC.getGenotypes()) {
            GenotypeLikelihoods.initializeAnyploidPLIndexToAlleleIndices(originalVC.getNAlleles() - 1, genotype.getPloidy());
        }

        final VariantContext selectedVCwithGTs = new VariantContextBuilder(selectedVC).genotypes(originalVC.getGenotypes()).make();

        final GenotypesContext oldGs = selectedVCwithGTs.getGenotypes();
        final GenotypesContext actual = selectedVCwithGTs.getNAlleles() == originalVC.getNAlleles() ? oldGs :
                AlleleSubsettingUtils.subsetAlleles(oldGs, originalVC.getAlleles(),
                        selectedVCwithGTs.getAlleles());

        Assert.assertEquals(actual.size(), expectedGenotypes.size());
        for (final Genotype expected : expectedGenotypes) {
            final Genotype actualGT = actual.get(expected.getSampleName());
            Assert.assertNotNull(actualGT);

            Assert.assertTrue(actualGT.sameGenotype(expected), String.format("expected same genotypes, found: %s and expected: %s", actual.toString(), expected.toString()));
        }
    }

    @DataProvider(name = "updatePLsAndADData")
    public Object[][] makeUpdatePLsAndADData() {
        List<Object[]> tests = new ArrayList<>();


        int i = 0;
        final VariantContext vcBase = new VariantContextBuilder("test", "20", 10, 10, AC).make();

        final double[] homRefPL = new double[]{0, 10, 100};
        final double[] hetPL = new double[]{10, 0, 100};
        final double[] homVarPL = new double[]{100, 10, 0};
        final double[] uninformative = new double[]{0, 0, 0};

        final Genotype base = new GenotypeBuilder("NA12878").DP(10).GQ(100).make();

        // the simple case where no selection occurs
        final Genotype aaGT = new GenotypeBuilder(base).alleles(AA).AD(new int[]{10, 2}).PL(homRefPL).GQ(8).make();
        final Genotype acGT = new GenotypeBuilder(base).alleles(AC).AD(new int[]{10, 2}).PL(hetPL).GQ(8).make();
        final Genotype ccGT = new GenotypeBuilder(base).alleles(CC).AD(new int[]{10, 2}).PL(homVarPL).GQ(8).make();

        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(aaGT).make(), new VariantContextBuilder(vcBase).alleles(AC).make(), Collections.singletonList(new GenotypeBuilder(aaGT).make())});
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(acGT).make(), new VariantContextBuilder(vcBase).alleles(AC).make(), Collections.singletonList(new GenotypeBuilder(acGT).make())});
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(ccGT).make(), new VariantContextBuilder(vcBase).alleles(AC).make(), Collections.singletonList(new GenotypeBuilder(ccGT).make())});

        // uninformative test cases
        final Genotype uninformativeGT = new GenotypeBuilder(base).alleles(CC).noAD().PL(uninformative).GQ(0).make();
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(uninformativeGT).make(), new VariantContextBuilder(vcBase).alleles(AC).make(), Collections.singletonList(uninformativeGT)});
        final Genotype emptyGT = new GenotypeBuilder(base).alleles(Collections.nCopies(2, Allele.NO_CALL)).noAD().noPL().noGQ().make();
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(emptyGT).make(), new VariantContextBuilder(vcBase).alleles(AC).make(), Collections.singletonList(emptyGT)});

        // actually subsetting down from multiple alt values
        final double[] homRef3AllelesPL = new double[]{0, -10, -20, -30, -40, -50};
        final double[] hetRefC3AllelesPL = new double[]{-10, 0, -20, -30, -40, -50};
        final double[] homC3AllelesPL = new double[]{-20, -10, 0, -30, -40, -50};
        final double[] hetRefG3AllelesPL = new double[]{-20, -10, -30, 0, -40, -50};
        final double[] hetCG3AllelesPL = new double[]{-20, -10, -30, -40, 0, -50}; // AA, AC, CC, AG, CG, GG
        final double[] homG3AllelesPL = new double[]{-20, -10, -30, -40, -50, 0};  // AA, AC, CC, AG, CG, GG

        final int[] homRef3AllelesAD = new int[]{20, 0, 1};
        final int[] hetRefC3AllelesAD = new int[]{10, 10, 1};
        final int[] homC3AllelesAD = new int[]{0, 20, 1};
        final int[] hetRefG3AllelesAD = new int[]{10, 0, 11};
        final int[] hetCG3AllelesAD = new int[]{0, 12, 11}; // AA, AC, CC, AG, CG, GG
        final int[] homG3AllelesAD = new int[]{0, 1, 21};  // AA, AC, CC, AG, CG, GG

        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).source("test-" + i++).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).AD(homRef3AllelesAD).PL(homRef3AllelesPL).make()).make(),
                new VariantContextBuilder(vcBase).source("test-" + i).alleles(AC).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(AA).PL(new double[]{0, -10, -20}).AD(new int[]{20, 0}).GQ(100).make())});

        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).source("test-" + i++).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).AD(hetRefC3AllelesAD).PL(hetRefC3AllelesPL).make()).make(),
                new VariantContextBuilder(vcBase).source("test-" + i).alleles(AC).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(AC).PL(new double[]{-10, 0, -20}).AD(new int[]{10, 10}).GQ(100).make())});

        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).source("test-" + i++).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).AD(homC3AllelesAD).PL(homC3AllelesPL).make()).make(),
                new VariantContextBuilder(vcBase).source("test-" + i).alleles(AC).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(CC).PL(new double[]{-20, -10, 0}).AD(new int[]{0, 20}).GQ(100).make())});

        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).source("test-" + i++).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).AD(hetRefG3AllelesAD).PL(hetRefG3AllelesPL).make()).make(),
                new VariantContextBuilder(vcBase).source("test-" + i).alleles(AG).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(AG).PL(new double[]{-20, 0, -50}).AD(new int[]{10, 11}).GQ(200).make())});

        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).source("test-" + i++).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).AD(hetCG3AllelesAD).PL(hetCG3AllelesPL).make()).make(),
                new VariantContextBuilder(vcBase).source("test-" + i).alleles(AG).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(AA).PL(new double[]{0, -20, -30}).AD(new int[]{0, 11}).GQ(200).make())});

        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).source("test-" + i++).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).AD(homG3AllelesAD).PL(homG3AllelesPL).make()).make(),
                new VariantContextBuilder(vcBase).source("test-" + i).alleles(AG).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(GG).PL(new double[]{-20, -40, 0}).AD(new int[]{0, 21}).GQ(200).make())});

        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).source("test-" + i++).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).AD(homG3AllelesAD).noPL().make()).make(),
                new VariantContextBuilder(vcBase).source("test-" + i).alleles(AG).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(AlleleSubsettingUtils.DIPLOID_NO_CALL).noPL().AD(new int[]{0, 21}).noGQ().make())});

        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name="swapAlleleData")
    Iterator<Object[]> swapAlleleData(){
        List<Object[]> tests = new ArrayList<>();

        final VariantContext vcBase = new VariantContextBuilder("test", "20", 10, 10, AC).make();
        final Genotype base = new GenotypeBuilder("NA12878").DP(10).GQ(100).make();

        tests.add(new Object[]{vcBase, C, G, new VariantContextBuilder(vcBase).alleles(AG).make(), true});
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(new GenotypeBuilder(base).alleles(AC).make()).make(), C, G,
                new VariantContextBuilder(vcBase).alleles(AG).genotypes(new GenotypeBuilder(base).alleles(AG).make()).make(), true});
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(new GenotypeBuilder(base).alleles(AA).make()).make(), C, G,
                new VariantContextBuilder(vcBase).alleles(AG).genotypes(new GenotypeBuilder(base).alleles(AA).make()).make(), true});
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(new GenotypeBuilder(base).alleles(CC).make()).make(), C, G,
                new VariantContextBuilder(vcBase).alleles(AG).genotypes(new GenotypeBuilder(base).alleles(GG).make()).make(), true});

        tests.add(new Object[]{vcBase, AlleleSubsettingUtils.NON_REF_ALLELE, G, null, false});


        return tests.iterator();
    }

    @Test(dataProvider = "swapAlleleData")
    void testSwapAlleles(final VariantContext ctx, final Allele oldAllele, final Allele newAllele, final VariantContext result, final boolean shouldSucceed){
        try {
            final VariantContext newVc = AlleleSubsettingUtils.swapAlleles(ctx, oldAllele, newAllele);
            VcfTestUtils.assertEquals(newVc, result);
        } catch (Exception e){
            if (shouldSucceed){
                throw e;
            } else return;
        }

    }
}
