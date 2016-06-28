package picard.fingerprint;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.QualityUtil;
import htsjdk.variant.variantcontext.Allele;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.util.MathUtil;

import java.util.Collections;
import java.util.List;

import static picard.util.TestNGUtil.assertEqualDoubleArrays;
import static java.lang.Math.log10;

/**
 * Basic tests for HaplotypeProbabilities and derived classes
 *
 * @author yossi farjoun
 */
public class HaplotypeProbabilitiesTest {

    static Snp snp1, snp2;
    static HaplotypeBlock hb1, hb2;

    @BeforeTest
    public static void initializeHaplotypeBlock() {
        snp1 = new Snp("SNP1", "test", 1, (byte) 'A', (byte) 'T', 0.25, Collections.<String>emptyList());
        snp2 = new Snp("SNP2", "test", 2, (byte) 'A', (byte) 'G', 0.5, Collections.<String>emptyList());

        hb1 = new HaplotypeBlock(.25);
        hb1.addSnp(snp1);

        hb2 = new HaplotypeBlock(.4);
        hb2.addSnp(snp1);
        hb2.addSnp(snp2);
    }

    @DataProvider(name = "dataTestpEvidenceGivenPriorFromGLs")
    public Object[][] dataTestpEvidenceGivenPriorFromGLs() {
        return new Object[][]{
                new Object[]{
                        new HaplotypeProbabilitiesFromGenotypeLikelihoods(hb1),
                        Collections.singletonList(snp1),
                        Collections.singletonList(false),
                        Collections.singletonList(new double[]{0, -1, -2})},

                new Object[]{
                        new HaplotypeProbabilitiesFromGenotypeLikelihoods(hb1),
                        Collections.singletonList(snp1),
                        Collections.singletonList(true),
                        Collections.singletonList(new double[]{0, -1, -2})},

                new Object[]{
                        new HaplotypeProbabilitiesFromGenotypeLikelihoods(hb2),
                        CollectionUtil.makeList(snp1, snp2),
                        CollectionUtil.makeList(false, false),
                        CollectionUtil.makeList(new double[]{0, -1, -2}, new double[]{0, -1, -2})},

                new Object[]{
                        new HaplotypeProbabilitiesFromGenotypeLikelihoods(hb2),
                        CollectionUtil.makeList(snp1, snp2),
                        CollectionUtil.makeList(false, false),
                        CollectionUtil.makeList(new double[]{0, -1, -2}, new double[]{-1, 0, -1})},

                new Object[]{
                        new HaplotypeProbabilitiesFromGenotypeLikelihoods(hb2),
                        CollectionUtil.makeList(snp1, snp2),
                        CollectionUtil.makeList(false, false),
                        CollectionUtil.makeList(new double[]{0, -1, -2}, new double[]{-2, -1, 0})},
        };
    }

    @Test(dataProvider = "dataTestpEvidenceGivenPriorFromGLs")
    public void testpEvidenceGivenPriorFromGLs(final HaplotypeProbabilitiesFromGenotypeLikelihoods hp, final List<Snp> snps, final List<Boolean> swaps, final List<double[]> GLs) throws Exception {

        for (int i = 0; i < snps.size(); ++i) {
            final Allele a = Allele.create(swaps.get(i) ? snps.get(i).getAllele2() : snps.get(i).getAllele1());
            final Allele b = Allele.create(swaps.get(i) ? snps.get(i).getAllele1() : snps.get(i).getAllele2());

            hp.addToLogLikelihoods(snps.get(i), CollectionUtil.makeList(a, b), GLs.get(i));
        }

        final double[] logLikelihood = new double[3];
        for (int genotype = 0; genotype < 3; genotype++) {
            logLikelihood[genotype] = log10(hp.getHaplotype().getHaplotypeFrequency(genotype));
            for (int i = 0; i < GLs.size(); i++) {
                final double[] genotypeLogLikelihoods = GLs.get(i);
                if (swaps.get(i))
                    logLikelihood[genotype] += genotypeLogLikelihoods[2 - genotype];
                else
                    logLikelihood[genotype] += genotypeLogLikelihoods[genotype];
            }

        }
        assertEqualDoubleArrays(hp.getPosteriorProbabilities(), MathUtil.pNormalizeLogProbability(logLikelihood), 1e-10);
    }

    @DataProvider(name = "dataTestHaplotypeProbabilitiesFromSequenceAddToProbs")
    public Object[][] dataTestHaplotypeProbabilitiesFromSequenceAddToProbs() {
        return new Object[][]{
                {new HaplotypeProbabilitiesFromSequence(hb1), snp1, new byte[]{}, 7},
                {new HaplotypeProbabilitiesFromSequence(hb1), snp1, new byte[]{'A'}, 7},
                {new HaplotypeProbabilitiesFromSequence(hb1), snp1, new byte[]{'G'}, 7},
                {new HaplotypeProbabilitiesFromSequence(hb1), snp1, new byte[]{'T'}, 7},
                {new HaplotypeProbabilitiesFromSequence(hb1), snp1, new byte[]{'A', 'T', 'A', 'A', 'A', 'A', 'A', 'A'}, 7},
                {new HaplotypeProbabilitiesFromSequence(hb1), snp1, new byte[]{'A', 'T', 'A', 'A', 'G', 'A', 'A', 'A'}, 7},
                {new HaplotypeProbabilitiesFromSequence(hb1), snp1, new byte[]{'T', 'T', 'A', 'A', 'A', 'T', 'T', 'A', 'T'}, 7},
                {new HaplotypeProbabilitiesFromSequence(hb1), snp1, new byte[]{'T', 'A', 'T', 'T', 'T', 'T', 'T', 'T', 'A'}, 7}
        };
    }

    @Test(dataProvider = "dataTestHaplotypeProbabilitiesFromSequenceAddToProbs")
    public void testHaplotypeProbabilitiesFromSequenceAddToProbs(final HaplotypeProbabilitiesFromSequence hp, final Snp snp, final byte[] bases, final int qual) throws Exception {

        for (final byte base : bases) {
            hp.addToProbs(snp, base, (byte) qual);
        }

        final double pError = QualityUtil.getErrorProbabilityFromPhredScore(qual);
        final double[] logLikelihood = new double[3];

        for (int genotype = 0; genotype < 3; genotype++) {
            logLikelihood[genotype] = log10(hp.getHaplotype().getHaplotypeFrequency(genotype));
            for (final byte a : bases) {
                final double theta = 0.5 * genotype;
                if (a == snp.getAllele1())
                    logLikelihood[genotype] += log10((1 - theta) * (1 - pError) + theta * pError);
                if (a == snp.getAllele2())
                    logLikelihood[genotype] += log10((1 - theta) * (pError) + theta * (1 - pError));
            }
        }
        final double[] posterior = MathUtil.pNormalizeLogProbability(logLikelihood);
        assertEqualDoubleArrays(hp.getPosteriorProbabilities(), posterior, 1e-10);
    }

    @DataProvider(name = "dataTestHaplotypeProbabilitiesFromContaminatorSequenceAddToProbs")
    public Object[][] dataTestHaplotypeProbabilitiesFromContaminatorSequenceAddToProbs() {
        return new Object[][]{
                {new HaplotypeProbabilitiesFromContaminatorSequence(hb1, .1), snp1, 0, 0},
                {new HaplotypeProbabilitiesFromContaminatorSequence(hb1, .1), snp1, 0, 1},
                {new HaplotypeProbabilitiesFromContaminatorSequence(hb1, .1), snp1, 0, 76},
                {new HaplotypeProbabilitiesFromContaminatorSequence(hb1, .1), snp1, 3, 76},
                {new HaplotypeProbabilitiesFromContaminatorSequence(hb1, .1), snp1, 7, 76},
                {new HaplotypeProbabilitiesFromContaminatorSequence(hb1, .1), snp1, 35, 76},
                {new HaplotypeProbabilitiesFromContaminatorSequence(hb1, .1), snp1, 40, 76},
                {new HaplotypeProbabilitiesFromContaminatorSequence(hb1, .1), snp1, 45, 76},
                {new HaplotypeProbabilitiesFromContaminatorSequence(hb1, .1), snp1, 69, 76},
                {new HaplotypeProbabilitiesFromContaminatorSequence(hb1, .1), snp1, 73, 76},
                {new HaplotypeProbabilitiesFromContaminatorSequence(hb1, .1), snp1, 76, 76}
        };
    }

    static final int[] genotypes = {0, 1, 2};

    @Test(dataProvider = "dataTestHaplotypeProbabilitiesFromContaminatorSequenceAddToProbs")
    public void testHaplotypeProbabilitiesFromContaminatorSequenceAddToProbs(final HaplotypeProbabilitiesFromContaminatorSequence hp, final Snp snp, final int nAlt, final int nTotal) throws Exception {

        final byte qual = 7;
        for (int i = 0; i < nAlt; i++) {
            hp.addToProbs(snp, snp.getAllele2(), qual);
        }
        for (int i = nAlt; i < nTotal; i++) {
            hp.addToProbs(snp, snp.getAllele1(), qual);
        }

        final double pError = QualityUtil.getErrorProbabilityFromPhredScore(qual);
        final double[] unnormalizedLikelihood = {0d, 0d, 0d};

        for (final int contG : genotypes) {
            for (final int mainG : genotypes) {
                final double pAlt = (hp.contamination * contG + (1 - hp.contamination) * mainG) / 2;
                double l = hp.getHaplotype().getHaplotypeFrequency(mainG);
                for (int i = 0; i < nAlt; i++) {
                    l *= pAlt * (1 - pError) + (1 - pAlt) * pError;
                }
                for (int i = nAlt; i < nTotal; i++) {
                    l *= pAlt * (pError) + (1 - pAlt) * (1 - pError);
                }
                unnormalizedLikelihood[contG] += l;
            }
        }
        final double[] likelihood = MathUtil.pNormalizeVector(unnormalizedLikelihood);

        assertEqualDoubleArrays(hp.getLikelihoods(), likelihood, 1e-10);
    }
}
