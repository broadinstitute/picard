package picard.fingerprint;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.util.TestNGUtil;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * Created by farjoun on 5/29/15.
 */
public class HaplotypeProbabilityOfNormalGivenTumorTest {

    private final double maf = 0.4;
    private final Snp snp = new Snp("test", "chr1", 1, (byte) 'A', (byte) 'C', maf, Collections.singletonList("dummy"));
    private final HaplotypeBlock hb = new HaplotypeBlock(maf);

    @DataProvider(name = "testGetLikelihoodsData")
    public Iterator<Object[]> testGetLikelihoodsData() {
        final List<Object[]> testData = new ArrayList<>();

        //make sure that giving 0 pLoH doesn't change the underlying likelihoods:
        testData.add(new Object[]{0.0, new double[]{1, 0, 0}, new double[]{1, 0, 0}});
        testData.add(new Object[]{0.0, new double[]{0, 1, 0}, new double[]{0, 1, 0}});
        testData.add(new Object[]{0.0, new double[]{0, 0, 1}, new double[]{0, 0, 1}});
        testData.add(new Object[]{0.0, new double[]{0, 0.4, 0.6}, new double[]{0, 0.4, 0.6}});
        testData.add(new Object[]{0.0, new double[]{0.3, 0.7, 0}, new double[]{0.3, 0.7, 0}});

        final double pLoh = 0.1;
        //see that non zero pLoH changes the likelihood of a HET site as expected:
        testData.add(new Object[]{pLoh, new double[]{0, 1, 0}, new double[]{0, 1 - pLoh, 0}});

        //HOMs will change a little
        testData.add(new Object[]{pLoh, new double[]{1, 0, 0}, new double[]{1, pLoh/2, 0}});
        testData.add(new Object[]{pLoh, new double[]{0, 0, 1}, new double[]{0, pLoh/2, 1}});
        testData.add(new Object[]{pLoh, new double[]{.3, 0, .7}, new double[]{.3, pLoh/2, .7}});

        // check that the calculation is linear
        testData.add(new Object[]{pLoh, new double[]{0, 0.5, 0.5}, new double[]{0, 0.5 * (1 - pLoh/2), 0.5}});
        testData.add(new Object[]{pLoh, new double[]{0.5, 0.5, 0}, new double[]{0.5, 0.5 * (1 - pLoh/2), 0}});

        return testData.iterator();
    }

    @Test(dataProvider = "testGetLikelihoodsData")
    public void testGetLikelihoods(final double pLoH, final double[] tumorLikelihood, final double[] normalLikelihood) throws Exception {
        final HaplotypeProbabilities hp = new HaplotypeProbabilitiesFromGenotype(snp, hb, tumorLikelihood[0], tumorLikelihood[1], tumorLikelihood[2]);
        final HaplotypeProbabilities hpTumor = new HaplotypeProbabilityOfNormalGivenTumor(hp, pLoH);

        TestNGUtil.assertEqualDoubleArrays(hpTumor.getLikelihoods(), normalLikelihood, 0.0001);
    }
}