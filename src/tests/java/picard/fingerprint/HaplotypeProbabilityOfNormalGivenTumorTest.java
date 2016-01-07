package picard.fingerprint;

import picard.util.TestNGUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * Created by farjoun on 5/29/15.
 */
public class HaplotypeProbabilityOfNormalGivenTumorTest {

    private double maf = 0.4;
    private Snp snp = new Snp("test", "chr1", 1, (byte) 'A', (byte) 'C', maf, Collections.singletonList("dummy"));
    private HaplotypeBlock hb = new HaplotypeBlock(maf);

    @DataProvider(name = "testGetLikelihoodsData")
    public Iterator<Object[]> testGetLikelihoodsData() {
        List<Object[]> testData = new ArrayList<>();

        //make sure that giving 0 pLoH doesn't change the underlying likelihoods:
        testData.add(new Object[]{0.0, new double[]{1, 0, 0}, new double[]{1, 0, 0}});
        testData.add(new Object[]{0.0, new double[]{0, 1, 0}, new double[]{0, 1, 0}});
        testData.add(new Object[]{0.0, new double[]{0, 0, 1}, new double[]{0, 0, 1}});
        testData.add(new Object[]{0.0, new double[]{0, 0.4, 0.6}, new double[]{0, 0.4, 0.6}});
        testData.add(new Object[]{0.0, new double[]{0.3, 0.7, 0}, new double[]{0.3, 0.7, 0}});

        //make sure that pLoH will not affect HOM likelihoods:
        testData.add(new Object[]{0.1, new double[]{1, 0, 0}, new double[]{1, 0, 0}});
        testData.add(new Object[]{0.2, new double[]{0, 0, 1}, new double[]{0, 0, 1}});
        testData.add(new Object[]{0.3, new double[]{.3, 0, .7}, new double[]{.3, 0, .7}});

        //see that non zero pLoH changes the likelihood of a HET site as expected:
        testData.add(new Object[]{0.1, new double[]{0, 1, 0}, new double[]{.1/2, 1-0.1, .1/2}});
        testData.add(new Object[]{0.1, new double[]{0, .5, .5}, new double[]{0.5*0.1*0.5, 0.5*(1-0.1), 0.5*1+0.5*0.1/2}});
        testData.add(new Object[]{0.1, new double[]{0.5, 0.5, 0}, new double[]{.5+0.5*0.1*0.5, 0.5*(1-0.1), 0.5*0.1*0.5}});

        return testData.iterator();
    }

    @Test(dataProvider = "testGetLikelihoodsData")
    public void testGetLikelihoods(double pLoH, double[] underlyingLikelihood, double[] tumorLikelihood) throws Exception {
        HaplotypeProbabilities hp = new HaplotypeProbabilitiesFromGenotype(snp, hb, underlyingLikelihood[0], underlyingLikelihood[1], underlyingLikelihood[2]);

        HaplotypeProbabilities hpTumor = new HaplotypeProbabilityOfNormalGivenTumor(hp, pLoH);

        TestNGUtil.assertEqualDoubleArrays(hpTumor.getLikelihoods(), tumorLikelihood, 0.0001);

    }
}