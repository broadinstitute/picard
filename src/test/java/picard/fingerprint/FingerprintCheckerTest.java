package picard.fingerprint;

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * Created by farjoun on 8/27/15.
 */
public class FingerprintCheckerTest {

    private final double maf = 0.4;
    private final Snp snp = new Snp("test", "chr1", 1, (byte) 'A', (byte) 'C', maf, Collections.singletonList("dummy"));
    private final HaplotypeBlock hb = new HaplotypeBlock(maf);

    @BeforeClass
    public void setup() {
        hb.addSnp(snp);
    }

    @Test
    public void testRandomSublist() throws Exception {

        final List<Integer> list = new ArrayList<>();
        list.add(1);
        list.add(2);
        list.add(3);

        Assert.assertEquals(list, FingerprintChecker.randomSublist(list, 3));
        Assert.assertEquals(list, FingerprintChecker.randomSublist(list, 4));
        Assert.assertEquals(FingerprintChecker.randomSublist(list, 2).size(), 2);
    }

    @DataProvider(name = "pLoH")
    public Iterator<Object[]> pLohData() {
        final List<Object[]> listOfDoubles = new ArrayList<>();

        for (int i = 1; i < 20; i++) {
            listOfDoubles.add(new Object[]{i / 40D});
        }
        return listOfDoubles.iterator();
    }

    @Test(dataProvider = "pLoH")
    public void testMatchResults(final double pLoH) {

        final Fingerprint fpObserved = new Fingerprint("test", null, "noop");
        final Fingerprint fpExpected = new Fingerprint("test", null, "noop");

        final HaplotypeProbabilities hpHet = new HaplotypeProbabilitiesFromGenotype(snp, hb, 0.0001, 1.0, 0.0001);
        final HaplotypeProbabilities hpHomRef = new HaplotypeProbabilitiesFromGenotype(snp, hb, 1.0, 0.00001, 0.000000001);

        // Expected is a het
        fpExpected.add(hpHet);

        // Observed is a hom, so possible scenario is that observed is tumor, and expected is normal
        fpObserved.add(hpHomRef);

        // get match results using pLOD
        final MatchResults mr = FingerprintChecker.calculateMatchResults(fpObserved, fpExpected, 0.01, pLoH);

        // make sure that it's more likely to be the same sample, if the observed is "tumor" and the expected is "normal"
        Assert.assertTrue(mr.getLodTN() > mr.getLOD());

        // make sure that the regular LOD is negative (we're comparing a HET to a HOM)
        Assert.assertTrue(mr.getLOD() < 0);

        // make sure that it's more likely to be tumor/normal rather than normal/tumor
        // (a hom normal isn't expected to be measured as a het in the tumor)
        Assert.assertTrue(mr.getLodTN() > mr.getLodNT());
    }
}