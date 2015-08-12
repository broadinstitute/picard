package picard.vcf;

import htsjdk.variant.variantcontext.*;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Test class for LiftoverVcf.
 *
 * Created by ebanks on 8/11/15.
 */
public class LiftoverVcfTest {

    @Test
    public void testFixReverseComplementedGenotypes() {

        final Allele refA = Allele.create("A", true);
        final Allele altC = Allele.create("C", false);
        final GenotypesContext originalGenotypes = GenotypesContext.create(3);
        originalGenotypes.add(new GenotypeBuilder("homref").alleles(Arrays.asList(refA, refA)).make());
        originalGenotypes.add(new GenotypeBuilder("het").alleles(Arrays.asList(refA, altC)).make());
        originalGenotypes.add(new GenotypeBuilder("homvar").alleles(Arrays.asList(altC, altC)).make());

        final Allele refT = Allele.create("T", true);
        final Allele altG = Allele.create("G", false);
        final GenotypesContext expectedGenotypes = GenotypesContext.create(3);
        expectedGenotypes.add(new GenotypeBuilder("homref").alleles(Arrays.asList(refT, refT)).make());
        expectedGenotypes.add(new GenotypeBuilder("het").alleles(Arrays.asList(refT, altG)).make());
        expectedGenotypes.add(new GenotypeBuilder("homvar").alleles(Arrays.asList(altG, altG)).make());

        final Map<Allele, Allele> reverseComplementAlleleMap = new HashMap<Allele, Allele>(2);
        reverseComplementAlleleMap.put(refA, refT);
        reverseComplementAlleleMap.put(altC, altG);
        final GenotypesContext actualGenotypes = LiftoverVcf.fixGenotypes(originalGenotypes, reverseComplementAlleleMap);

        for ( final String sample : Arrays.asList("homref", "het", "homvar") ) {
            final List<Allele> expected = expectedGenotypes.get(sample).getAlleles();
            final List<Allele> actual = actualGenotypes.get(sample).getAlleles();
            Assert.assertEquals(expected.get(0), actual.get(0));
            Assert.assertEquals(expected.get(1), actual.get(1));
        }
    }
}
