package picard.util;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.vcf.VcfTestUtils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by farjoun on 1/16/18.
 */
public class LiftoverUtilsTest {

    static final Allele A = Allele.create("A", false);
    static final Allele C = Allele.create("C", false);
    static final Allele RefA = Allele.create("A", true);
    static final Allele RefC = Allele.create("C", true);

    static String AF = "AF";
    static String MAX_AF = "MAX_AF";

    static final List<String> annotationsToDrop = Collections.singletonList(MAX_AF);
    static final List<String> annotationsToSwap = Collections.singletonList(AF);

    @DataProvider
    public Object[][] swapRefAltData() {

        String testName = "test1";

        final List<Object[]> tests = new ArrayList<>();

        VariantContextBuilder builder = new VariantContextBuilder(testName, "chr1", 2, 2, Arrays.asList(RefA, C));
        VariantContextBuilder resultBuilder = new VariantContextBuilder(testName, "chr1", 2, 2, Arrays.asList(RefC, A));
        resultBuilder.attribute("SwappedAlleles", true);

        tests.add(new Object[]{builder.make(), resultBuilder.make()});

        testName = "test2";
        builder.source(testName);
        builder.attribute(AF, .2);

        resultBuilder.source(testName);
        resultBuilder.attribute(AF, .8);
        tests.add(new Object[]{builder.make(), resultBuilder.make()});

        testName = "test3";
        builder.source(testName);
        builder.attribute(AF, ".");

        resultBuilder.source(testName);
        resultBuilder.attribute(AF, ".");
        tests.add(new Object[]{builder.make(), resultBuilder.make()});

        testName = "test4";
        builder.source(testName);
        builder.attribute(AF, 1);

        resultBuilder.source(testName);
        resultBuilder.attribute(AF, 0D);
        tests.add(new Object[]{builder.make(), resultBuilder.make()});

        testName = "test5";
        builder.source(testName);
        builder.attribute(MAX_AF, 1);

        resultBuilder.source(testName);
        tests.add(new Object[]{builder.make(), resultBuilder.make()});

        testName = "test5";
        builder.source(testName);
        builder.attribute("AF_MAX_shouldnt_be_dropped", 1);

        resultBuilder.source(testName);
        resultBuilder.attribute("AF_MAX_shouldnt_be_dropped", 1);
        tests.add(new Object[]{builder.make(), resultBuilder.make()});

        return tests.toArray(new Object[0][]);
    }

    @Test(dataProvider = "swapRefAltData")
    public void testSwapRefAlt(final VariantContext swapMe, final VariantContext expected) {
        VcfTestUtils.assertEquals(LiftoverUtils.swapRefAlt(swapMe, annotationsToSwap, annotationsToDrop), expected);
    }

    @DataProvider
    public Iterator<Object[]> allelesToStringData() {

        final List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{Arrays.asList(Allele.create("A", true), Allele.create("G")), Arrays.asList("A", "G")});
        tests.add(new Object[]{Arrays.asList(Allele.create("A", false), Allele.create("G")), Arrays.asList("A", "G")});
        tests.add(new Object[]{Arrays.asList(Allele.create("T", true), Allele.create("T")), Arrays.asList("T", "T")});
        tests.add(new Object[]{Arrays.asList(Allele.create("*", false), Allele.create("G")), Arrays.asList("*", "G")});
        tests.add(new Object[]{Arrays.asList(Allele.create("ATG", true), Allele.create("G")), Arrays.asList("ATG", "G")});
        tests.add(new Object[]{Arrays.asList(Allele.create("A", true), Allele.create("GGGTGT")), Arrays.asList("A", "GGGTGT")});
        tests.add(new Object[]{Arrays.asList(Allele.create("A", false), Allele.create("A")), Arrays.asList("A", "A")});

        tests.add(new Object[]{Arrays.asList(Allele.NO_CALL, Allele.NO_CALL), Arrays.asList(Allele.NO_CALL_STRING, Allele.NO_CALL_STRING)});
        tests.add(new Object[]{Arrays.asList(Allele.SPAN_DEL, Allele.create("A", true)), Arrays.asList(Allele.SPAN_DEL_STRING, "A")});
        tests.add(new Object[]{Arrays.asList(Allele.create("A", true), Allele.create(".")), Arrays.asList("A", Allele.NO_CALL_STRING)});
        tests.add(new Object[]{Arrays.asList(Allele.create("TT", true), Allele.create(".")), Arrays.asList("TT", Allele.NO_CALL_STRING)});

        return tests.iterator();
    }

    @Test(dataProvider = "allelesToStringData")
    public void testAllelesToString(final List<Allele> input, final List<String> output) {
        Assert.assertEquals(LiftoverUtils.allelesToStringList(input), output);

        //these should back-convert into the same alleles.
        List<Allele> restoredAlleles = output.stream().map(Allele::create).collect(Collectors.toList());
        for (int i = 0;i<input.size();i++){
            Assert.assertTrue(restoredAlleles.get(i).basesMatch(input.get(i)));
        }
    }
}