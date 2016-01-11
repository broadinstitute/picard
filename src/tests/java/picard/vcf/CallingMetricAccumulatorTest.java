package picard.vcf;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.FastGenotype;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import static org.testng.Assert.*;

/**
 * Created by farjoun on 12/26/15.
 */
public class CallingMetricAccumulatorTest {

    @DataProvider(name = "getSingletonSampleData")
    public Object[][] getSingletonSampleData() {
        final List<Object[]> retval = new ArrayList<>(10);

        final Allele ARef = Allele.create("A", true);
        final Allele G = Allele.create("G", false);
        final Allele C = Allele.create("C", false);

        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder().alleles(CollectionUtil.makeList(ARef, C));
        final VariantContextBuilder builder = new VariantContextBuilder();

        // one het
        final Genotype het = genotypeBuilder.name("het").make();
        builder.chr("1").start(1).stop(1).alleles(CollectionUtil.makeList(ARef, C, G)).genotypes(Collections.singletonList(het));
        retval.add(new Object[]{builder.make(), "het"});

        //a het and a hom ref
        final Genotype homref = genotypeBuilder.name("homref").alleles(CollectionUtil.makeList(ARef)).make();
        builder.genotypes(CollectionUtil.makeList(het, homref));
        retval.add(new Object[]{builder.make(), "het"});

        // a het, a homvar and a homref
        final Genotype homvar = genotypeBuilder.name("homvar").alleles(CollectionUtil.makeList(C)).make();
        builder.genotypes(CollectionUtil.makeList(het, homref, homvar));
        retval.add(new Object[]{builder.make(), null});

        // two hets and a homref
        final Genotype het2 = genotypeBuilder.name("het2").alleles(CollectionUtil.makeList(ARef, G)).make();
        builder.genotypes(CollectionUtil.makeList(het, homref, het2));
        retval.add(new Object[]{builder.make(), null});

        // a homvar
        builder.genotypes(CollectionUtil.makeList(homvar));
        retval.add(new Object[]{builder.make(), null});

        // a homvar, and a homref
        builder.genotypes(CollectionUtil.makeList(homvar, homref));
        retval.add(new Object[]{builder.make(), null});

        // two homrefs
        final Genotype homref2 = genotypeBuilder.name("homref2").alleles(CollectionUtil.makeList(ARef)).make();
        builder.genotypes(CollectionUtil.makeList(homref, homref2));
        retval.add(new Object[]{builder.make(), null});

        return retval.toArray(new Object[retval.size()][]);
    }

    @Test(dataProvider = "getSingletonSampleData")
    public void testGetSingletonSample(final VariantContext vc, final String sample) throws Exception {
        Assert.assertEquals(CallingMetricAccumulator.getSingletonSample(vc), sample);
    }
}
