package picard.vcf;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFFileReader;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.util.*;

import static picard.vcf.LiftoverVcf.leftAlignVariant;

/**
 * Test class for LiftoverVcf.
 * <p>
 * Created by ebanks on 8/11/15.
 */
public class LiftoverVcfTest extends CommandLineProgramTest {

    private static final File TEST_DATA_PATH = new File("testdata/picard/vcf/");
    private static final File CHAIN_FILE = new File(TEST_DATA_PATH, "test.over.chain");
    private static final File CHAIN_FILE_WITH_BAD_CONTIG = new File(TEST_DATA_PATH, "test.over.badContig.chain");
    private static final File REFERENCE_FILE = new File(TEST_DATA_PATH, "dummy.reference.fasta");
    private static final File OUTPUT_DATA_PATH = IOUtil.createTempDir("LiftoverVcfsTest", null);

    public String getCommandLineProgramName() {
        return LiftoverVcf.class.getSimpleName();
    }

    @AfterClass
    public void teardown() {
        IOUtil.deleteDirectoryTree(OUTPUT_DATA_PATH);
    }

    @DataProvider(name="liftoverReverseStrand")
    public Object[][] liftoverReverseStrand(){
        return new Object[][]{
                {"testLiftoverBiallelicIndels.vcf",3,0},
                {"testLiftoverMultiallelicIndels.vcf",0,2},
        };
    }


    @Test(dataProvider = "liftoverReverseStrand" )
    public void testReverseComplementedIndels(String filename, int expectedPassing, int expectedFailing) {
        final File liftOutputFile = new File(OUTPUT_DATA_PATH, "lift-delete-me.vcf");
        final File rejectOutputFile = new File(OUTPUT_DATA_PATH, "reject-delete-me.vcf");
        final File input = new File(TEST_DATA_PATH, filename);

        liftOutputFile.deleteOnExit();
        rejectOutputFile.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + liftOutputFile.getAbsolutePath(),
                "REJECT=" + rejectOutputFile.getAbsolutePath(),
                "CHAIN=" + CHAIN_FILE,
                "REFERENCE_SEQUENCE=" + REFERENCE_FILE,
                "CREATE_INDEX=false"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final VCFFileReader liftReader = new VCFFileReader(liftOutputFile, false);
        Assert.assertEquals(liftReader.iterator().stream().count(), expectedPassing, "The wrong number of variants was lifted over");

        final VCFFileReader rejectReader = new VCFFileReader(rejectOutputFile, false);
        Assert.assertEquals(rejectReader.iterator().stream().count(), expectedFailing, "The wrong number of variants was lifted over");
    }

    @DataProvider(name = "dataTestMissingContigInReference")
    public Object[][] dataTestHaplotypeProbabilitiesFromSequenceAddToProbs() {
        return new Object[][]{
                {false, LiftoverVcf.EXIT_CODE_WHEN_CONTIG_NOT_IN_REFERENCE},
                {true, 0}
        };
    }

    @Test(dataProvider = "dataTestMissingContigInReference")
    public void testMissingContigInReference(boolean warnOnMissingContext, int expectedReturnCode) {
        final File liftOutputFile = new File(OUTPUT_DATA_PATH, "lift-delete-me.vcf");
        final File rejectOutputFile = new File(OUTPUT_DATA_PATH, "reject-delete-me.vcf");
        final File input = new File(TEST_DATA_PATH, "testLiftoverUsingMissingContig.vcf");

        liftOutputFile.deleteOnExit();
        rejectOutputFile.deleteOnExit();

        // Test using WMC option
        final String[] argsWithWarnOnMissingContig = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + liftOutputFile.getAbsolutePath(),
                "REJECT=" + rejectOutputFile.getAbsolutePath(),
                "CHAIN=" + CHAIN_FILE_WITH_BAD_CONTIG,
                "REFERENCE_SEQUENCE=" + REFERENCE_FILE,
                "CREATE_INDEX=false",
                "WMC=" + warnOnMissingContext
        };
        Assert.assertEquals(runPicardCommandLine(argsWithWarnOnMissingContig), expectedReturnCode);
    }

    @DataProvider(name = "dataTestWriteOriginalPosition")
    public Object[][] dataTestWriteOriginalPosition() {
        return new Object[][]{
                {false},
                {true}
        };
    }

    @Test(dataProvider = "dataTestWriteOriginalPosition")
    public void testWriteOriginalPosition(boolean shouldWriteOriginalPosition) {
        final File liftOutputFile = new File(OUTPUT_DATA_PATH, "lift-delete-me.vcf");
        final File rejectOutputFile = new File(OUTPUT_DATA_PATH, "reject-delete-me.vcf");
        final File input = new File(TEST_DATA_PATH, "testLiftoverBiallelicIndels.vcf");

        liftOutputFile.deleteOnExit();
        rejectOutputFile.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + liftOutputFile.getAbsolutePath(),
                "REJECT=" + rejectOutputFile.getAbsolutePath(),
                "CHAIN=" + CHAIN_FILE,
                "REFERENCE_SEQUENCE=" + REFERENCE_FILE,
                "CREATE_INDEX=false",
                "WRITE_ORIGINAL_POSITION=" + shouldWriteOriginalPosition
        };

        runPicardCommandLine(args);

        try (VCFFileReader liftReader = new VCFFileReader(liftOutputFile, false)) {
            for (VariantContext vc : liftReader) {
                if (shouldWriteOriginalPosition) {
                    Assert.assertNotNull(vc.getAttribute(LiftoverVcf.ORIGINAL_CONTIG));
                    Assert.assertNotNull(vc.getAttribute(LiftoverVcf.ORIGINAL_START));
                } else {
                    Assert.assertFalse(vc.hasAttribute(LiftoverVcf.ORIGINAL_CONTIG));
                    Assert.assertFalse(vc.hasAttribute(LiftoverVcf.ORIGINAL_START));
                }
            }
        }
    }

    @DataProvider(name = "indelFlipData")
    public Iterator<Object[]> indelFlipData() {

        ReferenceSequence reference = new ReferenceSequence("chr1", 0,
                "CAAAAAAAAAACGTACGTACTCTCTCTCTACGT".getBytes());
        //       123456789 123456789 123456789 123

        Allele RefAAA = Allele.create("AAA", true);
        Allele RefCAA = Allele.create("CAA", true);
        Allele RefGTT = Allele.create("GTT", true);
        Allele RefACGT = Allele.create("ACGT", true);
        Allele RefAACG = Allele.create("AACG", true);
        Allele RefTTT = Allele.create("TTT", true);
        Allele RefCG = Allele.create("CG", true);
        Allele RefT = Allele.create("T", true);
        Allele RefA = Allele.create("A", true);
        Allele RefC = Allele.create("C", true);

        Allele A = Allele.create("A", false);
        Allele T = Allele.create("T", false);
        Allele C = Allele.create("C", false);
        Allele CG = Allele.create("CG", false);
        Allele CAA = Allele.create("CAA", false);
        Allele ACT = Allele.create("ACT", false);
        Allele TTT = Allele.create("TTT", false);
        Allele ATT = Allele.create("ATT", false);
        Allele TAG = Allele.create("TAG", false);
        Allele ACGT = Allele.create("ACGT", false);
        Allele AACG = Allele.create("AACG", false);

        final List<Object[]> tests = new ArrayList<>();

        final int CHAIN_SIZE = 540; // the length of the single chain in CHAIN_FILE

        VariantContextBuilder builder = new VariantContextBuilder().source("test1").chr("chr1");
        VariantContextBuilder result_builder = new VariantContextBuilder().source("test1").chr("chr1");

        // simple deletion
        // TTT*/T -> AAA*/A turns into left-aligned CAA*/C at position 1
        int start = 537;
        int stop = start + 2;
        builder.start(start).stop(stop).alleles(CollectionUtil.makeList(RefTTT, T));
        result_builder.start(1).stop(3).alleles(CollectionUtil.makeList(RefCAA, C));
        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        //simple insertion
        // T*/TTT -> A*/AAA -> turns into left-aligned C*/CAA at position 1
        stop = start;
        builder.source("test2").alleles(CollectionUtil.makeList(RefT, TTT)).stop(stop);
        result_builder.alleles(CollectionUtil.makeList(RefC, CAA)).start(1).stop(1);
        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        // non-simple deletion
        //  ACGT(T)*/A(T) -> AACG(T)*/A(T)
        //       position of 13 ^ in result
        start = CHAIN_SIZE - 13;
        stop = start + 3;

        builder.source("test3").start(start).stop(stop).alleles(CollectionUtil.makeList(RefACGT, A));
        result_builder.start(CHAIN_SIZE - stop).stop(CHAIN_SIZE - start).alleles(CollectionUtil.makeList(RefAACG, A));
        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        //  "CAAAAAAAAAACG---CGTACTCTCTCTCTACGT".getBytes());
        //  "CAAAAAAAAAACGacgCGTACTCTCTCTCTACGT".getBytes());

        // equivalent to

        //  "CAAAAAAAAA---ACGCGTACTCTCTCTCTACGT".getBytes());
        //  "CAAAAAAAAAACGacgCGTACTCTCTCTCTACGT".getBytes());

        //non-simple insertion
        // A(C)*/ACGT(C) -> G(T)*/GACG(T) -> gets moves 3 bases back due to left-aligning...
        //   position of 13 ^ in result
        // ...and becomes A->AACG at position 10

        start = CHAIN_SIZE - 13;
        stop = start;
        builder.source("test4").stop(stop).start(start).alleles(CollectionUtil.makeList(RefA, ACGT));
        result_builder.start(10).stop(10).alleles(CollectionUtil.makeList(RefA, AACG));
        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

//      // outside of chain (due to shifting of anchor point)
        start = stop = CHAIN_SIZE;
        builder.source("test5").stop(stop).start(start).alleles(CollectionUtil.makeList(RefA, ACGT));
        tests.add(new Object[]{builder.make(), reference, null});

//      // outside of chain
        start = stop = CHAIN_SIZE + 1;
        builder.source("test6").stop(stop).start(start).alleles(CollectionUtil.makeList(RefA, ACGT));
        tests.add(new Object[]{builder.make(), reference, null});

//      // MNP
        // GTT*(T)/ACGT(T) -> AAA(C)*/AACG(T) -> which is then normalized to A*/CG at position 11
        //               pos 11 ^ in the result
        start = CHAIN_SIZE - 11;
        stop = start + 2;
        builder.source("test7").stop(stop).start(start).alleles(CollectionUtil.makeList(RefGTT, ACGT));
        result_builder.start(11).stop(11).alleles(CollectionUtil.makeList(RefA, CG));
        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

//      // MNP
        // ACGT*(T)/ATT*(T) -> AACG(T)*/AAA(T) -> by normalization CG(T)*/A(T)
        //                 pos 13 ^ in the result
        start = CHAIN_SIZE - 13;
        stop = start + 3;
        builder.source("test8").stop(stop).start(start).alleles(CollectionUtil.makeList(RefACGT, ATT));
        result_builder.start(CHAIN_SIZE - stop + 2).stop(CHAIN_SIZE - stop + 3).alleles(CollectionUtil.makeList(RefCG, A));
        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        // needs left-aligning
        //  T->TAG       --> T(A)/TCT(A) -> by normalization A/ACT @ 19
        //    position 29    ^
        start = CHAIN_SIZE - 29;
        stop = start;
        builder.source("test9").stop(stop).start(start).alleles(CollectionUtil.makeList(RefT, TAG));
        result_builder.start(19).stop(19).alleles(CollectionUtil.makeList(RefA, ACT));
        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        return tests.iterator();
    }

    @Test(dataProvider = "indelFlipData")
    public void testFlipIndel(final VariantContext source, final ReferenceSequence reference, final VariantContext result) {

        final LiftOver liftOver = new LiftOver(CHAIN_FILE);

        final VariantContext flipped = LiftoverVcf.flipIndel(source, liftOver, reference);

        assertVcAreEqual(flipped, result);
    }

    @DataProvider(name = "leftAlignAllelesData")
    public Iterator<Object[]> leftAlignAllelesData() {

        ReferenceSequence reference = new ReferenceSequence("chr1", 0,
                "CAAAAAAAAAACGTACGTACTCTCTCTCTACGT".getBytes());
        //       123456789 123456789 123456789 123

        Allele RefG = Allele.create("G", true);
        Allele A = Allele.create("A", false);

        Allele RefA = Allele.create("A", true);
        Allele RefC = Allele.create("C", true);
        Allele RefAA = Allele.create("AA", true);
        Allele RefCA = Allele.create("CA", true);
        Allele RefCT = Allele.create("CT", true);
        Allele RefTC = Allele.create("TC", true);
        Allele RefACT = Allele.create("ACT", true);
        Allele RefCTCT = Allele.create("CTCT", true);
        Allele RefTCTC = Allele.create("TCTC", true);

        Allele AA = Allele.create("AA", false);
        Allele C = Allele.create("C", false);
        Allele CA = Allele.create("CA", false);
        Allele ACT = Allele.create("ACT", false);
        Allele CTCT = Allele.create("CTCT", false);
        Allele TCTC = Allele.create("TCTC", false);
        Allele CT = Allele.create("CT", false);
        Allele TC = Allele.create("TC", false);

        Allele T = Allele.create("T", false);

        final List<Object[]> tests = new ArrayList<>();

        VariantContextBuilder builder = new VariantContextBuilder().source("test1").chr("chr1");
        VariantContextBuilder result_builder = new VariantContextBuilder().source("test1").chr("chr1");

        // simple SNP
        // G*/A -> G/A
        int start = 13;
        int stop = start;
        builder.start(start).stop(stop).alleles(CollectionUtil.makeList(RefG, A));
        result_builder.start(stop).stop(start).alleles(CollectionUtil.makeList(RefG, A));
        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        for (start = 1; start <= reference.getBases().length; start++) {
            builder.start(start).stop(start);
            builder.alleles(CollectionUtil.makeList(
                    Allele.create(reference.getBaseString().substring(start - 1, start), true),
                    reference.getBaseString().charAt(start - 1) == 'A' ? T : A));

            tests.add(new Object[]{builder.make(), reference, builder.make()});
        }

        // AA/A in initial polyA repeat -> CA/C at the beginning
        result_builder.start(1).stop(2).alleles(CollectionUtil.makeList(RefCA, C));
        for (start = 2; start <= 11; start++) {
            builder.start(start).stop(start + 1);
            builder.alleles(CollectionUtil.makeList(
                    RefAA, A));

            tests.add(new Object[]{builder.make(), reference, result_builder.make()});
        }

        // A/AA in initial polyA repeat -> C/CA at the beginning
        result_builder.start(1).stop(1).alleles(CollectionUtil.makeList(RefC, CA));
        for (start = 2; start <= 11; start++) {
            builder.start(start).stop(start);
            builder.alleles(CollectionUtil.makeList(
                    RefA, AA));

            tests.add(new Object[]{builder.make(), reference, result_builder.make()});
        }

        //CT/CTCT -> A/ACT in CT repeat region
        result_builder.start(19).stop(19).alleles(CollectionUtil.makeList(RefA, ACT));
        for (start = 20; start <= 27; start += 2) {
            builder.start(start).stop(start + 1);
            builder.alleles(CollectionUtil.makeList(
                    RefCT, CTCT));

            tests.add(new Object[]{builder.make(), reference, result_builder.make()});
        }

        //TC/TCTC -> A/ACT in CT repeat region
        for (start = 21; start <= 29; start += 2) {
            builder.start(start).stop(start + 1);
            builder.alleles(CollectionUtil.makeList(
                    RefTC, TCTC));

            tests.add(new Object[]{builder.make(), reference, result_builder.make()});
        }

        //CTCT/CT -> ACT/A in CT repeat region
        result_builder.start(19).stop(21).alleles(CollectionUtil.makeList(RefACT, A));
        for (start = 20; start <= 27; start += 2) {
            builder.start(start).stop(start + 3);
            builder.alleles(CollectionUtil.makeList(
                    RefCTCT, CT));

            tests.add(new Object[]{builder.make(), reference, result_builder.make()});
        }

        //TCTC/TC-> ACT/A in CT repeat region
        for (start = 21; start <= 29; start += 2) {
            builder.start(start).stop(start + 3);
            builder.alleles(CollectionUtil.makeList(
                    RefTCTC, TC));

            tests.add(new Object[]{builder.make(), reference, result_builder.make()});
        }

        // for ease of reading, here's the reference sequence
        // "CAAAAAAAAAACGTACGTACTCTCTCTCTACGT"
        //  123456789 123456789 123456789 123

        result_builder.alleles("AACGT", "A").start(10).stop(14);
        for (start = 10; start < 17; start++) {
            for (stop = start + 4; stop < 20; stop++) {
                builder.alleles(
                        // -1 here due to reference string being 0-based.
                        reference.getBaseString().substring(start - 1, stop + 1 - 1),
                        reference.getBaseString().substring(start - 1, stop - 3 - 1)).start(start).stop(stop);
                tests.add(new Object[]{builder.make(), reference, result_builder.make()});
            }
        }

        // test vc with genotypes:

        result_builder.start(10).stop(14).alleles("AACGT", "A", "AACG", "AACGTACGT");
        builder.start(9).stop(18).alleles("AAACGTACGT", "AAACGT", "AAACGACGT", "AAACGTACGTACGT");
        Collection<Genotype> genotypes = new ArrayList<>();
        Collection<Genotype> results_genotypes = new ArrayList<>();
        GenotypeBuilder gBuilder = new GenotypeBuilder();

        for (int sample = 1; sample < 4; sample++) {
            gBuilder.name("sample" + sample)
                    .alleles(CollectionUtil.makeList(builder.getAlleles().get(0), builder.getAlleles().get(sample)));
            genotypes.add(gBuilder.make());
            gBuilder.alleles(CollectionUtil.makeList(result_builder.getAlleles().get(0), result_builder.getAlleles().get(sample)));
            results_genotypes.add(gBuilder.make());
        }
        gBuilder.name("sample_non_ref_het")
                .alleles(CollectionUtil.makeList(builder.getAlleles().get(1), builder.getAlleles().get(2)));
        genotypes.add(gBuilder.make());
        gBuilder.alleles(CollectionUtil.makeList(result_builder.getAlleles().get(1), result_builder.getAlleles().get(2)));
        results_genotypes.add(gBuilder.make());

        builder.genotypes(genotypes);
        result_builder.genotypes(results_genotypes);

        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        return tests.iterator();
    }

    @Test(dataProvider = "leftAlignAllelesData")
    public void testLeftAlignVariants(final VariantContext source, final ReferenceSequence reference, final VariantContext result) {

        assertVcAreEqual(leftAlignVariant(source, reference), result);
    }

    private void assertVcAreEqual(final VariantContext actual, final VariantContext expected) {

        if (expected == null) {
            Assert.assertNull(actual);
            return;
        }

        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.getContig(), expected.getContig());
        Assert.assertEquals(actual.getStart(), expected.getStart());
        Assert.assertEquals(actual.getEnd(), expected.getEnd());

        expected.getAlleles().sort(new AlleleComparator());
        actual.getAlleles().sort(new AlleleComparator());

        Assert.assertEquals(expected.getAlleles(), actual.getAlleles());
        assertGenotypeContextsAreEquals(actual.getGenotypes(), expected.getGenotypes());
    }

    static class AlleleComparator implements Comparator<Allele> {

        @Override
        public int compare(Allele o1, Allele o2) {
            int ret = (o1.isReference() ? 1 : 0) - (o2.isReference() ? 1 : 0);
            if (ret != 0) return ret;

            return o1.getBaseString().compareTo(o2.getBaseString());
        }
    }

    private void assertGenotypeContextsAreEquals(final GenotypesContext actual, final GenotypesContext expected) {
        if (expected == null) {
            Assert.assertNull(actual);
            return;
        }
        Assert.assertEquals(expected.getSampleNamesOrderedByName(), expected.getSampleNamesOrderedByName());

        for (String name : expected.getSampleNamesOrderedByName()) {
            Assert.assertEquals(expected.get(name).getAlleles(), actual.get(name).getAlleles());
        }
    }
}
