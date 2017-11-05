package picard.util;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFFileReader;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.vcf.LiftoverVcf;
import picard.vcf.VcfTestUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

/**
 * Test class for LiftoverVcf.
 * <p>
 * Created by ebanks on 8/11/15.
 */
public class LiftoverVcfTest extends CommandLineProgramTest {

    private static final File TEST_DATA_PATH = new File("testdata/picard/vcf/");
    private static final File CHAIN_FILE = new File(TEST_DATA_PATH, "test.over.chain");
    private static final File TWO_INTERVAL_CHAIN_FILE = new File(TEST_DATA_PATH, "test.two.block.over.chain");

    private static final File CHAIN_FILE_WITH_BAD_CONTIG = new File(TEST_DATA_PATH, "test.over.badContig.chain");
    private static final File REFERENCE_FILE = new File(TEST_DATA_PATH, "dummy.reference.fasta");
    private static final File TWO_INTERVALS_REFERENCE_FILE = new File(TEST_DATA_PATH, "dummy.two.block.reference.fasta");

    private static final File OUTPUT_DATA_PATH = IOUtil.createTempDir("LiftoverVcfsTest", null);

    public String getCommandLineProgramName() {
        return LiftoverVcf.class.getSimpleName();
    }

    @AfterClass
    public void teardown() {
        IOUtil.deleteDirectoryTree(OUTPUT_DATA_PATH);
    }

    @DataProvider(name = "liftoverReverseStrand")
    public Object[][] liftoverReverseStrand() {
        return new Object[][]{
                {"testLiftoverBiallelicIndels.vcf", 3, 0},
                {"testLiftoverMultiallelicIndels.vcf", 0, 2},
                {"testLiftoverFailingVariants.vcf", 0, 3},
        };
    }

    @Test(dataProvider = "liftoverReverseStrand")
    public void testReverseComplementedIndels(final String filename, final int expectedPassing, final int expectedFailing) {
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
    public void testMissingContigInReference(final boolean warnOnMissingContext, final int expectedReturnCode) {
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
    public void testWriteOriginalPosition(final boolean shouldWriteOriginalPosition) {
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

        try (final VCFFileReader liftReader = new VCFFileReader(liftOutputFile, false)) {
            for (final VariantContext vc : liftReader) {
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

        final ReferenceSequence reference = new ReferenceSequence("chr1", 0,
                "CAAAAAAAAAACGTACGTACTCTCTCTCTACGT".getBytes());
        //       123456789 123456789 123456789 123

        final Allele RefCAA = Allele.create("CAA", true);
        final Allele RefGTT = Allele.create("GTT", true);
        final Allele RefACGT = Allele.create("ACGT", true);
        final Allele RefAACG = Allele.create("AACG", true);
        final Allele RefTTT = Allele.create("TTT", true);
        final Allele RefCG = Allele.create("CG", true);
        final Allele RefT = Allele.create("T", true);
        final Allele RefA = Allele.create("A", true);
        final Allele RefC = Allele.create("C", true);
        final Allele RefG = Allele.create("G", true);

        final Allele A = Allele.create("A", false);
        final Allele T = Allele.create("T", false);
        final Allele C = Allele.create("C", false);
        final Allele CG = Allele.create("CG", false);
        final Allele GA = Allele.create("GA", false);
        final Allele TC = Allele.create("TC", false);
        final Allele CAA = Allele.create("CAA", false);
        final Allele ACT = Allele.create("ACT", false);
        final Allele TTT = Allele.create("TTT", false);
        final Allele ATT = Allele.create("ATT", false);
        final Allele GTT = Allele.create("GTT", false);
        final Allele AAC = Allele.create("AAC", false);
        final Allele TAG = Allele.create("TAG", false);
        final Allele ACGT = Allele.create("ACGT", false);
        final Allele AACG = Allele.create("AACG", false);

        final List<Object[]> tests = new ArrayList<>();

        final int CHAIN_SIZE = 540; // the length of the single chain in CHAIN_FILE

        final VariantContextBuilder builder = new VariantContextBuilder().source("test1").chr("chr1");
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder("test1");
        final GenotypeBuilder resultGenotypeBuilder = new GenotypeBuilder("test1");
        final VariantContextBuilder result_builder = new VariantContextBuilder().source("test1").chr("chr1");

        // simple deletion
        // TTT*/T -> AAA*/A turns into left-aligned CAA*/C at position 1
        int start = CHAIN_SIZE - 3;
        int stop = start + 2;
        builder.start(start).stop(stop).alleles(CollectionUtil.makeList(RefTTT, T));
        result_builder.start(1).stop(3).alleles(CollectionUtil.makeList(RefCAA, C));

        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());

        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        //simple insertion
        // T*/TTT -> A*/AAA -> turns into left-aligned C*/CAA at position 1
        stop = start;
        builder.source("test2").alleles(CollectionUtil.makeList(RefT, TTT)).stop(stop);
        result_builder.alleles(CollectionUtil.makeList(RefC, CAA)).start(1).stop(1);
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());

        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        builder.noGenotypes();
        result_builder.noGenotypes();

        // non-simple deletion
        //  ACGT(T)*/A(T) -> AACG(T)*/A(T)
        //       position of 13 ^ in result
        start = CHAIN_SIZE - 13;
        stop = start + 3;

        builder.source("test3").start(start).stop(stop).alleles(CollectionUtil.makeList(RefACGT, A));
        result_builder.start(CHAIN_SIZE - stop).stop(CHAIN_SIZE - start).alleles(CollectionUtil.makeList(RefAACG, A));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        //  "CAAAAAAAAAACG---CGTACTCTCTCTCTACGT" -- Allele A
        //  "CAAAAAAAAAACGacgCGTACTCTCTCTCTACGT" -- Allele B

        // is equivalent to

        //  "CAAAAAAAAA---ACGCGTACTCTCTCTCTACGT" -- Allele A
        //  "CAAAAAAAAAACGacgCGTACTCTCTCTCTACGT" -- Allele B

        // non-simple insertion
        // A(C)*/ACGT(C) -> G(T)*/GACG(T) -> gets moves 3 bases back due to left-aligning...
        //   position of 13 ^ in result
        // ...and becomes A->AACG at position 10

        start = CHAIN_SIZE - 13;
        stop = start;
        builder.source("test4").stop(stop).start(start).alleles(CollectionUtil.makeList(RefA, ACGT));
        result_builder.start(10).stop(10).alleles(CollectionUtil.makeList(RefA, AACG));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        // just outside of chain & contig, testing that we do not read into negative indices
        // or reference
        start = stop = CHAIN_SIZE;
        builder.source("test5").stop(stop).start(start).alleles(CollectionUtil.makeList(RefG, GTT));
        result_builder.start(1).stop(1).alleles(CollectionUtil.makeList(RefC, AAC));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        // outside of chain
        start = stop = CHAIN_SIZE + 1;
        builder.source("test6").stop(stop).start(start).alleles(CollectionUtil.makeList(RefA, ACGT));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{builder.make(), reference, null});

        // MNP
        // GTT*(T)/ACGT(T) -> AAA(C)*/AACG(T) -> which is then normalized to A*/CG at position 11
        //               pos 11 ^ in the result
        start = CHAIN_SIZE - 11;
        stop = start + 2;
        builder.source("test7").stop(stop).start(start).alleles(CollectionUtil.makeList(RefGTT, ACGT));
        result_builder.start(11).stop(11).alleles(CollectionUtil.makeList(RefA, CG));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        // MNP
        // ACGT*(T)/ATT*(T) -> AACG(T)*/AAA(T) -> by normalization CG(T)*/A(T)
        //                 pos 13 ^ in the result
        start = CHAIN_SIZE - 13;
        stop = start + 3;
        builder.source("test8").stop(stop).start(start).alleles(CollectionUtil.makeList(RefACGT, ATT));
        result_builder.start(CHAIN_SIZE - stop + 2).stop(CHAIN_SIZE - stop + 3).alleles(CollectionUtil.makeList(RefCG, A));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        // needs left-aligning
        //  T->TAG       --> T(A)/TCT(A) -> by normalization A/ACT @ 19
        //    position 29    ^
        start = CHAIN_SIZE - 29;
        stop = start;
        builder.source("test9").stop(stop).start(start).alleles(CollectionUtil.makeList(RefT, TAG));
        result_builder.start(19).stop(19).alleles(CollectionUtil.makeList(RefA, ACT));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        // insertion at end of section
        // a test that converts the initial C to a AC which requires
        // code in LiftOver::extendOneBase(., false)
        //
        // G*(.)/GA(.) -> .(C)/.T(C) but that's not legal, so
        // -> C/TC
        start = CHAIN_SIZE;
        stop = CHAIN_SIZE;

        builder.start(start).stop(stop).alleles(CollectionUtil.makeList(RefG, GA));
        result_builder.start(1).stop(1).alleles(CollectionUtil.makeList(RefC, TC));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        // insertion at end of section
        // a test that converts the initial C to a AC which requires
        // code in LiftOver::extendOneBase(., false)
        //
        // G*(.)/GG(.) -> .(C)/.C(C) but that's not legal, so
        // -> C/CC
        start = CHAIN_SIZE;
        stop = CHAIN_SIZE;

        builder.start(start).stop(stop).alleles(CollectionUtil.makeList(RefG, Allele.create("GG")));
        result_builder.start(1).stop(1).alleles(CollectionUtil.makeList(RefC, Allele.create("CC")));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        // insertion at end of section
        // a test that converts the initial C to a AC which requires
        // code in LiftOver::extendOneBase(., false)
        //
        // improperly aligned
        // TTTT(G)->TTTTG(G) -> C/CC
        start = CHAIN_SIZE - 4;
        stop = CHAIN_SIZE - 1;

        builder.start(start).stop(stop).alleles(CollectionUtil.makeList(Allele.create("TTTT", true), Allele.create("TTTTG")));
        result_builder.start(1).stop(1).alleles(CollectionUtil.makeList(RefC, Allele.create("CC")));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        return tests.iterator();
    }

    @Test(dataProvider = "indelFlipData")
    public void testFlipIndel(final VariantContext source, final ReferenceSequence reference, final VariantContext result) {

        final LiftOver liftOver = new LiftOver(CHAIN_FILE);
        final Interval originalLocus = new Interval(source.getContig(), source.getStart(), source.getEnd());
        final Interval target = liftOver.liftOver(originalLocus);
        if (target != null && !target.isNegativeStrand()){
            throw new RuntimeException("not reversed");
        }

        final VariantContext flipped = LiftoverUtils.liftVariant(source, target, reference, false);

        VcfTestUtils.assertEquals(flipped, result);
    }

    @DataProvider(name = "leftAlignAllelesData")
    public Iterator<Object[]> leftAlignAllelesData() {

        final ReferenceSequence reference = new ReferenceSequence("chr1", 0,
                "CAAAAAAAAAACGTACGTACTCTCTCTCTACGT".getBytes());
        //       123456789 123456789 123456789 123

        final Allele RefG = Allele.create("G", true);
        final Allele A = Allele.create("A", false);

        final Allele RefA = Allele.create("A", true);
        final Allele RefC = Allele.create("C", true);
        final Allele RefAA = Allele.create("AA", true);
        final Allele RefCA = Allele.create("CA", true);
        final Allele RefCT = Allele.create("CT", true);
        final Allele RefTC = Allele.create("TC", true);
        final Allele RefACT = Allele.create("ACT", true);
        final Allele RefCTCT = Allele.create("CTCT", true);
        final Allele RefTCTC = Allele.create("TCTC", true);

        final Allele AA = Allele.create("AA", false);
        final Allele C = Allele.create("C", false);
        final Allele CA = Allele.create("CA", false);
        final Allele ACT = Allele.create("ACT", false);
        final Allele CTCT = Allele.create("CTCT", false);
        final Allele TCTC = Allele.create("TCTC", false);
        final Allele CT = Allele.create("CT", false);
        final Allele TC = Allele.create("TC", false);

        final Allele T = Allele.create("T", false);

        final List<Object[]> tests = new ArrayList<>();

        final VariantContextBuilder builder = new VariantContextBuilder().source("test1").chr("chr1");
        final VariantContextBuilder result_builder = new VariantContextBuilder().source("test1").chr("chr1");
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder("test1");
        final GenotypeBuilder resultGenotypeBuilder = new GenotypeBuilder("test1");

        // left aligning at the edge of the reference
        // simple SNP
        // CAAA*/CCAAA -> C->CC
        int start = 1;
        int stop = 4;
        builder.start(start).stop(stop).alleles(CollectionUtil.makeList(Allele.create("CAAA", true), Allele.create("CCAAA")));

        result_builder.start(1).stop(1).alleles(CollectionUtil.makeList(RefC, Allele.create("CC")));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());

        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        // simple SNP
        // G*/A -> G/A
        start = 13;
        stop = start;
        builder.source("test1_5");
        builder.start(start).stop(stop).alleles(CollectionUtil.makeList(RefG, A));
        result_builder.start(stop).stop(start).alleles(CollectionUtil.makeList(RefG, A));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        builder.source("test2");
        for (start = 1; start <= reference.getBases().length; start++) {
            builder.start(start).stop(start);
            builder.alleles(CollectionUtil.makeList(
                    Allele.create(reference.getBaseString().substring(start - 1, start), true),
                    reference.getBaseString().charAt(start - 1) == 'A' ? T : A));

            genotypeBuilder.alleles(builder.getAlleles());
            resultGenotypeBuilder.alleles(result_builder.getAlleles());
            builder.genotypes(genotypeBuilder.make());
            result_builder.genotypes(resultGenotypeBuilder.make());
            tests.add(new Object[]{builder.make(), reference, builder.make()});
        }

        // AA/A in initial polyA repeat -> CA/C at the beginning
        result_builder.start(1).stop(2).alleles(CollectionUtil.makeList(RefCA, C));
        builder.source("test3");
        for (start = 2; start <= 11; start++) {
            builder.start(start).stop(start + 1);
            builder.alleles(CollectionUtil.makeList(
                    RefAA, A));

            genotypeBuilder.alleles(builder.getAlleles());
            resultGenotypeBuilder.alleles(result_builder.getAlleles());
            builder.genotypes(genotypeBuilder.make());
            result_builder.genotypes(resultGenotypeBuilder.make());
            tests.add(new Object[]{builder.make(), reference, result_builder.make()});
        }

        // A/AA in initial polyA repeat -> C/CA at the beginning
        result_builder.start(1).stop(1).alleles(CollectionUtil.makeList(RefC, CA));
        builder.source("test4");
        for (start = 2; start <= 11; start++) {
            builder.start(start).stop(start);
            builder.alleles(CollectionUtil.makeList(
                    RefA, AA));

            genotypeBuilder.alleles(builder.getAlleles());
            resultGenotypeBuilder.alleles(result_builder.getAlleles());
            builder.genotypes(genotypeBuilder.make());
            result_builder.genotypes(resultGenotypeBuilder.make());
            tests.add(new Object[]{builder.make(), reference, result_builder.make()});
        }

        //CT/CTCT -> A/ACT in CT repeat region
        result_builder.start(19).stop(19).alleles(CollectionUtil.makeList(RefA, ACT));
        builder.source("test5");
        for (start = 20; start <= 27; start += 2) {
            builder.start(start).stop(start + 1);
            builder.alleles(CollectionUtil.makeList(
                    RefCT, CTCT));

            genotypeBuilder.alleles(builder.getAlleles());
            resultGenotypeBuilder.alleles(result_builder.getAlleles());
            builder.genotypes(genotypeBuilder.make());
            result_builder.genotypes(resultGenotypeBuilder.make());
            tests.add(new Object[]{builder.make(), reference, result_builder.make()});
        }

        builder.source("test6");
        //TC/TCTC -> A/ACT in CT repeat region
        for (start = 21; start <= 29; start += 2) {
            builder.start(start).stop(start + 1);
            builder.alleles(CollectionUtil.makeList(
                    RefTC, TCTC));

            genotypeBuilder.alleles(builder.getAlleles());
            resultGenotypeBuilder.alleles(result_builder.getAlleles());
            builder.genotypes(genotypeBuilder.make());
            result_builder.genotypes(resultGenotypeBuilder.make());
            tests.add(new Object[]{builder.make(), reference, result_builder.make()});
        }

        //CTCT/CT -> ACT/A in CT repeat region
        builder.source("test7");
        result_builder.start(19).stop(21).alleles(CollectionUtil.makeList(RefACT, A));
        for (start = 20; start <= 27; start += 2) {
            builder.start(start).stop(start + 3);
            builder.alleles(CollectionUtil.makeList(
                    RefCTCT, CT));

            genotypeBuilder.alleles(builder.getAlleles());
            resultGenotypeBuilder.alleles(result_builder.getAlleles());
            builder.genotypes(genotypeBuilder.make());
            result_builder.genotypes(resultGenotypeBuilder.make());
            tests.add(new Object[]{builder.make(), reference, result_builder.make()});
        }

        builder.source("test8");
        //TCTC/TC-> ACT/A in CT repeat region
        for (start = 21; start <= 29; start += 2) {
            builder.start(start).stop(start + 3);
            builder.alleles(CollectionUtil.makeList(
                    RefTCTC, TC));

            genotypeBuilder.alleles(builder.getAlleles());
            resultGenotypeBuilder.alleles(result_builder.getAlleles());
            builder.genotypes(genotypeBuilder.make());
            result_builder.genotypes(resultGenotypeBuilder.make());
            tests.add(new Object[]{builder.make(), reference, result_builder.make()});
        }

        // for ease of reading, here's the reference sequence
        // "CAAAAAAAAAACGTACGTACTCTCTCTCTACGT"
        //  123456789 123456789 123456789 123

        builder.source("test9");
        result_builder.alleles("AACGT", "A").start(10).stop(14);
        for (start = 10; start < 17; start++) {
            for (stop = start + 4; stop < 20; stop++) {
                builder.alleles(
                        // -1 here due to reference string being 0-based.
                        reference.getBaseString().substring(start - 1, stop + 1 - 1),
                        reference.getBaseString().substring(start - 1, stop - 3 - 1)).start(start).stop(stop);
                genotypeBuilder.alleles(builder.getAlleles());
                resultGenotypeBuilder.alleles(result_builder.getAlleles());
                builder.genotypes(genotypeBuilder.make());
                result_builder.genotypes(resultGenotypeBuilder.make());
                tests.add(new Object[]{builder.make(), reference, result_builder.make()});
            }
        }

        // test vc with genotypes:

        builder.source("test10");
        result_builder.start(10).stop(14).alleles("AACGT", "A", "AACG", "AACGTACGT");
        builder.start(9).stop(18).alleles("AAACGTACGT", "AAACGT", "AAACGACGT", "AAACGTACGTACGT");
        final Collection<Genotype> genotypes = new ArrayList<>();
        final Collection<Genotype> results_genotypes = new ArrayList<>();
        final GenotypeBuilder gBuilder = new GenotypeBuilder();

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

        builder.source("test12");
        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        return tests.iterator();
    }

    @Test(dataProvider = "leftAlignAllelesData")
    public void testLeftAlignVariants(final VariantContext source, final ReferenceSequence reference, final VariantContext result) {
        VariantContextBuilder vcb = new VariantContextBuilder(source);

        LiftoverUtils.leftAlignVariant(vcb, source.getStart(), source.getEnd(), source.getAlleles(), reference);
        vcb.genotypes(LiftoverUtils.fixGenotypes(source.getGenotypes(), source.getAlleles(), vcb.getAlleles()));

        VcfTestUtils.assertEquals(vcb.make(), result);
    }

    @DataProvider(name = "indelNoFlipData")
    public Iterator<Object[]> indelNoFlipData() {

        final VariantContextBuilder builder = new VariantContextBuilder().source("test1").chr("chr1");
        final VariantContextBuilder result_builder = new VariantContextBuilder().source("test1").chr("chr1");
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder("test1");
        final GenotypeBuilder resultGenotypeBuilder = new GenotypeBuilder("test1");
        final List<Object[]> tests = new ArrayList<>();

        // some more tests with a more complicated chain File. this one has 2 intervals in the relevant block
        // the cigar string would be 540M5I500M5D40M (if chains were written using cigar strings....)

        final Allele AAAARef = Allele.create("AAAA", true);
        final Allele AAA = Allele.create("AAA", false);
        final Allele CAAARef = Allele.create("CAAA", true);
        final Allele CAA = Allele.create("CAA", false);
        final Allele ARef = Allele.create("A", true);
        final Allele A = Allele.create("A", false);
        final Allele CRef = Allele.create("C", true);
        final Allele T = Allele.create("T", false);

        int start = 1;
        int stop = 4;
        int offset;

        final ReferenceSequence twoIntervalChainReference = new FastaSequenceFile(TWO_INTERVALS_REFERENCE_FILE, false).nextSequence();

        final LiftOver liftOver = new LiftOver(TWO_INTERVAL_CHAIN_FILE);

        // trivial snp
        builder.source("test1");
        builder.start(1).stop(1).alleles(CollectionUtil.makeList(CRef, A));
        result_builder.start(1).stop(1).alleles(CollectionUtil.makeList(CRef, A));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{liftOver, twoIntervalChainReference, builder.make(), result_builder.make()});

        // trivial case indel
        builder.source("test2");
        builder.start(start).stop(stop).alleles(CollectionUtil.makeList(CAAARef, CAA));
        result_builder.start(1).stop(4).alleles(CollectionUtil.makeList(CAAARef, CAA));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{liftOver, twoIntervalChainReference, builder.make(), result_builder.make()});

        // near end of interval indel
        builder.source("test3");
        builder.start(537).stop(540).alleles(CollectionUtil.makeList(AAAARef, AAA));
        result_builder.start(537).stop(540).alleles(CollectionUtil.makeList(AAAARef, AAA));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{liftOver, twoIntervalChainReference, builder.make(), result_builder.make()});

        // near end of interval snp
        builder.source("test4");
        builder.start(540).stop(540).alleles(CollectionUtil.makeList(ARef, T));
        result_builder.start(540).stop(540).alleles(CollectionUtil.makeList(ARef, T));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{liftOver, twoIntervalChainReference, builder.make(), result_builder.make()});

        // straddling chains indel
        builder.source("test5");
        builder.start(538).stop(541).alleles(CollectionUtil.makeList(AAAARef, AAA));
        result_builder.start(537).stop(540).alleles(CollectionUtil.makeList(AAAARef, AAA));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{liftOver, twoIntervalChainReference, builder.make(), null});

        // near start of second interval snp
        builder.source("test6");
        start = 541;
        offset = 5;
        builder.start(start).stop(start).alleles(CollectionUtil.makeList(ARef, T));
        result_builder.start(start + offset).stop(start + offset).alleles(CollectionUtil.makeList(ARef, T));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{liftOver, twoIntervalChainReference, builder.make(), result_builder.make()});

        // near start of second interval indel
        builder.source("test7");
        builder.start(start).stop(start).alleles(CollectionUtil.makeList(ARef, T));
        result_builder.start(start + offset).stop(start + offset).alleles(CollectionUtil.makeList(ARef, T));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{liftOver, twoIntervalChainReference, builder.make(), result_builder.make()});

        // near end of second interval snp
        builder.source("test8");
        start = 1040;
        offset = 5;
        builder.start(start).stop(start).alleles(CollectionUtil.makeList(ARef, T));
        result_builder.start(start + offset).stop(start + offset).alleles(CollectionUtil.makeList(ARef, T));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{liftOver, twoIntervalChainReference, builder.make(), result_builder.make()});

        // near end of second interval indel
        builder.source("test9");
        start = 1037;
        stop = 1040;
        builder.start(start).stop(stop).alleles(CollectionUtil.makeList(AAAARef, AAA));
        result_builder.start(start + offset).stop(stop + offset).alleles(CollectionUtil.makeList(AAAARef, AAA));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{liftOver, twoIntervalChainReference, builder.make(), result_builder.make()});

        // straddling interval indel
        builder.source("test9");
        start = 1038;
        stop = 1041;
        builder.start(start).stop(stop).alleles(CollectionUtil.makeList(AAAARef, AAA));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{liftOver, twoIntervalChainReference, builder.make(), null});

        // straddling interval indel
        builder.source("test10");
        start = 1045;
        stop = 1048;
        builder.start(start).stop(stop).alleles(CollectionUtil.makeList(AAAARef, AAA));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{liftOver, twoIntervalChainReference, builder.make(), null});

        // vanishing snp
        builder.source("test11");
        start = 1045;
        stop = 1045;
        builder.start(start).stop(stop).alleles(CollectionUtil.makeList(ARef, T));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{liftOver, twoIntervalChainReference, builder.make(), null});

        //  after second interval indel
        builder.source("test12");
        start = 1046;
        stop = 1049;
        offset = 0;
        builder.start(start).stop(stop).alleles(CollectionUtil.makeList(AAAARef, AAA));
        result_builder.start(start + offset).stop(stop + offset).alleles(CollectionUtil.makeList(AAAARef, AAA));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{liftOver, twoIntervalChainReference, builder.make(), result_builder.make()});

        // near start of second interval snp
        builder.source("test13");
        start = 1046;
        builder.start(start).stop(start).alleles(CollectionUtil.makeList(ARef, T));
        result_builder.start(start + offset).stop(start + offset).alleles(CollectionUtil.makeList(ARef, T));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{liftOver, twoIntervalChainReference, builder.make(), result_builder.make()});

        return tests.iterator();
    }

    @Test(dataProvider = "indelNoFlipData")
    public void testLiftOverSimpleIndels(final LiftOver liftOver, final ReferenceSequence reference, final VariantContext source, final VariantContext result) {

        final Interval target = liftOver.liftOver(new Interval(source.getContig(), source.getStart(), source.getEnd()), .95);

        VariantContextBuilder vcb = LiftoverUtils.liftSimpleVariantContext(source, target);
        VcfTestUtils.assertEquals(vcb == null ? null : vcb.make(), result);
    }

}
