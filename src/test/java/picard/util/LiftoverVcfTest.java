package picard.util;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFConstants;

import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.vcf.LiftoverVcf;
import picard.vcf.VcfTestUtils;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

/**
 * Test class for LiftoverVcf.
 * <p>
 * Created by ebanks on 8/11/15.
 */
public class LiftoverVcfTest extends CommandLineProgramTest {

    private static final File TEST_DATA_PATH = new File("testdata/picard/vcf/LiftOver");
    private static final File CHAIN_FILE = new File(TEST_DATA_PATH, "test.over.chain");
    private static final File POSITIVE_CHAIN_FILE = new File(TEST_DATA_PATH, "test.positive.over.chain");
    private static final File TWO_INTERVAL_CHAIN_FILE = new File(TEST_DATA_PATH, "test.two.block.over.chain");

    private static final File CHAIN_FILE_WITH_BAD_CONTIG = new File(TEST_DATA_PATH, "test.over.badContig.chain");
    private static final File REFERENCE_FILE = new File(TEST_DATA_PATH, "dummy.reference.fasta");
    private static final File TWO_INTERVALS_REFERENCE_FILE = new File(TEST_DATA_PATH, "dummy.two.block.reference.fasta");

    private static final File OUTPUT_DATA_PATH = IOUtil.createTempDir("LiftoverVcfsTest", null);

    private final int CHAIN_SIZE = 540; // the length of the single chain in CHAIN_FILE

    // comment can help find positions in string. not inline due to IDE shenanigans
    //                                       123456789 123456789 123456789 123
    private static final String refString = "CAAAAAAAAAACGTACGTACTCTCTCTCTACGT";
    private static final ReferenceSequence REFERENCE = new ReferenceSequence("chr1", 0, refString.getBytes());

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
                {"testLiftoverBiallelicIndels.vcf", 5, 0},
                {"testLiftoverMultiallelicIndels.vcf", 2, 0},
                {"testLiftoverFailingVariants.vcf", 3, 0},
                {"testLiftoverMixedVariants.vcf", 4, 0},
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
                "RECOVER_SWAPPED_REF_ALT=true",
                "CREATE_INDEX=false"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final VCFFileReader liftReader = new VCFFileReader(liftOutputFile, false);
        Assert.assertEquals(liftReader.iterator().stream().count(), expectedPassing, "The wrong number of variants were lifted over.");

        final VCFFileReader rejectReader = new VCFFileReader(rejectOutputFile, false);
        Assert.assertEquals(rejectReader.iterator().stream().count(), expectedFailing, "The wrong number of variants were rejected.");
    }

    @Test
    public void testChangingInfoFields() {
        final String filename = "testLiftoverFailingVariants.vcf";
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
                "RECOVER_SWAPPED_REF_ALT=true",
                "CREATE_INDEX=false"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final VCFFileReader liftReader = new VCFFileReader(liftOutputFile, false);
        final CloseableIterator<VariantContext> iterator = liftReader.iterator();
        while (iterator.hasNext()) {
            final VariantContext vc = iterator.next();
            Assert.assertTrue(vc.hasAttribute("AF"));

            Assert.assertEquals(vc.hasAttribute(LiftoverUtils.SWAPPED_ALLELES), vc.hasAttribute("FLIPPED_AF"), vc.toString());
            Assert.assertEquals(!vc.hasAttribute(LiftoverUtils.SWAPPED_ALLELES), vc.hasAttribute("MAX_AF") && vc.getAttribute("MAX_AF") != null, vc.toString());

            if (vc.hasAttribute(LiftoverUtils.SWAPPED_ALLELES)) {
                final double af = vc.getAttributeAsDouble("AF", -1.0);
                final double flippedAf = vc.getAttributeAsDouble("FLIPPED_AF", -2.0);

                Assert.assertTrue(TestNGUtil.compareDoubleWithAccuracy(af, flippedAf, 0.01),
                        "Error while comparing AF. expected " + flippedAf + " but found " + af);

            }
        }
    }

    @Test
    public void testChangingTagArguments() {
        final String filename = "testLiftoverFailingVariants.vcf";
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
                "TAGS_TO_REVERSE=null",
                "TAGS_TO_DROP=null",
                "CREATE_INDEX=false"
        };

        // we don't actually care what the results are here -- we just want to make sure that it doesn't fail
    }

    @Test
    public void testReverseComplementFailureDoesNotErrorOut() {
        final VariantContextBuilder builder = new VariantContextBuilder().source("test").loc("chr1", 1, 4);
        final Allele originalRef = Allele.create("CCCC", true);
        final Allele originalAlt = Allele.create("C", false);
        builder.alleles(Arrays.asList(originalRef, originalAlt));

        final Interval interval = new Interval("chr1", 1, 4, true, "test ");

        final String reference = "ATGATGATGA";
        final ReferenceSequence refSeq = new ReferenceSequence("chr1", 10, reference.getBytes());

        // we don't actually care what the results are here -- we just want to make sure that it doesn't fail
        final VariantContextBuilder result = LiftoverUtils.reverseComplementVariantContext(builder.make(), interval, refSeq);
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

        Assert.assertEquals(runPicardCommandLine(args), 0);

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

    @DataProvider(name = "dataTestSort")
    public Object[][] dataTestVcfSorted() {
        return new Object[][]{
                {false},
                {true}
        };
    }

    @Test(dataProvider = "dataTestSort")
    public void testVcfSorted(final boolean disableSort) {
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
                "CREATE_INDEX=" + (!disableSort),
                "DISABLE_SORT=" + disableSort
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);
        //check input/header
        try (final VCFFileReader inputReader = new VCFFileReader(input, false)) {
                final VCFHeader header = inputReader.getFileHeader();
                Assert.assertNull(header.getOtherHeaderLine("reference"));
                }
        //try to open with / without index
        try (final VCFFileReader liftReader = new VCFFileReader(liftOutputFile, !disableSort)) {
                final VCFHeader header = liftReader.getFileHeader();
                Assert.assertNotNull(header.getOtherHeaderLine("reference"));
                try (final CloseableIterator<VariantContext> iter = liftReader.iterator()) {
                Assert.assertEquals(iter.stream().count(), 5L);
                }
            }
        }
    

    
    
    @DataProvider
    Iterator<Object[]> testWriteVcfData() {

        List<Object[]> tests = new ArrayList<>();
        {
            final File input = new File(TEST_DATA_PATH, "testLiftoverMismatchingSnps.vcf");
            final File expectedVcf = new File(TEST_DATA_PATH, "vcfWithFlippedAlleles.lift.vcf");
            final File expectedRejectVcf = new File(TEST_DATA_PATH, "vcfWithFlippedAlleles.reject.vcf");

            tests.add(new Object[]{input, expectedVcf, expectedRejectVcf, TWO_INTERVALS_REFERENCE_FILE, TWO_INTERVAL_CHAIN_FILE});
        }

        {
            final File input = new File(TEST_DATA_PATH, "testLiftoverMismatchingSnps.vcf");
            final File expectedVcf = new File(TEST_DATA_PATH, "vcfWithFlippedAllelesNegativeChain.lift.vcf");
            final File expectedRejectVcf = new File(TEST_DATA_PATH, "vcfWithFlippedAllelesNegativeChain.reject.vcf");

            tests.add(new Object[]{input, expectedVcf, expectedRejectVcf, REFERENCE_FILE, CHAIN_FILE});
        }
        {
            final File input = new File(TEST_DATA_PATH, "testLiftoverMixedVariants.vcf");
            final File expectedVcf = new File(TEST_DATA_PATH, "vcfWithMixed.lift.vcf");
            final File expectedRejectVcf = new File(TEST_DATA_PATH, "vcfWithMixed.reject.vcf");

            tests.add(new Object[]{input, expectedVcf, expectedRejectVcf, REFERENCE_FILE, CHAIN_FILE});
        }

        return tests.iterator();
    }

    @Test(dataProvider = "testWriteVcfData")
    public void testWriteVcfWithFlippedAlleles(
            final File input,
            final File expectedVcf,
            final File expectedRejectVcf,
            final File reference,
            final File liftoverChain) throws IOException {

        final File liftOutputFile = new File(OUTPUT_DATA_PATH, "lift-delete-me.vcf");
        final File rejectOutputFile = new File(OUTPUT_DATA_PATH, "reject-delete-me.vcf");

        liftOutputFile.deleteOnExit();
        rejectOutputFile.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + liftOutputFile.getAbsolutePath(),
                "REJECT=" + rejectOutputFile.getAbsolutePath(),
                "CHAIN=" + liftoverChain.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
                "RECOVER_SWAPPED_REF_ALT=true",
                "CREATE_INDEX=false"
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        try (final VCFFileReader liftReader = new VCFFileReader(liftOutputFile, false)) {
            Assert.assertTrue(liftReader.getFileHeader().hasInfoLine(LiftoverUtils.SWAPPED_ALLELES));
            for (final VariantContext vc : liftReader) {
                Assert.assertFalse(vc.hasAttribute(LiftoverVcf.ORIGINAL_CONTIG));
                Assert.assertFalse(vc.hasAttribute(LiftoverVcf.ORIGINAL_START));

            }
        }

        VcfTestUtils.assertVcfFilesAreEqual(liftOutputFile, expectedVcf);
        VcfTestUtils.assertVcfFilesAreEqual(rejectOutputFile, expectedRejectVcf);
    }


    @DataProvider(name = "indelFlipData")
    public Iterator<Object[]> indelFlipData() {

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

        final Allele spanningDeletion = Allele.create(Allele.SPAN_DEL_STRING, false);

        final List<Object[]> tests = new ArrayList<>();

        final VariantContextBuilder builder = new VariantContextBuilder().source("test1").chr("chr1");
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder("test1");
        final GenotypeBuilder resultGenotypeBuilder = new GenotypeBuilder("test1");
        final VariantContextBuilder result_builder = new VariantContextBuilder().source("test1").chr("chr1");

        // simple deletion
        // TTT*/T -> AAA*/A turns into left-aligned CAA*/C at position 1
        int start = CHAIN_SIZE - 3;
        int stop = start + 2;
        builder.start(start).stop(stop).alleles(CollectionUtil.makeList(RefTTT, T));
        result_builder.start(1).stop(3).alleles(CollectionUtil.makeList(RefCAA, C)).attribute(LiftoverUtils.REV_COMPED_ALLELES, true);

        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());

        tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});

        //simple insertion
        // T*/TTT -> A*/AAA -> turns into left-aligned C*/CAA at position 1
        stop = start;
        builder.source("test2").alleles(CollectionUtil.makeList(RefT, TTT)).stop(stop);
        result_builder.alleles(CollectionUtil.makeList(RefC, CAA)).start(1).stop(1);
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());

        tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});

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
        tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});

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
        tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});

        // just outside of chain & contig, testing that we do not read into negative indices
        // or reference
        start = stop = CHAIN_SIZE;
        builder.source("test5").stop(stop).start(start).alleles(CollectionUtil.makeList(RefG, GTT));
        result_builder.start(1).stop(1).alleles(CollectionUtil.makeList(RefC, AAC));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});

        // outside of chain
        start = stop = CHAIN_SIZE + 1;
        builder.source("test6").stop(stop).start(start).alleles(CollectionUtil.makeList(RefA, ACGT));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{builder.make(), REFERENCE, null});

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
        tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});

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
        tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});

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
        tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});

        // insertion at end of section
        // a test that converts the initial C to a AC which requires
        // code in LiftOver::extendOneBase(., false)
        //
        // G*(.)/GA(.) -> .(C)/.T(C) but that's not legal, so
        // -> C/TC
        start = CHAIN_SIZE;
        stop = CHAIN_SIZE;

        builder.source("test10").start(start).stop(stop).alleles(CollectionUtil.makeList(RefG, GA));
        result_builder.start(1).stop(1).alleles(CollectionUtil.makeList(RefC, TC));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});

        // insertion at end of section
        // a test that converts the initial C to a AC which requires
        // code in LiftOver::extendOneBase(., false)
        //
        // G*(.)/GG(.) -> .(C)/.C(C) but that's not legal, so
        // -> C/CC
        start = CHAIN_SIZE;
        stop = CHAIN_SIZE;

        builder.source("test11").start(start).stop(stop).alleles(CollectionUtil.makeList(RefG, Allele.create("GG")));
        result_builder.start(1).stop(1).alleles(CollectionUtil.makeList(RefC, Allele.create("CC")));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});

        // insertion at end of section
        // a test that converts the initial C to a AC which requires
        // code in LiftOver::extendOneBase(., false)
        //
        // improperly aligned
        // TTTT(G)->TTTTG(G) -> C/CC
        start = CHAIN_SIZE - 4;
        stop = CHAIN_SIZE - 1;

        builder.source("test12").start(start).stop(stop).alleles(CollectionUtil.makeList(Allele.create("TTTT", true), Allele.create("TTTTG")));
        result_builder.start(1).stop(1).alleles(CollectionUtil.makeList(RefC, Allele.create("CC")));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});

        // insertion at end of section, with a spanning deletion
        // a test that converts the initial C to a AC which requires
        // code in LiftOver::extendOneBase(., false)
        //
        // improperly aligned with spanning deletion
        // TTTT(G)->TTTTG(G) -> C/CC

        start = CHAIN_SIZE - 4;
        stop = CHAIN_SIZE - 1;

        // Note that the spanning deletion prevents the indel from being left aligned.
        builder.source("test13").start(start).stop(stop).alleles(CollectionUtil.makeList(Allele.create("TTTT", true), Allele.create("TTTTG"), spanningDeletion));
        result_builder.start(1).stop(1).alleles(CollectionUtil.makeList(RefC, Allele.create("CC"), spanningDeletion));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});

        // non-simple deletion with spanning deletion
        //  ACGT(T)*/A(T) -> AACG(T)*/A(T)
        //       position of 13 ^ in result
        start = CHAIN_SIZE - 13;
        stop = start + 3;

        builder.source("test14").start(start).stop(stop).alleles(CollectionUtil.makeList(RefACGT, A, spanningDeletion));
        result_builder.start(CHAIN_SIZE - stop).stop(CHAIN_SIZE - start).alleles(CollectionUtil.makeList(RefAACG, A, spanningDeletion));
        genotypeBuilder.alleles(builder.getAlleles());
        resultGenotypeBuilder.alleles(result_builder.getAlleles());
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});
        
        builder.noGenotypes();
        result_builder.noGenotypes();

        return tests.iterator();
    }

    @Test(dataProvider = "indelFlipData")
    public void testFlipIndel(final VariantContext source, final ReferenceSequence reference, final VariantContext result) {

        final LiftOver liftOver = new LiftOver(CHAIN_FILE);
        final Interval originalLocus = new Interval(source.getContig(), source.getStart(), source.getEnd());
        final Interval target = liftOver.liftOver(originalLocus);
        if (target != null && !target.isNegativeStrand()) {
            throw new RuntimeException("not reversed");
        }

        final VariantContext flipped = LiftoverUtils.liftVariant(source, target, reference, false, false);

        VcfTestUtils.assertEquals(flipped, result);
    }

    @DataProvider(name = "indelFlipDataWithOriginalAllele")
    public Iterator<Object[]> indelFlipDataWithOriginalAllele() {
        final List<Object[]> tests = new ArrayList<>();
        indelFlipData().forEachRemaining(x -> {
            if (x[2] == null) {
                tests.add(x);
            } else {
                VariantContext source = (VariantContext) x[0];
                VariantContextBuilder result_builder = new VariantContextBuilder((VariantContext) x[2]);
                result_builder.attribute(LiftoverVcf.ORIGINAL_ALLELES, LiftoverUtils.allelesToStringList(source.getAlleles()));
                tests.add(new Object[]{x[0], x[1], result_builder.make()});
            }
        });

        return tests.iterator();
    }

    @Test(dataProvider = "indelFlipDataWithOriginalAllele")
    public void testFlipIndelWithOriginalAlleles(final VariantContext source, final ReferenceSequence reference, final VariantContext result) {

        final LiftOver liftOver = new LiftOver(CHAIN_FILE);
        final Interval originalLocus = new Interval(source.getContig(), source.getStart(), source.getEnd());
        final Interval target = liftOver.liftOver(originalLocus);
        if (target != null && !target.isNegativeStrand()) {
            throw new RuntimeException("not reversed");
        }

        final VariantContext flipped = LiftoverUtils.liftVariant(source, target, reference, false, true);

        VcfTestUtils.assertEquals(flipped, result);
    }

    @DataProvider(name = "snpWithChangedRef")
    public Iterator<Object[]> snpWithChangedRef() {

        final ReferenceSequence reference = new ReferenceSequence("chr1", 0,
                refString.getBytes());
        //       123456789 123456789 123456789 123

        final Allele RefT = Allele.create("T", true);
        final Allele RefA = Allele.create("A", true);
        final Allele RefC = Allele.create("C", true);
        final Allele RefG = Allele.create("G", true);

        final Allele A = Allele.create("A", false);
        final Allele T = Allele.create("T", false);
        final Allele C = Allele.create("C", false);
        final Allele G = Allele.create("G", false);

        final List<Object[]> tests = new ArrayList<>();

        final int CHAIN_SIZE = 540; // the length of the single chain in POSITIVE_CHAIN_FILE

        final VariantContextBuilder builder = new VariantContextBuilder().source("test1").chr("chr1");
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder("test1");
        final GenotypeBuilder resultGenotypeBuilder = new GenotypeBuilder("test1");
        final VariantContextBuilder result_builder = new VariantContextBuilder().source("test1").chr("chr1").attribute(LiftoverUtils.SWAPPED_ALLELES, true);

        // simple snp
        int start = 12;
        int stop = 12;
        builder.start(start).stop(stop).alleles(CollectionUtil.makeList(RefA, C));
        result_builder.start(start).stop(stop).alleles(CollectionUtil.makeList(A, RefC));

        genotypeBuilder.alleles(builder.getAlleles()).AD(new int[]{5, 4}).PL(new int[]{40, 0, 60});
        resultGenotypeBuilder.alleles(result_builder.getAlleles()).AD(new int[]{4, 5}).PL(new int[]{60, 0, 40});

        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());

        tests.add(new Object[]{builder.make(), reference, result_builder.make()});

        builder.start(start).stop(stop).alleles(CollectionUtil.makeList(RefA, C));
        result_builder.start(start).stop(stop).alleles(CollectionUtil.makeList(RefC, A));

        genotypeBuilder.alleles(Arrays.asList(C, C)).AD(new int[]{0, 10}).PL(new int[]{400, 40, 0});
        resultGenotypeBuilder.alleles(Arrays.asList(RefC, RefC)).AD(new int[]{10, 0}).PL(new int[]{0, 40, 400});

        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());

        tests.add(new Object[]{builder.make(), reference, result_builder.make()});
        return tests.iterator();
    }

    @Test(dataProvider = "snpWithChangedRef")
    public void snpWithChangedRef(final VariantContext source, final ReferenceSequence reference, final VariantContext result) {

        final LiftOver liftOver = new LiftOver(POSITIVE_CHAIN_FILE);
        final Interval originalLocus = new Interval(source.getContig(), source.getStart(), source.getEnd());
        final Interval target = liftOver.liftOver(originalLocus);

        final VariantContext flipped = LiftoverUtils.swapRefAlt(source, LiftoverUtils.DEFAULT_TAGS_TO_REVERSE, LiftoverUtils.DEFAULT_TAGS_TO_DROP);

        VcfTestUtils.assertEquals(flipped, result);
    }

    @DataProvider(name = "leftAlignAllelesData")
    public Iterator<Object[]> leftAlignAllelesData() {

        final Allele RefG = Allele.create("G", true);
        final Allele A = Allele.create("A", false);

        final Allele RefA = Allele.create("A", true);
        final Allele RefC = Allele.create("C", true);
        final Allele RefAA = Allele.create("AA", true);
        final Allele RefAC = Allele.create("AC", true);
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

        tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});

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
        tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});

        for (start = 1; start <= REFERENCE.getBases().length; start++) {
            builder.source("test2-" + start);
            builder.start(start).stop(start);
            builder.alleles(CollectionUtil.makeList(
                    Allele.create(REFERENCE.getBaseString().substring(start - 1, start), true),
                    REFERENCE.getBaseString().charAt(start - 1) == 'A' ? T : A));

            genotypeBuilder.alleles(builder.getAlleles());
            resultGenotypeBuilder.alleles(result_builder.getAlleles());
            builder.genotypes(genotypeBuilder.make());
            result_builder.genotypes(resultGenotypeBuilder.make());
            tests.add(new Object[]{builder.make(), REFERENCE, builder.make()});
        }

        // AA/A in initial polyA repeat -> CA/C at the beginning
        result_builder.start(1).stop(2).alleles(CollectionUtil.makeList(RefCA, C));
        for (start = 2; start <= 10; start++) {
            builder.source("test3-" + start);
            builder.start(start).stop(start + 1);
            builder.alleles(CollectionUtil.makeList(
                    RefAA, A));

            genotypeBuilder.alleles(builder.getAlleles());
            resultGenotypeBuilder.alleles(result_builder.getAlleles());
            builder.genotypes(genotypeBuilder.make());
            result_builder.genotypes(resultGenotypeBuilder.make());
            tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});
        }

        // A/AA in initial polyA repeat -> C/CA at the beginning
        result_builder.start(1).stop(1).alleles(CollectionUtil.makeList(RefC, CA));
        for (start = 2; start <= 11; start++) {
            builder.source("test4-" + start);
            builder.start(start).stop(start);
            builder.alleles(CollectionUtil.makeList(
                    RefA, AA));

            genotypeBuilder.alleles(builder.getAlleles());
            resultGenotypeBuilder.alleles(result_builder.getAlleles());
            builder.genotypes(genotypeBuilder.make());
            result_builder.genotypes(resultGenotypeBuilder.make());
            tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});
        }

        //CT/CTCT -> A/ACT in CT repeat region
        result_builder.start(19).stop(19).alleles(CollectionUtil.makeList(RefA, ACT));
        for (start = 20; start <= 27; start += 2) {
            builder.source("test5-" + start);
            builder.start(start).stop(start + 1);
            builder.alleles(CollectionUtil.makeList(
                    RefCT, CTCT));

            genotypeBuilder.alleles(builder.getAlleles());
            resultGenotypeBuilder.alleles(result_builder.getAlleles());
            builder.genotypes(genotypeBuilder.make());
            result_builder.genotypes(resultGenotypeBuilder.make());
            tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});
        }

        //TC/TCTC -> A/ACT in CT repeat region
        for (start = 21; start < 29; start += 2) {
            builder.source("test6-" + start);
            builder.start(start).stop(start + 1);
            builder.alleles(CollectionUtil.makeList(
                    RefTC, TCTC));

            genotypeBuilder.alleles(builder.getAlleles());
            resultGenotypeBuilder.alleles(result_builder.getAlleles());
            builder.genotypes(genotypeBuilder.make());
            result_builder.genotypes(resultGenotypeBuilder.make());
            tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});
        }

        //CTCT/CT -> ACT/A in CT repeat region
        result_builder.start(19).stop(21).alleles(CollectionUtil.makeList(RefACT, A));
        for (start = 20; start <= 27; start += 2) {
            builder.source("test7-" + start);
            builder.start(start).stop(start + 3);
            builder.alleles(CollectionUtil.makeList(
                    RefCTCT, CT));

            genotypeBuilder.alleles(builder.getAlleles());
            resultGenotypeBuilder.alleles(result_builder.getAlleles());
            builder.genotypes(genotypeBuilder.make());
            result_builder.genotypes(resultGenotypeBuilder.make());
            tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});
        }

        //TCTC/TC-> ACT/A in CT repeat region
        for (start = 21; start < 27; start += 2) {
            builder.start(start).stop(start + 3);
            builder.alleles(CollectionUtil.makeList(
                    RefTCTC, TC));
            builder.source("test8-" + start);

            genotypeBuilder.alleles(builder.getAlleles());
            resultGenotypeBuilder.alleles(result_builder.getAlleles());
            builder.genotypes(genotypeBuilder.make());
            result_builder.genotypes(resultGenotypeBuilder.make());
            tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});
        }

        // for ease of reading, here's the reference sequence
        // "CAAAAAAAAAACGTACGTACTCTCTCTCTACGT"
        //  123456789 123456789 123456789 123

        result_builder.alleles("AACGT", "A").start(10).stop(14);
        for (start = 10; start < 17; start++) {
            for (stop = start + 4; stop < 20; stop++) {
                builder.source("test9-" + start + "-" + stop);
                builder.alleles(
                        // -1 here due to reference string being 0-based.
                        REFERENCE.getBaseString().substring(start - 1, stop + 1 - 1),
                        REFERENCE.getBaseString().substring(start - 1, stop - 3 - 1)).start(start).stop(stop);
                genotypeBuilder.alleles(builder.getAlleles());
                resultGenotypeBuilder.alleles(result_builder.getAlleles());
                builder.genotypes(genotypeBuilder.make());
                result_builder.genotypes(resultGenotypeBuilder.make());
                tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});
            }
        }

        // test vc with genotypes:

        result_builder.start(10).stop(14).alleles("AACGT", "A", "AACG", "AACGTACGT");
        builder.start(9).stop(18).alleles("AAACGTACGT", "AAACGT", "AAACGACGT", "AAACGTACGTACGT");
        final Collection<Genotype> genotypes = new ArrayList<>();
        final Collection<Genotype> results_genotypes = new ArrayList<>();
        final GenotypeBuilder gBuilder = new GenotypeBuilder();

        for (int sample = 1; sample < 4; sample++) {
            builder.source("test10-" + start);
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
        tests.add(new Object[]{builder.make(), REFERENCE, result_builder.make()});

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

    @DataProvider(name = "noCallAndSymbolicData")
    public Iterator<Object[]> noCallAndSymbolicData() {

        final VariantContextBuilder builder = new VariantContextBuilder().source("test1").chr("chr1").attribute("TEST_ATTRIBUTE",1).attribute("TEST_ATTRIBUTE2","Hi").attribute("TEST_ATTRIBUTE3",new double[] {1.2,2.1});
        final VariantContextBuilder result_builder = new VariantContextBuilder().source("test1").chr("chr1").attribute("TEST_ATTRIBUTE",1).attribute("TEST_ATTRIBUTE2","Hi").attribute("TEST_ATTRIBUTE3",new double[] {1.2,2.1});
        result_builder.attribute(LiftoverVcf.ORIGINAL_CONTIG, "chr1");
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder("test1");
        final GenotypeBuilder resultGenotypeBuilder = new GenotypeBuilder("test1");
        final List<Object[]> tests = new ArrayList<>();

        final Allele CRef = Allele.create("C", true);
        final Allele GRef = Allele.create("G", true);
        final Allele T = Allele.create("T", false);
        final Allele A = Allele.create("A", false);
        final Allele DEL = Allele.create("*", false);

        final LiftOver liftOver = new LiftOver(TWO_INTERVAL_CHAIN_FILE);
        final LiftOver liftOverRC = new LiftOver(CHAIN_FILE);

        builder.source("test1");
        int start = 10;
        builder.start(start).stop(start).alleles(CollectionUtil.makeList(CRef, T));
        result_builder.start(start).stop(start).alleles(CollectionUtil.makeList(CRef, T));
        genotypeBuilder.alleles(CollectionUtil.makeList(Allele.create("."), Allele.create(".")));
        resultGenotypeBuilder.alleles(CollectionUtil.makeList(Allele.create("."), Allele.create(".")));
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());
        result_builder.attribute(LiftoverVcf.ORIGINAL_START, start);

        tests.add(new Object[]{liftOver, builder.make(), result_builder.make(), false});

        builder.source("test2");
        builder.start(start).stop(start).alleles(CollectionUtil.makeList(CRef, T, DEL));
        result_builder.start(start).stop(start).alleles(CollectionUtil.makeList(CRef, T, DEL));
        result_builder.attribute(LiftoverVcf.ORIGINAL_START, start);
        genotypeBuilder.alleles(CollectionUtil.makeList(T, DEL));
        resultGenotypeBuilder.alleles(CollectionUtil.makeList(T, DEL));
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());

        tests.add(new Object[]{liftOver, builder.make(), result_builder.make(), false});

        //reverse complement
        builder.source("test3");
        int offset = 3;
        start = CHAIN_SIZE - offset;
        int liftedStart = 1 + offset;
        builder.start(start).stop(start).alleles(CollectionUtil.makeList(CRef, T, DEL));
        result_builder.start(liftedStart).stop(liftedStart).alleles(CollectionUtil.makeList(GRef, A, DEL));
        result_builder.attribute(LiftoverVcf.ORIGINAL_START, start);
        result_builder.attribute(LiftoverVcf.ORIGINAL_ALLELES, LiftoverUtils.allelesToStringList(builder.getAlleles()));
        result_builder.attribute(LiftoverUtils.REV_COMPED_ALLELES, true);


        genotypeBuilder.alleles(CollectionUtil.makeList(T, DEL));
        resultGenotypeBuilder.alleles(CollectionUtil.makeList(A, DEL));
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());

        tests.add(new Object[]{liftOverRC, builder.make(), result_builder.make(), true});

        builder.source("test4");
        offset = 4;
        start = CHAIN_SIZE - offset;
        liftedStart = 1 + offset;
        builder.start(start).stop(start).alleles(CollectionUtil.makeList(CRef, T));
        result_builder.start(liftedStart).stop(liftedStart).alleles(CollectionUtil.makeList(GRef, A));
        result_builder.attribute(LiftoverVcf.ORIGINAL_START, start);
        result_builder.attribute(LiftoverVcf.ORIGINAL_ALLELES, LiftoverUtils.allelesToStringList(builder.getAlleles()));
        genotypeBuilder.alleles(CollectionUtil.makeList(T, Allele.NO_CALL));
        resultGenotypeBuilder.alleles(CollectionUtil.makeList(A, Allele.NO_CALL));
        builder.genotypes(genotypeBuilder.make());
        result_builder.genotypes(resultGenotypeBuilder.make());

        tests.add(new Object[]{liftOverRC, builder.make(), result_builder.make(), true});

        return tests.iterator();
    }

    @Test(dataProvider = "noCallAndSymbolicData")
    public void testLiftOverNoCallAndSymbolic(final LiftOver liftOver, final VariantContext source, final VariantContext result, final boolean expectReversed) {

        final Interval target = liftOver.liftOver(new Interval(source.getContig(), source.getStart(), source.getEnd()), .95);

        Assert.assertEquals(target.isNegativeStrand(), expectReversed);

        VariantContext vc = LiftoverUtils.liftVariant(source, target, REFERENCE, true, true);
        VcfTestUtils.assertEquals(vc, result);

        Assert.assertEquals(vc.getAttribute(LiftoverVcf.ORIGINAL_CONTIG), source.getContig());
        Assert.assertEquals(vc.getAttribute(LiftoverVcf.ORIGINAL_START), source.getStart());
        Assert.assertTrue(source.getAlleles().equals(result.getAlleles()) != result.hasAttribute(LiftoverVcf.ORIGINAL_ALLELES));
        if (!source.getAlleles().equals(result.getAlleles())) {
            List<String> resultAlleles = new ArrayList<>();
            source.getAlleles().forEach(a -> resultAlleles.add(a.getDisplayString()));
            Assert.assertEquals(resultAlleles, result.getAttributeAsStringList(LiftoverVcf.ORIGINAL_ALLELES, null));
        }
        result.getAttributes().entrySet().stream().filter(e->e.getKey().matches("TEST_ATTRIBUTE.*")).forEach(e -> {

            Assert.assertTrue(vc.hasAttribute(e.getKey()));
            Assert.assertEquals(vc.getAttribute(e.getKey()), e.getValue());
        });
    }

    @Test()
    public void testNoDictionary() throws IOException {
        final Path liftOutput = Files.createTempFile("tmpoutput", ".vcf");
        liftOutput.toFile().deleteOnExit();
        final Path rejectOutput = Files.createTempFile("tmpreject", ".vcf");
        rejectOutput.toFile().deleteOnExit();
        final Path input = TEST_DATA_PATH.toPath().resolve("testLiftoverBiallelicIndels.vcf");
        final Path tmpCopyDir = Files.createTempDirectory("copy");
        tmpCopyDir.toFile().deleteOnExit();
        final Path referenceCopy = tmpCopyDir.resolve("refCopy.fasta");
        referenceCopy.toFile().deleteOnExit();
        Files.copy(REFERENCE_FILE.toPath(), referenceCopy);
        final String[] args = new String[]{
                "INPUT=" + input,
                "OUTPUT=" + liftOutput,
                "REJECT=" + rejectOutput,
                "CHAIN=" + CHAIN_FILE,
                "REFERENCE_SEQUENCE=" + referenceCopy,
        };
        Assert.assertEquals(runPicardCommandLine(args), 1);
    }

    @Test(expectedExceptions = SAMException.class, expectedExceptionsMessageRegExp = "File exists but is not readable:.*")
    public void testUnreadableDictionary() throws IOException {
        final Path liftOutput = Files.createTempFile("tmpouput", ".vcf");
        liftOutput.toFile().deleteOnExit();
        final Path rejectOutput = Files.createTempFile("tmpreject", ".vcf");
        rejectOutput.toFile().deleteOnExit();
        final Path input = TEST_DATA_PATH.toPath().resolve("testLiftoverBiallelicIndels.vcf");
        final Path tmpCopyDir = Files.createTempDirectory("copy");
        tmpCopyDir.toFile().deleteOnExit();
        final Path referenceCopy = tmpCopyDir.resolve("refCopy.fasta");
        referenceCopy.toFile().deleteOnExit();
        Files.copy(REFERENCE_FILE.toPath(), referenceCopy);
        final Path dictCopy = referenceCopy.resolveSibling(referenceCopy.toFile().getName().replaceAll("fasta$", "dict"));
        dictCopy.toFile().deleteOnExit();
        final Path dictionary = REFERENCE_FILE.toPath().resolveSibling(REFERENCE_FILE.getName().replaceAll("fasta$", "dict"));
        Files.copy(dictionary, dictCopy);
        dictCopy.toFile().setReadable(false);
        final String[] args = new String[]{
                "INPUT=" + input,
                "OUTPUT=" + liftOutput,
                "REJECT=" + rejectOutput,
                "CHAIN=" + CHAIN_FILE,
                "REFERENCE_SEQUENCE=" + referenceCopy,
        };
        runPicardCommandLine(args);
    }
}
