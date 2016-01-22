package picard.vcf;

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

/**
 * Test class for LiftoverVcf.
 *
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

    @Test
    public void testDoNotFixReverseComplementedIndels() {
        final File liftOutputFile = new File(OUTPUT_DATA_PATH, "lift-delete-me.vcf");
        final File rejectOutputFile = new File(OUTPUT_DATA_PATH, "reject-delete-me.vcf");
        final File input = new File(TEST_DATA_PATH, "testLiftover.vcf");

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
        for (final VariantContext inputContext : liftReader) {
            Assert.fail("there should be no passing indels in the liftover");
        }
        final VCFFileReader rejectReader = new VCFFileReader(rejectOutputFile, false);
        int counter = 0;
        for (final VariantContext inputContext : rejectReader) {
            counter++;
        }
        Assert.assertEquals(counter, 2, "the wrong number of rejected indels faile the liftover");
    }

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
}
