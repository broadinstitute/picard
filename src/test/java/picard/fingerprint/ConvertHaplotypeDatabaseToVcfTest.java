package picard.fingerprint;

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class ConvertHaplotypeDatabaseToVcfTest extends CommandLineProgramTest {
    private final File TEST_DATA_DIR = new File("testdata/picard/fingerprint/");
    private final File TEST_MAP =
            new File(TEST_DATA_DIR, "haplotypeMap_small.txt");

    private final File TEST_FASTA =
            new File(TEST_DATA_DIR, "reference.fasta");

    @Override
    public String getCommandLineProgramName() {
        return ConvertHaplotypeDatabaseToVcf.class.getSimpleName();
    }

    @Test
    public void testConvertHaplotypeDatabaseToVcf() throws IOException {
        final File outfile = getTempOutputFile("testConvertHaplotypeDatabase", ".vcf");
        final String[] args = new String[]{
                "I=" + TEST_MAP.getAbsolutePath(),
                "R=" + TEST_FASTA.getAbsolutePath(),
                "O=" + outfile.getAbsolutePath()
        };

        runPicardCommandLine(args);

        final HaplotypeMap hapDBFromTxt = new HaplotypeMap(TEST_MAP);
        final HaplotypeMap hapDBFromVcf = new HaplotypeMap(outfile);

        assertHaplotypeMapsAreEquivalent(hapDBFromTxt, hapDBFromVcf);
    }

    private void assertHaplotypeMapsAreEquivalent(final HaplotypeMap hapMap1, final HaplotypeMap hapMap2) {
        final List<HaplotypeBlock> hapBlocks1 = new ArrayList<>(hapMap1.getHaplotypes());
        final List<HaplotypeBlock> hapBlocks2 = new ArrayList<>(hapMap2.getHaplotypes());

        Assert.assertEquals(hapBlocks1.size(), hapBlocks2.size());

        Collections.sort(hapBlocks1);
        Collections.sort(hapBlocks2);

        for (int i = 0; i < hapBlocks1.size(); i++) {
            assertHaplotypeBlocksAreEquivalent(hapBlocks1.get(i), hapBlocks2.get(i));
        }
    }

    private void assertHaplotypeBlocksAreEquivalent(final HaplotypeBlock hapBlock1, final HaplotypeBlock hapBlock2) {
        //equals method only checks for equality of domain
        Assert.assertEquals(hapBlock1, hapBlock2);

        final List<Snp> snps1 = new ArrayList<>(hapBlock1.getSnps());
        final List<Snp> snps2 = new ArrayList<>(hapBlock2.getSnps());

        Assert.assertEquals(snps1.size(), snps2.size());

        Collections.sort(snps1);
        Collections.sort(snps2);

        for (int i = 0; i < snps1.size(); i++) {
            assertSnpsAreEquivalent(snps1.get(i), snps2.get(i));
        }
    }

    private void assertSnpsAreEquivalent(final Snp snp1, final Snp snp2) {
        //equals method only checks for equality of domain
        Assert.assertEquals(snp1, snp2);

        Assert.assertEquals(snp1.getName(), snp2.getName());

        //allele order may have been swapped
        if (snp1.getAllele1() == snp2.getAllele1() && snp1.getAllele2() == snp2.getAllele2()) {
            Assert.assertEquals(snp1.getMaf(), snp2.getMaf(), 0.0005); //equal rounded to third decimal
        } else if (snp1.getAllele1() == snp2.getAllele2() && snp1.getAllele2() == snp2.getAllele1()) {
            Assert.assertEquals(snp1.getMaf(), 1 - snp2.getMaf(), 0.0005); //equal rounded to third decimal
        } else {
            throw new AssertionError("Snps " + snp1.getName() + " have different alleles");
        }
    }
}