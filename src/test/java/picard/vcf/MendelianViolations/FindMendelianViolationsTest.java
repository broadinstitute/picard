package picard.vcf.MendelianViolations;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.TestUtil;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.variant.vcf.VCFCodec;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.pedigree.Sex;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Pattern;

import static picard.vcf.MendelianViolations.MendelianViolationDetector.MendelianViolation.*;

/**
 * Created by farjoun on 5/13/16.
 */
public class FindMendelianViolationsTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/vcf");

    @Test
    public void testFindMedelianViolations() throws IOException {
        final File vcfFile = new File(TEST_DATA_DIR, "CEUTrio_plus_FAKE.vcf");

        final File vcfIndexFile = new File(TEST_DATA_DIR, "CEUTrio_plus_FAKE.vcf.idx");

        if (vcfFile.lastModified() > vcfIndexFile.lastModified()) {
            if (vcfIndexFile.exists()) {
                System.err.println("Deleting " + vcfIndexFile);
                vcfIndexFile.delete();
            }
            IndexFactory.createDynamicIndex(vcfFile, new VCFCodec()).writeBasedOnFeatureFile(vcfFile);
        }
        final File pedFile = new File(TEST_DATA_DIR, "NA12878.ped");

        final File directoryForViolations = TestUtil.getTempDirectory("MendelianViolations", "temp");
        directoryForViolations.mkdir();
        final File violationsNA12878 = new File(directoryForViolations, "NA12878.vcf");
        final File violationsFAKE = new File(directoryForViolations, "FAKE.vcf");
        violationsNA12878.deleteOnExit();
        violationsFAKE.deleteOnExit();
        directoryForViolations.deleteOnExit();

        final File resultantMetrics = File.createTempFile("MendelianViolations", "file");
        resultantMetrics.deleteOnExit();

        final FindMendelianViolations program = new FindMendelianViolations();
        program.INPUT = vcfFile;
        program.TRIOS = pedFile;
        program.OUTPUT = resultantMetrics;
        program.VCF_DIR = directoryForViolations;
        program.MIN_DP = 10;

        Assert.assertEquals(program.doWork(), 0);

        final MetricsFile<MendelianViolationMetrics, Comparable<?>> summary = new MetricsFile<>();
        summary.read(new FileReader(resultantMetrics));

        Assert.assertEquals(summary.getMetrics().size(), 2);

        MendelianViolationMetrics metrics = summary.getMetrics().get(0);

        final MendelianViolationMetrics mv = new MendelianViolationMetrics();

        mv.FAMILY_ID = "NA12878";
        mv.OFFSPRING = "NA12878";
        mv.OFFSPRING_SEX = Sex.Female;
        mv.FATHER = "NA12891";
        mv.MOTHER = "NA12892";
        mv.NUM_VARIANT_SITES =     grep(vcfFile, "GoodVariant") +         grep(vcfFile, "GoodVariantNA12878");
        mv.NUM_DIPLOID_DENOVO =    grep(vcfFile, "DiploidDenovoBOTH") +   grep(vcfFile, "DiploidDenovoNA12878");
        mv.NUM_HOMREF_HOMVAR_HOM = grep(vcfFile, "HomrefHomvarHomBOTH") + grep(vcfFile, "HomrefHomvarHomNA12878");
        mv.NUM_HOMVAR_HOMVAR_HET = grep(vcfFile, "HomvarHomvarHetBOTH") + grep(vcfFile, "HomrefHomvarHetNA12878");

        mv.NUM_HOM_HET_HOM =    grep(vcfFile, "HomHetHomBOTH") +     grep(vcfFile, "HomHetHomNA12878");
        mv.NUM_HAPLOID_DENOVO = grep(vcfFile, "HaploidDenovoBOTH") + grep(vcfFile, "HaploidDenovoNA12878");
        mv.NUM_HAPLOID_OTHER =  grep(vcfFile, "HaploidOtherBOTH") +  grep(vcfFile, "HaploidOtherNA12878");
        mv.NUM_OTHER =          grep(vcfFile, "OtherBOTH") +         grep(vcfFile, "OtherNA12878");

        mv.calculateDerivedFields();

        Assert.assertEquals(mv.TOTAL_MENDELIAN_VIOLATIONS, mv.NUM_DIPLOID_DENOVO +
                mv.NUM_HOMVAR_HOMVAR_HET + mv.NUM_HOMREF_HOMVAR_HOM + mv.NUM_HOM_HET_HOM +
                mv.NUM_HAPLOID_DENOVO + mv.NUM_HAPLOID_OTHER + mv.NUM_OTHER );

        Assert.assertEquals(metrics, mv);

        Assert.assertEquals(grepMv(violationsNA12878, Diploid_Denovo.name()), mv.NUM_DIPLOID_DENOVO);
        Assert.assertEquals(grepMv(violationsNA12878, HomVar_HomVar_Het.name()), mv.NUM_HOMVAR_HOMVAR_HET);
        Assert.assertEquals(grepMv(violationsNA12878, HomRef_HomVar_Hom.name()), mv.NUM_HOMREF_HOMVAR_HOM);
        Assert.assertEquals(grepMv(violationsNA12878, Hom_Het_Hom.name()), mv.NUM_HOM_HET_HOM);
        Assert.assertEquals(grepMv(violationsNA12878, Haploid_Denovo.name()), mv.NUM_HAPLOID_DENOVO);
        Assert.assertEquals(grepMv(violationsNA12878, Haploid_Other.name()), mv.NUM_HAPLOID_OTHER);
        Assert.assertEquals(grepMv(violationsNA12878, Other.name()), mv.NUM_OTHER);

        metrics = summary.getMetrics().get(1);

        mv.FAMILY_ID = "FAKE";
        mv.OFFSPRING = "FAKE";
        mv.OFFSPRING_SEX = Sex.Male;
        mv.NUM_DIPLOID_DENOVO =    grep(vcfFile, "DiploidDenovoBOTH") +   grep(vcfFile, "DiploidDenovoFAKE");
        mv.NUM_HOMREF_HOMVAR_HOM = grep(vcfFile, "HomrefHomvarHomBOTH") + grep(vcfFile, "HomrefHomvarHomFAKE");
        mv.NUM_HOMVAR_HOMVAR_HET = grep(vcfFile, "HomvarHomvarHetBOTH") + grep(vcfFile, "HomvarHomvarHetFAKE");
        mv.NUM_HAPLOID_DENOVO =    grep(vcfFile, "HaploidDenovoBOTH") +   grep(vcfFile, "HaploidDenovoFAKE");
        mv.NUM_HAPLOID_OTHER =     grep(vcfFile, "HaploidOtherBOTH") +    grep(vcfFile, "HaploidOtherFAKE");
        mv.NUM_OTHER =             grep(vcfFile, "OtherBOTH") +           grep(vcfFile, "OtherFAKE");


        mv.calculateDerivedFields();
        Assert.assertEquals(metrics, mv);

        Assert.assertEquals(mv.TOTAL_MENDELIAN_VIOLATIONS, mv.NUM_DIPLOID_DENOVO +
                mv.NUM_HOMVAR_HOMVAR_HET + mv.NUM_HOMREF_HOMVAR_HOM + mv.NUM_HOM_HET_HOM +
                mv.NUM_HAPLOID_DENOVO + mv.NUM_HAPLOID_OTHER + mv.NUM_OTHER );

        Assert.assertEquals(grepMv(violationsFAKE, Diploid_Denovo.name()), mv.NUM_DIPLOID_DENOVO);
        Assert.assertEquals(grepMv(violationsFAKE, HomVar_HomVar_Het.name()), mv.NUM_HOMVAR_HOMVAR_HET);
        Assert.assertEquals(grepMv(violationsFAKE, HomRef_HomVar_Hom.name()), mv.NUM_HOMREF_HOMVAR_HOM);
        Assert.assertEquals(grepMv(violationsFAKE, Hom_Het_Hom.name()), mv.NUM_HOM_HET_HOM);
        Assert.assertEquals(grepMv(violationsFAKE, Haploid_Denovo.name()), mv.NUM_HAPLOID_DENOVO);
        Assert.assertEquals(grepMv(violationsFAKE, Haploid_Other.name()), mv.NUM_HAPLOID_OTHER);
        Assert.assertEquals(grepMv(violationsFAKE, Other.name()), mv.NUM_OTHER);

        TestUtil.recursiveDelete(directoryForViolations);
    }

    /** returns the number of lines in the file that contain a regular expression (decorated with "MV=" and
     * expected to be in an INFO field in a vcf)
     *
     * @param file File to examine
     * @param regex String containing a regular expression to look for in the file
     * @return the number of lines that contain regex
     */
    private int grep(final File file, final String regex) {

        int results = 0;
        final Pattern pattern = Pattern.compile(".*"+regex+".*");
        try (final LineIteratorImpl li = new LineIteratorImpl(new AsciiLineReader(IOUtil.openFileForReading(file)))) {

            while (li.hasNext()) {
                final String line = li.next();
                if (pattern.matcher(line).matches()) {
                    results++;
                }
            }
        } catch (final IOException e) {
            e.printStackTrace();
        }
        return results;
    }

    private int grepMv(final File file, final String regex) {
        return grep(file, "[;\t]MV="+regex+"[;\t]");
    }
}