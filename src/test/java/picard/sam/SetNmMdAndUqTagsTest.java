package picard.sam;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class SetNmMdAndUqTagsTest {

    private static final File fasta = new File("testdata/picard/sam/merger.fasta");
    @DataProvider(name="filesToFix")
    Object[][] testValidSortData() {
        return new Object[][]{
                new Object[]{new File("testdata/picard/sam/aligned.sam"), fasta},
                new Object[]{new File("testdata/picard/sam/aligned_queryname_sorted.sam"), fasta},
                new Object[]{new File("testdata/picard/sam/aligned_queryname_sorted.sam"), fasta},
        };
    }

    @Test(dataProvider = "filesToFix")
    public void testValidSort(final File input, final File reference) throws IOException {
        final File sortOutput = File.createTempFile("Sort", ".bam");
        sortOutput.deleteOnExit();
        final File fixOutput = File.createTempFile("Fix", ".bam");
        fixOutput.deleteOnExit();
        final File validateOutput = File.createTempFile("Sort", ".validation_report");
        validateOutput.deleteOnExit();

        sort(input, sortOutput);
        fixFile(sortOutput, fixOutput, reference);
        validate(fixOutput, validateOutput, reference, false);
    }

    @Test(dataProvider = "filesToFix")
    public void testUQValidSort(final File input, final File reference) throws IOException {
        final File sortOutput = File.createTempFile("Sort", ".bam");
        sortOutput.deleteOnExit();
        final File fixOutput = File.createTempFile("Fix", ".bam");
        fixOutput.deleteOnExit();
        final File validateOutput = File.createTempFile("Sort", ".validation_report");
        validateOutput.deleteOnExit();

        sort(input, sortOutput);
        calcUQ(sortOutput, fixOutput, reference);
        //ignore warnings because a bam with no NM tag throws warnings
        validate(fixOutput, validateOutput, reference, true);
    }

    private void validate(final File input, final File output, final File reference, final boolean ignoreWarnings) {

        final String[] args = {
                "INPUT=" + input,
                "OUTPUT=" + output,
                "MODE=VERBOSE",
                "REFERENCE_SEQUENCE=" + reference,
                "IGNORE_WARNINGS=" + ignoreWarnings};

        ValidateSamFile validateSam = new ValidateSamFile();
        Assert.assertEquals(validateSam.instanceMain(args), 0, "Validate did not succeed");
    }

    private void sort(final File input, final File output) {

        final String[] args = {
                "INPUT=" + input,
                "OUTPUT=" + output,
                "SORT_ORDER=coordinate"
               };

        SortSam sortSam = new SortSam();
        Assert.assertEquals(sortSam.instanceMain(args), 0, "Sort did not succeed");
    }

    private void fixFile(final File input, final File output, final File reference) throws IOException {

        final String[] args = {
                "INPUT=" + input,
                "OUTPUT=" + output,
                "REFERENCE_SEQUENCE="+reference };

        SetNmMdAndUqTags setNmMdAndUqTags = new SetNmMdAndUqTags();
        Assert.assertEquals(setNmMdAndUqTags.instanceMain(args), 0, "Fix did not succeed");
    }

    private void calcUQ(final File input, final File output, final File reference) throws IOException {

        final String[] args = {
                "INPUT=" + input,
                "OUTPUT=" + output,
                "REFERENCE_SEQUENCE=" + reference };

        CalculateUqTag calculateUqTag = new CalculateUqTag();
        Assert.assertEquals(calculateUqTag.instanceMain(args), 0, "Uq calc did not succeed");
    }
}
