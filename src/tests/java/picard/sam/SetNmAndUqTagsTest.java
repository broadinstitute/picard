package picard.sam;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class SetNmAndUqTagsTest {

    private static final File fasta = new File("testdata/picard/sam/merger.fasta");
    @DataProvider(name="filesToFix")
    Object[][] TestValidSortData() {
        return new Object[][]{
                new Object[]{new File("testdata/picard/sam/aligned.sam"), fasta},
                new Object[]{new File("testdata/picard/sam/aligned_queryname_sorted.sam"), fasta},
                new Object[]{new File("testdata/picard/sam/aligned_queryname_sorted.sam"), fasta},
        };
    }

    @Test(dataProvider = "filesToFix")
    public void TestValidSort(final File input, final File reference) throws IOException {
        final File sortOutput = File.createTempFile("Sort", ".bam");
        sortOutput.deleteOnExit();
        final File fixOutput = File.createTempFile("Fix", ".bam");
        fixOutput.deleteOnExit();
        final File validateOutput = File.createTempFile("Sort", ".validation_report");
        validateOutput.deleteOnExit();

        sort(input, sortOutput);
        fixFile(sortOutput, fixOutput, reference);
        validate(fixOutput,validateOutput, reference);
    }

    private void validate(final File input, final File output, final File reference) {

        final String[] args = new String[] {
                "INPUT="+input,
                "OUTPUT="+output,
                "MODE=VERBOSE",
                "REFERENCE_SEQUENCE="+reference };

        ValidateSamFile validateSam = new ValidateSamFile();
        Assert.assertEquals(validateSam.instanceMain(args), 0, "validate did not succeed");
    }

    private void sort(final File input, final File output) {

        final String[] args = new String[] {
                "INPUT=" + input,
                "OUTPUT=" + output,
                "SORT_ORDER=coordinate"
               };

        SortSam sortSam = new SortSam();
        Assert.assertEquals(sortSam.instanceMain(args), 0, "Sort did not succeed");
    }

    private void fixFile(final File input, final File output, final File reference) throws IOException {

        final String[] args = new String[] {
                "INPUT="+input,
                "OUTPUT="+output,
                "REFERENCE_SEQUENCE="+reference };

        SetNmAndUqTags setNmAndUqTags = new SetNmAndUqTags();
        Assert.assertEquals(setNmAndUqTags.instanceMain(args), 0, "Fix did not succeed");
    }
}
