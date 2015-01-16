package picard.sam.testers;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SAMValidationError;
import htsjdk.samtools.SamFileValidator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.TestUtil;
import org.testng.Assert;
import picard.cmdline.CommandLineProgram;
import picard.sam.CleanSam;

import java.io.PrintWriter;
import java.util.Arrays;

/**
 * This class is the extension of the SamFileTester to test CleanSam with SAM files generated on the fly.
 */
public class CleanSamTester extends SamFileTester {
    private final String expectedCigar;

    @Override
    public String getCommandLineProgramName() { return CleanSam.class.getSimpleName(); }

    public CleanSamTester(final String expectedCigar, final int readLength, final int defaultChromosomeLength) {
        super(readLength, true, defaultChromosomeLength);
        this.expectedCigar = expectedCigar;
    }


    protected void test() {
        try {
            final SamFileValidator validator = new SamFileValidator(new PrintWriter(System.out), 8000);

            // Validate it has the expected cigar
            validator.setIgnoreWarnings(true);
            validator.setVerbose(true, 1000);
            validator.setErrorsToIgnore(Arrays.asList(SAMValidationError.Type.MISSING_READ_GROUP));
            SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
            SamReader samReader = factory.open(getOutput());
            final SAMRecordIterator iterator = samReader.iterator();
            while (iterator.hasNext()) {
                final SAMRecord rec = iterator.next();
                Assert.assertEquals(rec.getCigarString(), expectedCigar);
                if (SAMUtils.hasMateCigar(rec)) {
                    Assert.assertEquals(SAMUtils.getMateCigarString(rec), expectedCigar);
                }
            }
            CloserUtil.close(samReader);

            // Run validation on the output file
            samReader = factory.open(getOutput());
            final boolean validated = validator.validateSamFileVerbose(samReader, null);
            CloserUtil.close(samReader);

            Assert.assertTrue(validated, "ValidateSamFile failed");
        } finally {
            TestUtil.recursiveDelete(getOutputDir());
        }
    }
}
