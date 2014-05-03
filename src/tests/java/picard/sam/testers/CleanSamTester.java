package net.sf.picard.sam.testers;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.sam.CleanSam;
import net.sf.picard.sam.SamFileValidator;
import net.sf.picard.sam.testers.SamFileTester;
import net.sf.samtools.*;
import net.sf.samtools.util.TestUtil;
import org.testng.Assert;

import java.io.PrintWriter;
import java.util.Arrays;

/**
 * This class is the extension of the SamFileTester to test CleanSam with SAM files generated on the fly.
 */
public class CleanSamTester extends SamFileTester {
    private final String expectedCigar;
    private final CleanSam program = new CleanSam();

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
            SAMFileReader samReader = new SAMFileReader(getOutput());
            samReader.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
            final SAMRecordIterator iterator = samReader.iterator();
            while (iterator.hasNext()) {
                final SAMRecord rec = iterator.next();
                Assert.assertEquals(rec.getCigarString(), expectedCigar);
                if (SAMUtils.hasMateCigar(rec)) {
                    Assert.assertEquals(SAMUtils.getMateCigarString(rec), expectedCigar);
                }
            }
            samReader.close();

            // Run validation on the output file
            samReader = new SAMFileReader(getOutput());
            final boolean validated = validator.validateSamFileVerbose(samReader, null);
            samReader.close();
            Assert.assertTrue(validated, "ValidateSamFile failed");
        } finally {
            TestUtil.recursiveDelete(getOutputDir());
        }
    }

    @Override
    protected CommandLineProgram getProgram() {
        return program;
    }

}
