package picard.sam.testers;

import htsjdk.samtools.SAMValidationError;
import org.testng.Assert;
import picard.cmdline.CommandLineProgramTest;
import picard.sam.ValidateSamFile;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * Created by ggrant on 6/12/14.
 */
public class ValidateSamTester extends CommandLineProgramTest {
    private Collection<SAMValidationError.Type> ignoreErrors = Collections.emptyList();
    public String getCommandLineProgramName() {
        return ValidateSamFile.class.getSimpleName();
    }

    public void setIgnoreError(final Collection<SAMValidationError.Type> errorsToIgnore){
        ignoreErrors = new ArrayList<>(errorsToIgnore);
    }

    public void assertSamValid(final File samFile) {
        final List<String> args = new ArrayList<>();
        args.add("I=" + samFile.getAbsolutePath());
        ignoreErrors.forEach(s -> args.add("IGNORE=" + s));
        final int validateExitStatus = runPicardCommandLine(args.toArray(new String[0]));
        Assert.assertEquals(validateExitStatus, 0, "There were problems in file: " + samFile.getAbsolutePath());
    }
}
