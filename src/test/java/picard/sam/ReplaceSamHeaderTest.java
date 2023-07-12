package picard.sam;

import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.util.IntervalListTools;

import static org.testng.Assert.*;

public class ReplaceSamHeaderTest extends CommandLineProgramTest {
    @Test
    public void testCloud(){

    }

    @Override
    public String getCommandLineProgramName() {
        return ReplaceSamHeader.class.getSimpleName();
    }
}