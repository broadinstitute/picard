package picard.cmdline;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

/**
 * Test that Picard CLPs don't log an ERROR just from starting up, which started happening
 * after we enabled the Intel inflater / deflater because of a missing configuration file.
 */
public class CommandLineProgramStartupErrorLogTest {

    @CommandLineProgramProperties(
            usage = "No-op CLP for testing initialization behavior.",
            usageShort = "No-op CLP for testing initialization behavior."
    )
    private class MockCLP extends CommandLineProgram {
        @Override
        public int doWork() { return 0; }
    }

    @Test
    public void testNoStartupErrorLog() {
        ByteArrayOutputStream stderrStream = new ByteArrayOutputStream();
        PrintStream newStderr = new PrintStream(stderrStream);
        PrintStream oldStderr = System.err;
        System.setErr(newStderr);

        CommandLineProgram clp = new MockCLP();

        try {
            clp.instanceMain(new String[]{});
            Assert.assertFalse(stderrStream.toString().contains("ERROR"), stderrStream.toString());
        } finally {
            System.setErr(oldStderr);
        }
    }
}
