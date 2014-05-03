package picard.sam.testers;


import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.TestUtil;
import org.testng.Assert;
import picard.cmdline.CommandLineProgram;
import picard.sam.MarkDuplicates;

import java.io.File;

/**
 * This class is an extension of SamFileTester used to test MarkDuplicates with SAM files generated on the fly.
 */
public class MarkDuplicatesTester extends SamFileTester {

    private final MarkDuplicates program = new MarkDuplicates();

    public MarkDuplicatesTester() {
        super(50, true);

        final File metrics = new File(getOutputDir(), "metrics.txt");
        addArg("METRICS_FILE=" + metrics);
    }

    @Override
    public void test() {
        try {
            // Read the output and check the duplicate flag
            final SAMFileReader reader = new SAMFileReader(getOutput());
            for (final SAMRecord record : reader) {
                final String key = samRecordToDuplicatesFlagsKey(record);
                Assert.assertTrue(this.duplicateFlags.containsKey(key));
                final boolean value = this.duplicateFlags.get(key);
                this.duplicateFlags.remove(key);
                if (value != record.getDuplicateReadFlag()) {
                    System.err.println("Mismatching read:");
                    System.err.print(record.getSAMString());
                }
                Assert.assertEquals(record.getDuplicateReadFlag(), value);
            }
            reader.close();
        } finally {
            TestUtil.recursiveDelete(getOutputDir());
        }
    }

    @Override
    protected CommandLineProgram getProgram() {
        return program;
    }
}

