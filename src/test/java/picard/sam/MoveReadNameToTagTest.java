package picard.sam;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;

public class MoveReadNameToTagTest extends CommandLineProgramTest {
    private static final File INPUT_FILE = new File("testdata/picard/sam/MoveReadNameToTag/sorted-pair.sam");

    public String getCommandLineProgramName() {
        return MoveReadNameToTag.class.getSimpleName();
    }

    @Test
    public void testOkFile() throws IOException {
        final String tmpDir = IOUtil.getDefaultTmpDir().getAbsolutePath();

        final String[] args = new String[]{
                "INPUT=" + INPUT_FILE.getAbsolutePath(),
                //"IN_RN=true",
                "TILE_TAG=true",
                "XY_FULL=true",
                "OUTPUT=" + tmpDir + "/out.sam"
        };

        runPicardCommandLine(args);
        final File out = new File(tmpDir, "out.sam");
        out.deleteOnExit();
        final SamReader reader = SamReaderFactory.makeDefault().open(out);

        for (SAMRecord rec : reader) {
            Assert.assertEquals(rec.getAttribute("XT"), 1101);
            Assert.assertEquals(rec.getAttribute("XX"), 10021);
            Assert.assertEquals(rec.getAttribute("XY"), 75823);
        }
    }
}
