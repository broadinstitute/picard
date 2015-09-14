package picard.cmdline;

import org.testng.annotations.Test;

import java.util.Collections;

/**
 * Created by farjoun on 9/10/15.
 */
public class PicardCommandLineTest {

    @Test
    public void TestPicardPublic() { // this tests fails if any CLP in picard is missing its @CommandLineProgramProperties annotation
        PicardCommandLine picardCommandLine = new PicardCommandLine();
        picardCommandLine.instanceMain(new String[]{""});
    }

}