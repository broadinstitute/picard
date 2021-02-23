package picard.illumina.parser;

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.illumina.parser.fakers.MultiTileBclFileFaker;
import picard.illumina.parser.readers.BclReader;

import java.io.File;

public class MultiTileBclFileFakerTest {
    private static final File TEST_DATA_DIRECTORY = new File("testdata/picard/illumina/parserTests/bciParser");

    @Test
    public void testFileLengthMatchesHeaderLength() throws Exception {
        final File fakeFile = File.createTempFile("MultiTileBclFileFakerTest", ".bcl.bgzf");
        final File fakeBciFile = new File(fakeFile.getAbsolutePath() + ".bci");
        fakeFile.deleteOnExit();
        fakeBciFile.deleteOnExit();

        new MultiTileBclFileFaker().fakeMultiTileBclFile(fakeFile, new TileIndex(new File(TEST_DATA_DIRECTORY, "0001.bcl.bgzf.bci")));

        Assert.assertEquals(1, BclReader.getNumberOfClusters(fakeFile));
    }
}