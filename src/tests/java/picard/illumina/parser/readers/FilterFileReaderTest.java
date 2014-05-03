package net.sf.picard.illumina.parser.readers;

import net.sf.picard.PicardException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.NoSuchElementException;

public class FilterFileReaderTest {
    public static File TEST_DATA_DIR = new File("testdata/net/sf/picard/illumina/readerTests");
    public static final File PASSING_FILTER_FILE = new File(TEST_DATA_DIR, "pf_passing.filter");

    public static final boolean [] expectedPfs = {
            false, false, false, false,     true,  false, false, false,      true,  true,  false, false,   true,  false, true,  false,
            true,  false, false, true,      true,  true,  true,  false,      true,  true,  false, true,    true,  true,  true,  false,
            true,  true,  true,  true,      false, true,  false, false,      false, true,  true,  false,   false, true,  false, true,
            false, true,  true,  true,      false, false, true,  false,      false, false, true,  true,    false, false, false, true
    };

    @Test
    public void readValidFile() {
        final FilterFileReader reader = new FilterFileReader(PASSING_FILTER_FILE);
        Assert.assertEquals(reader.numClusters, expectedPfs.length);
        for(int i = 0; i < reader.numClusters; i++) {
            Assert.assertEquals(reader.hasNext(), true);
            Assert.assertEquals(reader.next().booleanValue(),    expectedPfs[i]);
        }

        Assert.assertEquals(false, reader.hasNext());
    }

    @Test(expectedExceptions = NoSuchElementException.class)
    public void readPastEnd() {
        final FilterFileReader reader = new FilterFileReader(PASSING_FILTER_FILE);
        for(int i = 0; i < reader.numClusters; i++) {
            reader.next();
        }

        Assert.assertEquals(false, reader.hasNext());
        reader.next();
    }

    @DataProvider(name="failingFiles")
    public Object[][] failingFiles() {
        return new Object[][] {
            {"pf_failing1.filter"},
            {"pf_failing2.filter"},
            {"pf_tooLarge.filter"},
            {"pf_tooShort.filter"},
            {"pf_badOpeningBytes.filter"},
            {"pf_badVersionBytes.filter"},
            {"pf_nonExistentFile.filter"}
        };
    }

    @Test(dataProvider = "failingFiles", expectedExceptions = PicardException.class)
    public void readInvalidValues(final String failingFile) {
        final FilterFileReader reader = new FilterFileReader(new File(TEST_DATA_DIR, failingFile));
        while(reader.hasNext()) {
            reader.next();
        }
    }
}