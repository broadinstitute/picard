package net.sf.picard.illumina.parser.readers;

import net.sf.picard.PicardException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class BclReaderTest {

    public static File TestDataDir = new File("testdata/net/sf/picard/illumina/readerTests");
    public static final File PASSING_BCL_FILE  = new File(TestDataDir, "bcl_passing.bcl");
    public static final File QUAL_0FAILING_BCL_FILE  = new File(TestDataDir, "bcl_failing.bcl");
    public static final File QUAL_1FAILING_BCL_FILE  = new File(TestDataDir, "bcl_failing2.bcl");
    public static final File FILE_TOO_LARGE  = new File(TestDataDir, "bcl_tooLarge.bcl");
    public static final File FILE_TOO_SHORT  = new File(TestDataDir, "bcl_tooShort.bcl");

    public static final char [] expectedBases = new char[]{
            'C', 'A', 'A', 'A',     'T', 'C', 'T', 'G',     'T', 'A', 'A', 'G',     'C', 'C', 'A', 'A', 
            'C', 'A', 'C', 'C',     'A', 'A', 'C', 'G',     'A', 'T', 'A', 'C',     'A', 'A', 'C', 'A', 
            'T', 'G', 'C', 'A',     'C', 'A', 'A', 'C',     'G', 'C', 'A', 'A',     'G', 'T', 'G', 'C', 
            'A', 'C', 'G', 'T',     'A', 'C', 'A', 'A',     'C', 'G', 'C', 'A',     'C', 'A', 'T', 'T', 
            'T', 'A', 'A', 'G',     'C', 'G', 'T', 'C',     'A', 'T', 'G', 'A',     'G', 'C', 'T', 'C', 
            'T', 'A', 'C', 'G',     'A', 'A', 'C', 'C',     'C', 'A', 'T', 'A',     'T', 'G', 'G', 'G', 
            'C', 'T', 'G', 'A',     'A', '.', '.', 'G',     'A', 'C', 'C', 'G',     'T', 'A', 'C', 'A', 
            'G', 'T', 'G', 'T',     'A', '.'
    };

    public static final int [] expectedQuals = new int[]{
        18, 29,  8, 17,     27, 25, 28, 27,      9, 29,  8, 20,     25, 24, 27, 27,
        30,  8, 19, 24,     29, 29, 25, 28,      8, 29, 26, 24,     29,  8, 18,  8,
        29, 28, 26, 29,     25, 8,  26, 25,     28, 25,  8, 28,     28, 27, 29, 26,
        25, 26, 27, 25,      8, 18,  8, 26,     24, 29, 25,  8,     24,  8, 25, 27,
        27, 25,  8, 28,     24, 27, 25, 25,      8, 27, 25,  8,     16, 24, 28, 25,
        28,  8, 24, 27,     25,  8, 20, 29,     24, 27, 28,  8,     23, 10, 23, 11,
        15, 11, 10, 12,     12,  2,  2, 31,     24,  8,  4, 36,     12, 17, 21,  4,
         8, 12, 18, 23,     27,  2
    };

    public byte [] qualsAsBytes() {
        final byte [] byteVals = new byte[expectedQuals.length];
        for(int i = 0; i < byteVals.length; i++) {
            byteVals[i] = (byte)expectedQuals[i];
        }
        return byteVals;
    }

    @Test
    public void readValidFile() {
        final BclReader reader = new BclReader(PASSING_BCL_FILE);
        final byte [] quals = qualsAsBytes();

        Assert.assertEquals(reader.numClusters, expectedBases.length);

        int readNum = 0;
        while(readNum < reader.numClusters) {
            final BclReader.BclValue bv = reader.next();
            Assert.assertEquals(bv.base,    expectedBases[readNum], " On num cluster: " + readNum);
            Assert.assertEquals(bv.quality, quals[readNum], " On num cluster: " + readNum);
            ++readNum;
        }
    }

    @DataProvider(name="failingFiles")
    public Object[][] failingFiles() {
        return new Object[][] {
            {QUAL_0FAILING_BCL_FILE},
            {QUAL_1FAILING_BCL_FILE},
            {new File(TestDataDir, "SomeNoneExistantFile.bcl")},
            {FILE_TOO_LARGE},
            {FILE_TOO_SHORT}
        };
    }

    @Test(expectedExceptions = PicardException.class, dataProvider = "failingFiles")
    public void failingFileTest(final File failingFile) {
        final BclReader reader = new BclReader(failingFile);
        Assert.assertEquals(reader.numClusters, expectedBases.length);
        while(reader.hasNext()) {
            reader.next();
        }
    }
}
