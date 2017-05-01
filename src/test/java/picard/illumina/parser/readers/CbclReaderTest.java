package picard.illumina.parser.readers;

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.illumina.parser.BclData;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class CbclReaderTest {

    public static final File TestDataDir = new File("testdata/picard/illumina/readerTests/cbcls");
    public static final File PASSING_CBCL_C99_1 = new File(TestDataDir + "/C1.1", "L001_1.cbcl");
    public static final File PASSING_CBCL_C100_1 = new File(TestDataDir + "/C2.1", "L001_1.cbcl");
    public static final File TILE_1101_FILTER = new File(TestDataDir, "tile_1101.filter");
    public static final File TILE_1102_FILTER = new File(TestDataDir, "tile_1102.filter");

    public static final char[] expectedBases = new char[]{
            'G', 'G', 'C', 'C', 'G', 'A', 'A', 'G'
    };

    public static final int[] expectedQuals = new int[]{
            37, 37, 37, 37, 37, 37, 37, 12
    };

    @Test
    public void testReadValidFile() {
        Map<Integer, File> filters = new HashMap<>();
        filters.put(1101, TILE_1101_FILTER);
        LocsFileReader locsFileReader = new LocsFileReader(new File("testdata/picard/illumina/readerTests/s_1_6.locs"));
        List<AbstractIlluminaPositionFileReader.PositionInfo> locs = locsFileReader.toList();

        final CbclReader reader = new CbclReader(Arrays.asList(PASSING_CBCL_C99_1, PASSING_CBCL_C100_1), filters, new int[]{2}, 1101, locs);

        int i = 0;
        while (reader.hasNext()) {
            final BclData bv = reader.next();
            for (int cluster = 0; cluster < bv.bases.length; cluster++) {
                for (int cycle = 0; cycle < bv.bases[cluster].length; cycle++) {
                    String actual = new String(new byte[]{bv.bases[cluster][cycle]});
                    String expected = new String(new char[]{expectedBases[i]});
                    Assert.assertEquals(actual, expected, "For cluster " + cluster + " cycle " + cycle + ",");
                    Assert.assertEquals(bv.qualities[cluster][cycle], expectedQuals[i], "For cluster " + cluster + " cycle " + cycle + ",");
                    i++;
                }
            }
        }
        Assert.assertEquals(i, expectedBases.length);
        reader.close();
    }
}
