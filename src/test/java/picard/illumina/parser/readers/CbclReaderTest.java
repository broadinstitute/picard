package picard.illumina.parser.readers;

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.illumina.parser.BclData;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class CbclReaderTest {

    public static final File TestDataDir = new File("testdata/picard/illumina/readerTests");
    public static final File PASSING_CBCL_C99_1 = new File(TestDataDir, "C99.1.cbcl");
    public static final File PASSING_CBCL_C100_1 = new File(TestDataDir, "C100.1.cbcl");
    public static final File TILE_1101_FILTER = new File(TestDataDir, "tile_1101.filter");
    public static final File TILE_1102_FILTER = new File(TestDataDir, "tile_1102.filter");

    public static final char[] expectedBases = new char[]{
            'G','G','C','C','G','A','A','G','T','A','C','C','C','T','G','A'
    };

    public static final int[] expectedQuals = new int[]{
            37,37,37,37,37,37,37,12,37,37,37,37,37,37,12,37
    };

    @Test
    public void testReadValidFile() {
        final CbclReader reader = new CbclReader(Arrays.asList(PASSING_CBCL_C99_1, PASSING_CBCL_C100_1), new int[] {2});
        Map<Integer, File> filters = new HashMap<>();
        filters.put(1101, TILE_1101_FILTER);
        filters.put(1102, TILE_1102_FILTER);
        reader.addFilters(filters);

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
