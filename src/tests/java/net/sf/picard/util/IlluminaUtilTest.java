package net.sf.picard.util;

import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;
import org.testng.Assert;

/**
 * Just some simple tests for the IlluminaUtil.getTileFromReadName() method for now.
 *
 * @author Tim Fennell
 */
public class IlluminaUtilTest {
    @Test(dataProvider="readNames") public void testFindTileInReadName(final String readName, final Integer tile) {
        final Integer otherTile = IlluminaUtil.getTileFromReadName(readName);
        Assert.assertEquals(otherTile, tile, "Tile numbers do not match for read name: " + readName);
    }

    @Test public void performanceTestGetTileFromReadName() {
        final int ITERATIONS = 5000000;

        final long startTime = System.currentTimeMillis();
        for (int i=0; i<ITERATIONS; ++i) {
            final Integer tile = IlluminaUtil.getTileFromReadName("300WFAAXX090909:1:1:1024:978#0/1");
            if (tile == null || tile != 1) throw new RuntimeException("WTF?");
        }
        final long endTime = System.currentTimeMillis();

        System.out.println("Time taken: " + (endTime-startTime) + "ms.");
    }

    @DataProvider(name="readNames")
    public Object[][] readNames() {
        return new Object[][] {
                new Object[] {"300WFAAXX:1:119:1024:978#0/1", 119},
                new Object[] {"300WFAAXX090909:1:1:1024:978#0/1", 1},
                new Object[] {"FOO", null},
                new Object[] {"FOO:BAR_splat", null}
        };
    }
}
