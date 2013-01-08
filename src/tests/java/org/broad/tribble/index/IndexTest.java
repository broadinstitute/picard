package org.broad.tribble.index;

import org.broad.tribble.TestUtils;
import org.broad.tribble.index.linear.LinearIndex;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


public class IndexTest {
    private final static String CHR = "1";
    private final static File MassiveIndexFile = new File(TestUtils.DATA_DIR + "Tb.vcf.idx");

    @DataProvider(name = "StartProvider")
    public Object[][] makeStartProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();

//        for ( int mid = 0; mid <= end; mid += 1000000 ) {
//            tests.add(new Object[]{0, mid, mid+1000000, end});
//        }

        tests.add(new Object[]{1226943, 1226943, 1226943, 2000000});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "StartProvider")
    public void testMassiveQuery(final int start, final int mid, final int mid2, final int end) throws IOException {
        LinearIndex index = (LinearIndex)IndexFactory.loadIndex(MassiveIndexFile.getAbsolutePath());

        final List<Block> leftBlocks = index.getBlocks(CHR, start, mid);
        final List<Block> rightBlocks = index.getBlocks(CHR, mid2, end); // gap must be big to avoid overlaps
        final List<Block> allBlocks = index.getBlocks(CHR, start, end);

        final long leftSize = leftBlocks.isEmpty() ? 0 : leftBlocks.get(0).getSize();
        final long rightSize = rightBlocks.isEmpty() ? 0 : rightBlocks.get(0).getSize();
        final long allSize = allBlocks.isEmpty() ? 0 : allBlocks.get(0).getSize();

        Assert.assertTrue(leftSize >= 0, "Expected leftSize to be positive " + leftSize);
        Assert.assertTrue(rightSize >= 0, "Expected rightSize to be positive " + rightSize);
        Assert.assertTrue(allSize >= 0, "Expected allSize to be positive " + allSize);

        Assert.assertTrue(allSize >= Math.max(leftSize,rightSize), "Expected size of joint query " + allSize + " to be at least >= max of left " + leftSize + " and right queries " + rightSize);
    }
}
