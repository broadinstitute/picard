package picard.analysis;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class CollectFlowBaisedMetricsTest {


    @DataProvider
    Object[][] readLengthGivenFlowData(){
        return new Object[][]{
                new Object[]{"ACGTACGTACGT", "ACGT", 5, 5},
                new Object[]{"ACGTACGTACGT", "AGCT", 5, 2},
                new Object[]{"ACGTACGTACGT", "TGCA", 20, 6},
                new Object[]{"AACGGGTATTCGAATACGT", "ACGT", 8, 10},
                new Object[]{"AACGGGTATTCGAATACGT", "ACGTAGCT", 8, 10},
        };
    }
    @Test(dataProvider = "readLengthGivenFlowData")
    public void testReadLengthGivenFlow(String reference, String flow, int flowLength, int readLength) {
        Assert.assertEquals(CollectFlowBiasMetrics.readLengthGivenFlow(reference,flow.getBytes(),flowLength),readLength);
    }
}