package picard.analysis;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class AlignmentSummaryMetricsCollectorTest {

    @DataProvider
    public Object[][] get3PrimeSoftClipData() {
        return new Object[][]{
                {"20M", true, 0},
                {"20M", false, 0},

                {"5H2S", true, 0},
                {"5H2S", false, 0},
                {"5H2S6H", false, 0},
                {"5H2S6H", true, 0},

                {"2H20M", true, 0},
                {"2H20M", false, 0},

                {"2S20M", false, 0},
                {"2S20M", true, 2},
                {"20M2S", true, 0},
                {"20M2S", false, 2},

                {"20M2S4S", false, 6},
                {"20M2S4S", true, 0},

                {"2S20M5S", false, 5},
                {"2S20M4S", true, 2},
                {"4S20M2S", true, 4},
                {"4S20M2S", false, 2},
                {"2H3S10M4S5H", false,4},
                {"2H3S10M4S5H", true,3},
        };
    }

    @Test(dataProvider = "get3PrimeSoftClipData")
    public void testGet3PrimeSoftClip(final String cigarString, final boolean negativeStrand, final int expectedSoftClips) {
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        final int resultant5PrimeSoftclips = AlignmentSummaryMetricsCollector.get3PrimeSoftClippedBases(cigar, negativeStrand);
        Assert.assertEquals(resultant5PrimeSoftclips, expectedSoftClips, cigarString + "-" + negativeStrand);
    }
}
