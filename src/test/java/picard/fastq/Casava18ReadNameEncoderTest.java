package picard.fastq;

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.illumina.parser.ClusterData;

public class Casava18ReadNameEncoderTest {
    @Test public void testGenerateReadName() {
        final Casava18ReadNameEncoder encoder = new Casava18ReadNameEncoder("INST", "run", "12345ACXX");
        final ClusterData data = new ClusterData();
        data.setLane(1);
        data.setTile(123);
        data.setX(1046);
        data.setY(1149);
        data.setMatchedBarcode("ACGTAATT");
        data.setPf(true);

        final String unpaired = encoder.generateReadName(data, null);
        final String r1 = encoder.generateReadName(data, 1);
        final String r2 = encoder.generateReadName(data, 2);

        Assert.assertEquals(unpaired, "INST:run:12345ACXX:1:123:1046:1149 :N:0:ACGTAATT");
        Assert.assertEquals(r1,       "INST:run:12345ACXX:1:123:1046:1149 1:N:0:ACGTAATT");
        Assert.assertEquals(r2,       "INST:run:12345ACXX:1:123:1046:1149 2:N:0:ACGTAATT");
    }
}
