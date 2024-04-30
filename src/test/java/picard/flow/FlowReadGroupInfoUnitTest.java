package picard.flow;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.FileSystems;
import java.nio.file.Path;


public class FlowReadGroupInfoUnitTest {
    @Test
    void testReadGroupParsing(){
        final String    testResourceDir =  "testdata/picard/flow/reads/";
        final String inputDir = testResourceDir + "/input/";

        final Path inputFile = FileSystems.getDefault().getPath(inputDir, "sample_mc.bam");
        final SamReader reader = SamReaderFactory.makeDefault().open(new File(inputFile.toString()));
        SAMFileHeader header = reader.getFileHeader();
        SAMReadGroupRecord rg1 = header.getReadGroup("UGAv3-72");
        FlowReadGroupInfo frg1 = new FlowReadGroupInfo(rg1);
        assert(frg1.maxClass==12);
        assert(frg1.flowOrder.startsWith("TGCA"));
        assert(frg1.flowOrder.startsWith("TGCA"));
        assert(frg1.isFlowPlatform);

        SAMReadGroupRecord rg2 = header.getReadGroup("UGAv3-73");
        FlowReadGroupInfo frg2 = new FlowReadGroupInfo(rg2);
        assert(frg2.maxClass==20);
        assert(frg2.flowOrder.startsWith("TGCA"));
        assert(frg2.isFlowPlatform);
        SAMReadGroupRecord rg3 = header.getReadGroup("UGAv3-74");
        FlowReadGroupInfo frg3 = new FlowReadGroupInfo(rg3);
        assert(!frg3.isFlowPlatform);
    }
}
