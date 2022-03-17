package picard.vcf;

import htsjdk.variant.vcf.VCFFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgram;
import picard.nio.PicardHtsPath;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by bradt on 9/3/14.
 */
public class MergeVcfsTest extends AbstractVcfMergingClpTester {

    @Override
    protected CommandLineProgram getProgram() {
        return new MergeVcfs();
    }
    
    @Test
    public void TestComments() throws IOException {
        final String comment1 = "This is a comment";
        final List<String> args = new ArrayList<>();
        final File output = VcfTestUtils.createTemporaryIndexedFile("result", ".vcf");

        args.add("INPUT=" + new PicardHtsPath(TEST_DATA_PATH + "mini.vcf"));
        args.add("OUTPUT=" + output.getAbsolutePath());
        args.add("COMMENT=" + comment1);

        Assert.assertEquals(new MergeVcfs().instanceMain(args.toArray(new String[]{})), 0);

       try (VCFFileReader reader = new VCFFileReader(output, false)) {
            Assert.assertTrue(reader.getFileHeader().getMetaDataInInputOrder().stream().anyMatch(H->H.getKey().equals("MergeVcfs.comment") && H.getValue().equals(comment1)));
       }
    }
}
