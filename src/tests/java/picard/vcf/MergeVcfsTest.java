package picard.vcf;

import picard.cmdline.CommandLineProgram;

/**
 * Created by bradt on 9/3/14.
 */
public class MergeVcfsTest extends AbstractVcfMergingTest {

    @Override
    protected CommandLineProgram getProgram() {
        return new MergeVcfs();
    }
}
