/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package picard.sam;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;

public class CompareSAMsTest extends CommandLineProgramTest {
    private static final File TEST_FILES_DIR = new File("testdata/picard/sam/CompareSAMs");

    public String getCommandLineProgramName() {
        return CompareSAMs.class.getSimpleName();
    }

    @DataProvider(name="compareSams")
    public Object[][] compareSamsTestData() {
        final ArrayList<String> lenientArgs = new ArrayList<>(Arrays.asList("LENIENT_DUP=true", "MQ0_MATCH=true"));
        return new Object[][] {
                {"genomic_sorted.sam", "unsorted.sam", null, false},
                {"genomic_sorted.sam", "chr21.sam", null, false},
                {"genomic_sorted.sam", "bigger_seq_dict.sam", null, false},
                {"bigger_seq_dict.sam", "bigger_seq_dict.sam", null, true},
                {"genomic_sorted.sam", "genomic_sorted.sam", null, true},
                {"genomic_sorted.sam", "has_non_primary.sam", null, true},
                {"genomic_sorted_5.sam", "genomic_sorted_5_plus.sam", null, false},
                {"group_same_coord.sam", "group_same_coord_diff_order.sam", null, false},
                {"genomic_sorted_same_position.sam", "genomic_sorted_same_position.sam", null, true},
                {"group_same_coord.sam", "diff_coords.sam", null, false},
                {"genomic_sorted.sam", "unmapped_first.sam", null, false},
                {"genomic_sorted.sam", "unmapped_second.sam", null, false},
                {"unmapped_first.sam", "unmapped_second.sam", null, false},
                {"unmapped_first.sam", "unmapped_first.sam", null, true},
                {"genomic_sorted.sam", "genomic_sorted_sam_v1.6.sam", null, false},
                {"unsorted.sam", "unsorted.sam", null, true},
                {"unsorted.sam", "unsorted2.sam", null, false},
                {"dup1.sam", "dup2.sam", lenientArgs, true},
                {"dup1.sam", "dup3.sam", lenientArgs, false},
                {"dup1.sam", "dup4.sam", lenientArgs, false},
                {"dup1.sam", "dup5.sam", lenientArgs, false},
                {"dup1.sam", "dup2.sam", null, false},
                {"dup1_queryname.sam", "dup2_queryname.sam", lenientArgs, true},
                {"dup1_queryname.sam", "dup3_queryname.sam", lenientArgs, false},
                {"dup1_queryname.sam", "dup4_queryname.sam", lenientArgs, false},
                {"dup1_queryname.sam", "dup5_queryname.sam", lenientArgs, false},
                {"dup1_queryname.sam", "dup2_queryname.sam", null, false},
                {"genomic_sorted.sam", "mq0_2.sam", lenientArgs, false},
                {"mq0_1.sam", "mq0_2.sam", lenientArgs, true},
                {"mq0_1.sam", "mq0_2.sam", null, false}
        };
    }

    @Test(dataProvider = "compareSams")
    public void testComparisons(final String f1, final String f2, final ArrayList<String> args, final boolean areEqual) throws IOException {
        final Path tmpOutput = Files.createTempFile("compareSam", ".tsv");
        final String in1 = new File(TEST_FILES_DIR, f1).getAbsolutePath();
        final String in2 = new File(TEST_FILES_DIR, f2).getAbsolutePath();
        ArrayList<String> commandArgs = new ArrayList<>(
                Arrays.asList(
                        "I=" + in1,
                        "I=" + in2,
                        "O=" + tmpOutput
                )
        );
        if (args != null) {
            commandArgs.addAll(args);
        }
        Assert.assertEquals(runPicardCommandLine(commandArgs) == 0, areEqual);
        final MetricsFile<SamComparisonMetric, Comparable<?>> metricsOutput = new MetricsFile<>();
        metricsOutput.read(new FileReader(tmpOutput.toFile()));

        //swap order of input files
        commandArgs = new ArrayList<>(
                Arrays.asList(
                        "I=" + in2,
                        "I=" + in1,
                        "O=" + tmpOutput
                )
        );
        if (args != null) {
            commandArgs.addAll(args);
        }
        Assert.assertEquals(runPicardCommandLine(commandArgs) == 0, areEqual);
        metricsOutput.read(new FileReader(tmpOutput.toFile()));

        Assert.assertEquals(metricsOutput.getMetrics().get(0).leftFile, in1);
        Assert.assertEquals(metricsOutput.getMetrics().get(0).rightFile, in2);

        Assert.assertEquals(metricsOutput.getMetrics().get(1).leftFile, in2);
        Assert.assertEquals(metricsOutput.getMetrics().get(1).rightFile, in1);
    }
}
