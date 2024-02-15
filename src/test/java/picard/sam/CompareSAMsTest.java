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
import htsjdk.samtools.util.Histogram;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class CompareSAMsTest extends CommandLineProgramTest {
    private static final File TEST_FILES_DIR = new File("testdata/picard/sam/CompareSAMs");

    public String getCommandLineProgramName() {
        return CompareSAMs.class.getSimpleName();
    }

    @DataProvider(name="compareSams")
    public Object[][] compareSamsTestData() {
        final ArrayList<String> lenientArgs = new ArrayList<>(Arrays.asList("LENIENT_DUP=true", "LENIENT_LOW_MQ_ALIGNMENT=true"));
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
                {"duplicate_base.sam", "duplicate_four_mismatch_strict.sam", lenientArgs, true},
                {"duplicate_base.sam", "duplicate_four_mismatch_lenient_one_align_differ.sam", lenientArgs, false},
                {"duplicate_base.sam", "duplicate_two_mismatch_lenient.sam", lenientArgs, false},
                {"duplicate_base.sam", "duplicate_four_mismatch_lenient.sam", lenientArgs, false},
                {"duplicate_base.sam", "duplicate_four_mismatch_strict.sam", null, false},
                {"duplicate_base_queryname.sam", "duplicate_four_mismatch_strict_queryname.sam", lenientArgs, true},
                {"duplicate_base_queryname.sam", "duplicate_four_mismatch_lenient_one_align_differ_queryname.sam", lenientArgs, false},
                {"duplicate_base_queryname.sam", "duplicate_two_mismatch_lenient_queryname.sam", lenientArgs, false},
                {"duplicate_base_queryname.sam", "duplicate_four_mismatch_lenient_queryname.sam", lenientArgs, false},
                {"duplicate_base_queryname.sam", "duplicate_four_mismatch_strict_queryname.sam", null, false},
                {"genomic_sorted.sam", "mq0_2.sam", lenientArgs, false},
                {"mq0_1.sam", "mq0_2.sam", lenientArgs, true},
                {"mq0_1.sam", "mq0_2.sam", null, false}
        };
    }

    @Test(dataProvider = "compareSams")
    public void testComparisons(final String f1, final String f2, final ArrayList<String> args, final boolean areEqual) throws IOException {
        final File tmpOutput = File.createTempFile("compareSam", ".tsv");
        tmpOutput.deleteOnExit();
        final String in1 = new File(TEST_FILES_DIR, f1).getAbsolutePath();
        final String in2 = new File(TEST_FILES_DIR, f2).getAbsolutePath();
        ArrayList<String> commandArgs = new ArrayList<>(
                Arrays.asList(
                        in1,
                        in2,
                        "O=" + tmpOutput
                )
        );
        if (args != null) {
            commandArgs.addAll(args);
        }
        Assert.assertEquals(runPicardCommandLine(commandArgs) == 0, areEqual);
        final MetricsFile<SamComparisonMetric, Comparable<?>> metricsOutput = new MetricsFile<>();
        metricsOutput.read(new FileReader(tmpOutput));
        Assert.assertEquals(metricsOutput.getNumHistograms(), 0);

        //swap order of input files
        commandArgs = new ArrayList<>(
                Arrays.asList(
                        in2,
                        in1,
                        "O=" + tmpOutput
                )
        );
        if (args != null) {
            commandArgs.addAll(args);
        }
        Assert.assertEquals(runPicardCommandLine(commandArgs) == 0, areEqual);
        metricsOutput.read(new FileReader(tmpOutput));
        Assert.assertEquals(metricsOutput.getNumHistograms(), 0);

        Assert.assertEquals(metricsOutput.getMetrics().get(0).LEFT_FILE, in1);
        Assert.assertEquals(metricsOutput.getMetrics().get(0).RIGHT_FILE, in2);

        Assert.assertEquals(metricsOutput.getMetrics().get(1).LEFT_FILE, in2);
        Assert.assertEquals(metricsOutput.getMetrics().get(1).RIGHT_FILE, in1);
    }

    @DataProvider(name="compareSamsMQConcordance")
    public Object[][] compareSamsMQConcordanceTestData() {
        return new Object[][] {
                {"genomic_sorted.sam", "unsorted.sam", null},
                {"genomic_sorted.sam", "chr21.sam", null},
                {"genomic_sorted.sam", "bigger_seq_dict.sam", null},
                {"genomic_sorted.sam", "genomic_sorted.sam", new Object[][] { {"20,20", 1}, {"30,30", 1}}},
                {"genomic_sorted.sam", "has_non_primary.sam", new Object[][] { {"20,20", 1}, {"30,30", 1}}},
                {"genomic_sorted_5.sam", "genomic_sorted_5_plus.sam", new Object[][] { {"20,20", 1}, {"30,30", 4}}},
                {"group_same_coord.sam", "group_same_coord_diff_order.sam", new Object[][] { {"20,20", 1}, {"30,30", 2}}},
                {"genomic_sorted_same_position.sam", "genomic_sorted_same_position.sam", new Object[][] { {"0,0", 2}}},
                {"group_same_coord.sam", "diff_coords.sam", new Object[][] { {"20,20", 1}, {"30,30", 4}}},
                {"genomic_sorted.sam", "unmapped_first.sam", new Object[][] { {"20,0", 1}, {"30,30", 1}}},
                {"genomic_sorted.sam", "unmapped_second.sam", new Object[][] { {"30,0", 1}, {"20,20", 1}}},
                {"unmapped_first.sam", "unmapped_second.sam", new Object[][] { {"0,20", 1}, {"30,0", 1}}},
                {"unmapped_first.sam", "unmapped_first.sam", new Object[][] { {"0,0", 1}, {"30,30", 1}}},
                {"genomic_sorted.sam", "genomic_sorted_sam_v1.6.sam", new Object[][] { {"20,20", 1}, {"30,30", 1}}},
                {"unsorted.sam", "unsorted.sam", new Object[][] { {"20,20", 1}, {"30,30", 1}}},
                {"unsorted.sam", "unsorted2.sam", new Object[][] { {"20,20", 1}}},
                {"duplicate_base.sam", "duplicate_four_mismatch_strict.sam", new Object[][] { {"20,20", 2}, {"30,30", 12}}},
                {"duplicate_base.sam", "duplicate_four_mismatch_lenient_one_align_differ.sam", new Object[][] { {"20,20", 2}, {"30,30", 12}}},
                {"duplicate_base.sam", "duplicate_two_mismatch_lenient.sam", new Object[][] { {"20,20", 2}, {"30,30", 12}}},
                {"duplicate_base.sam", "duplicate_four_mismatch_lenient.sam", new Object[][] { {"20,20", 2}, {"30,30", 12}}},
                {"duplicate_base.sam", "duplicate_four_mismatch_strict.sam", new Object[][] { {"20,20", 2}, {"30,30", 12}}},
                {"duplicate_base_queryname.sam", "duplicate_four_mismatch_strict_queryname.sam", new Object[][] { {"20,20", 2}, {"30,30", 12}}},
                {"duplicate_base_queryname.sam", "duplicate_four_mismatch_lenient_one_align_differ_queryname.sam", new Object[][] { {"20,20", 2}, {"30,30", 12}}},
                {"duplicate_base_queryname.sam", "duplicate_two_mismatch_lenient_queryname.sam", new Object[][] { {"20,20", 2}, {"30,30", 12}}},
                {"duplicate_base_queryname.sam", "duplicate_four_mismatch_lenient_queryname.sam", new Object[][] { {"20,20", 2}, {"30,30", 12}}},
                {"duplicate_base_queryname.sam", "duplicate_four_mismatch_strict_queryname.sam", new Object[][] { {"20,20", 2}, {"30,30", 12}}},
                {"genomic_sorted.sam", "mq0_2.sam", new Object[][] { {"20,0", 1}, {"30,30", 1}}},
                {"mq0_1.sam", "mq0_2.sam", new Object[][] { {"0,0", 1}, {"30,30", 1}}}
        };
    }

    @Test(dataProvider = "compareSamsMQConcordance")
    public void testMQConcordance(final String f1, final String f2, final Object[][] expectedMQConcordance) throws IOException {
        final File tmpOutput = File.createTempFile("compareSam", ".tsv");
        tmpOutput.deleteOnExit();
        final String in1 = new File(TEST_FILES_DIR, f1).getAbsolutePath();
        final String in2 = new File(TEST_FILES_DIR, f2).getAbsolutePath();
        final ArrayList<String> commandArgs = new ArrayList<>(
                Arrays.asList(
                        in1,
                        in2,
                        "O=" + tmpOutput,
                        "COMPARE_MQ=true"
                )
        );

        runPicardCommandLine(commandArgs);

        final MetricsFile<SamComparisonMetric, String> metricsOutput = new MetricsFile<>();
        metricsOutput.read(new FileReader(tmpOutput));

        // If the files cannot be compared (e.g. if their sort order differs) then there should be no histogram in the metrics file.
        if (expectedMQConcordance == null) {
            Assert.assertEquals(metricsOutput.getNumHistograms(), 0);
        } else {
            Assert.assertEquals(metricsOutput.getNumHistograms(), 1);

            final Map<String, Integer> expectedMQConcordanceMap = Stream.of(expectedMQConcordance).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1]));
            final Histogram<String> expectedMQConcordanceHistogram = new Histogram<>();
            expectedMQConcordanceMap.forEach(expectedMQConcordanceHistogram::increment);

            final Histogram<String> mqConcordanceHistogram = metricsOutput.getHistogram();
            Assert.assertEquals(mqConcordanceHistogram, expectedMQConcordanceHistogram);
        }
    }
}
