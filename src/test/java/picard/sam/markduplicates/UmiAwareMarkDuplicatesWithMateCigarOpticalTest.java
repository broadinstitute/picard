/*
 * The MIT License
 *
 * Copyright (c) 2026 The Broad Institute
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

package picard.sam.markduplicates;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import static picard.sam.markduplicates.MateCigarAndUmiTestSupport.*;

/**
 * Tests optical duplicate detection in {@link UmiAwareMarkDuplicatesWithMateCigar}: optical duplicates are counted
 * within the duplicate sets that remain after splitting by UMI, against the pair kept in each. All reads are in one
 * read group of library libA, and all pairs have the first end forward at chr1:1000 and the second end reverse at
 * chr1:1200 (50M); with equal scores the smallest name is kept.
 *
 * Read names give tile 1101 and x positions k1 1000, k2 1020, o1 1050, o2 1080 and f1 3000, all at y 1000. k1, o1 and f1
 * have UMI AAAA; k2 and o2 have UMI AAAT.
 */
public class UmiAwareMarkDuplicatesWithMateCigarOpticalTest {
    private static final String K1 = "RUN:1:1101:1000:1000";
    private static final String K2 = "RUN:1:1101:1020:1000";
    private static final String O1 = "RUN:1:1101:1050:1000";
    private static final String O2 = "RUN:1:1101:1080:1000";
    private static final String F1 = "RUN:1:1101:3000:1000";

    private static File input() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        for (final String name : List.of(K1, O1, F1)) {
            addPair(builder, name, 0, 1000, 1200).forEach(r -> r.setAttribute("RX", "AAAA"));
        }
        for (final String name : List.of(K2, O2)) {
            addPair(builder, name, 0, 1000, 1200).forEach(r -> r.setAttribute("RX", "AAAT"));
        }
        return writeCoordinateSorted(builder);
    }

    @DataProvider(name = "editDistances")
    public Object[][] editDistances() {
        return new Object[][]{
                // With 0 the UMIs form two sets, {k1, o1, f1} keeping k1 and {k2, o2} keeping k2: o1 and o2 are optical
                // duplicates, although k2 lies within 100 pixels of k1 and o1.
                {0, List.of(O1, F1, O2), 2, Map.of(3.0, 1.0, 2.0, 1.0), Map.of(2.0, 2.0), Map.of(2.0, 1.0, 1.0, 1.0)},
                // With 1 they form one set keeping k1, in which k2, o1 and o2 are optical duplicates.
                {1, List.of(K2, O1, O2, F1), 3, Map.of(5.0, 1.0), Map.of(4.0, 1.0), Map.of(2.0, 1.0)},
        };
    }

    @Test(dataProvider = "editDistances")
    public void testOpticalDuplicatesWithinUmiSets(final int maxEditDistance, final List<String> duplicatePairs,
                                                  final int opticalDuplicates, final Map<Double, Double> allSets,
                                                  final Map<Double, Double> opticalSets,
                                                  final Map<Double, Double> nonOpticalSets) throws IOException {
        final File output = newTempFile("UmiAwareMarkDuplicatesWithMateCigarOpticalTest", ".sam");
        final File metrics = newTempFile("UmiAwareMarkDuplicatesWithMateCigarOpticalTest", ".metrics");
        final File umiMetrics = newTempFile("UmiAwareMarkDuplicatesWithMateCigarOpticalTest", ".umi_metrics");
        final String[] args = {
                "INPUT=" + input().getAbsolutePath(),
                "OUTPUT=" + output.getAbsolutePath(),
                "METRICS_FILE=" + metrics.getAbsolutePath(),
                "UMI_METRICS_FILE=" + umiMetrics.getAbsolutePath(),
                "DUPLICATE_SCORING_STRATEGY=TOTAL_MAPPED_REFERENCE_LENGTH",
                "MAX_EDIT_DISTANCE_TO_JOIN=" + maxEditDistance};
        Assert.assertEquals(new UmiAwareMarkDuplicatesWithMateCigar().instanceMain(args), 0);

        final List<SAMRecord> records = readRecords(output);
        final List<String> expectedDuplicates = new ArrayList<>();
        duplicatePairs.forEach(name -> expectedDuplicates.addAll(List.of(name + "/1", name + "/2")));
        assertDuplicates(records, expectedDuplicates);
        Assert.assertEquals(records.size(), 10);

        assertDuplicationMetrics(readDuplicationMetrics(metrics).get("libA"), "libA", 0, 5, 0, 0, 0,
                duplicatePairs.size(), opticalDuplicates);
        final Map<String, Map<Double, Double>> histograms = readHistograms(metrics);
        Assert.assertEquals(histograms.get("all_sets"), allSets);
        Assert.assertEquals(histograms.get("optical_sets"), opticalSets);
        Assert.assertEquals(histograms.get("non_optical_sets"), nonOpticalSets);
    }
}
