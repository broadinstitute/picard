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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordSetBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.DuplicationMetrics;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import static picard.sam.markduplicates.MateCigarAndUmiTestSupport.*;

/**
 * End-to-end tests of {@link SimpleMarkDuplicatesWithMateCigar}: on the inputs used for
 * {@link MarkDuplicatesWithMateCigar} it finds the same duplicate sets, keepers and metrics, including optical
 * duplicates, separate libraries, removal of duplicates, and spilling its sorts to disk.
 */
public class SimpleMarkDuplicatesWithMateCigarBehaviorTest {

    private static class Result {
        final List<SAMRecord> records;
        final SAMFileHeader header;
        final Map<String, DuplicationMetrics> metrics;
        final Map<String, Map<Double, Double>> histograms;

        Result(final File output, final File metrics) throws IOException {
            this.records = readRecords(output);
            this.header = readHeader(output);
            this.metrics = readDuplicationMetrics(metrics);
            this.histograms = readHistograms(metrics);
        }
    }

    private static Result run(final File input, final String... extraArgs) throws IOException {
        final File output = newTempFile("SimpleMarkDuplicatesWithMateCigarBehaviorTest", ".sam");
        final File metrics = newTempFile("SimpleMarkDuplicatesWithMateCigarBehaviorTest", ".metrics");
        final List<String> args = new ArrayList<>(List.of(
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + output.getAbsolutePath(),
                "METRICS_FILE=" + metrics.getAbsolutePath(),
                "DUPLICATE_SCORING_STRATEGY=TOTAL_MAPPED_REFERENCE_LENGTH"));
        args.addAll(Arrays.asList(extraArgs));
        Assert.assertEquals(new SimpleMarkDuplicatesWithMateCigar().instanceMain(args.toArray(new String[0])), 0);
        return new Result(output, metrics);
    }

    @DataProvider(name = "settings")
    public Object[][] settings() {
        return new Object[][]{
                {false, 500000},
                {true, 500000},
                {false, 2},
                {true, 3},
        };
    }

    @Test(dataProvider = "settings")
    public void testMixedInputAgreesWithExpectedDuplicateSets(final boolean removeDuplicates, final int maxRecordsInRam) throws IOException {
        final File input = writeCoordinateSorted(mixedInput());
        final List<String> inputOrder = keys(readRecords(input));

        final Result result = run(input, "REMOVE_DUPLICATES=" + removeDuplicates, "MAX_RECORDS_IN_RAM=" + maxRecordsInRam);

        // Output is re-sorted, so records tied on position may swap once their duplicate flags change.
        Assert.assertEquals(result.header.getSortOrder(), SAMFileHeader.SortOrder.coordinate);
        final SAMRecordCoordinateComparator comparator = new SAMRecordCoordinateComparator();
        for (int i = 1; i < result.records.size(); i++) {
            Assert.assertTrue(comparator.compare(result.records.get(i - 1), result.records.get(i)) <= 0);
        }
        final List<String> expectedKeys = inputOrder.stream()
                .filter(key -> !removeDuplicates || !MIXED_DUPLICATES.contains(key))
                .collect(Collectors.toList());
        Assert.assertEqualsNoOrder(keys(result.records).toArray(), expectedKeys.toArray());
        assertDuplicates(result.records, removeDuplicates ? List.of() : MIXED_DUPLICATES);
        Assert.assertEquals(result.metrics.size(), 1);
        assertMixedMetrics(result.metrics.get(MIXED_LIBRARY));
    }

    @Test
    public void testFivePrimeEndsIncludeClipsOnBothEnds() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        addPair(builder, "cl1", 0, 1000, false, "50M", 1300, true, "50M");
        addPair(builder, "cl2", 0, 1005, false, "5S45M", 1300, true, "45M5S");
        addPair(builder, "cl3", 0, 1003, false, "3H47M", 1302, true, "2H45M3H");
        addPair(builder, "cl4", 0, 1000, false, "50M", 1300, true, "45M");
        addPair(builder, "cl5", 0, 1010, false, "5S45M", 1300, true, "50M");

        final Result result = run(writeCoordinateSorted(builder));
        assertDuplicates(result.records, List.of("cl2/1", "cl2/2", "cl3/1", "cl3/2"));
        assertDuplicationMetrics(result.metrics.get("libA"), "libA", 0, 5, 0, 0, 0, 2, 0);
    }

    @DataProvider(name = "opticalSettings")
    public Object[][] opticalSettings() {
        return new Object[][]{
                {new String[]{}, 1, Map.of(2.0, 1.0), Map.of(3.0, 1.0, 2.0, 1.0)},
                {new String[]{"OPTICAL_DUPLICATE_PIXEL_DISTANCE=99"}, 0, Map.of(), Map.of(4.0, 1.0, 2.0, 1.0)},
                {new String[]{"OPTICAL_DUPLICATE_PIXEL_DISTANCE=101"}, 3, Map.of(3.0, 1.0, 2.0, 1.0), Map.of(2.0, 1.0, 1.0, 1.0)},
        };
    }

    @Test(dataProvider = "opticalSettings")
    public void testOpticalDuplicates(final String[] args, final int expectedOpticalDuplicates,
                                      final Map<Double, Double> opticalSets, final Map<Double, Double> nonOpticalSets) throws IOException {
        final Result result = run(opticalInput(false), args);
        assertDuplicates(result.records, opticalDuplicates(false));
        assertDuplicationMetrics(result.metrics.get("libA"), "libA", 0, 6, 0, 0, 0, 4, expectedOpticalDuplicates);
        Assert.assertEquals(result.histograms.get("all_sets"), Map.of(4.0, 1.0, 2.0, 1.0));
        Assert.assertEquals(result.histograms.getOrDefault("optical_sets", Map.of()), opticalSets);
        Assert.assertEquals(result.histograms.get("non_optical_sets"), nonOpticalSets);
    }

    @Test
    public void testOpticalDuplicatesWithoutReadNameRegex() throws IOException {
        final Result result = run(opticalInput(true), "READ_NAME_REGEX=null");
        assertDuplicates(result.records, opticalDuplicates(true));
        assertDuplicationMetrics(result.metrics.get("libA"), "libA", 0, 7, 0, 0, 0, 5, 0);
    }

    @Test
    public void testLibrariesAreMarkedSeparately() throws IOException {
        final Result result = run(writeCoordinateSorted(multiLibraryInput()));
        assertMultiLibraryResult(result.records, result.metrics);
    }

    @Test
    public void testPairsWithMatesOnOtherContigs() throws IOException {
        final Result result = run(writeCoordinateSorted(interContigInput()));
        assertDuplicates(result.records, List.of("x2/1", "x2/2"));
        assertDuplicationMetrics(result.metrics.get("libA"), "libA", 0, 3, 0, 0, 0, 1, 0);
    }

    @Test
    public void testCommentsAreAddedToHeader() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        addPair(builder, "p1", 0, 1000, 1200);
        final Result result = run(writeCoordinateSorted(builder), "COMMENT=first note", "COMMENT=second");
        Assert.assertEquals(result.header.getComments(), List.of("@CO\tfirst note", "@CO\tsecond"));
    }

    @Test
    public void testEmptyInput() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");

        final Result result = run(writeCoordinateSorted(builder));
        Assert.assertTrue(result.records.isEmpty());
        Assert.assertEquals(result.metrics.size(), 1);
        assertDuplicationMetrics(result.metrics.get("libA"), "libA", 0, 0, 0, 0, 0, 0, 0);
    }
}
