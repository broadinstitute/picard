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
import htsjdk.samtools.SAMTag;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.DuplicationMetrics;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import static picard.sam.markduplicates.MateCigarAndUmiTestSupport.*;

/**
 * Tests of {@link CollectDuplicateMetrics}: the per-library duplication metrics it computes from the duplicate flags
 * already in a file, which records count as examined, unmapped, secondary or duplicate, the libraries it reports, and
 * the optical duplicate count it never estimates.
 */
public class CollectDuplicateMetricsBehaviorTest {

    private static Map<String, DuplicationMetrics> run(final File input) throws IOException {
        final File metrics = newTempFile("CollectDuplicateMetricsBehaviorTest", ".metrics");
        final String[] args = {"INPUT=" + input.getAbsolutePath(), "METRICS_FILE=" + metrics.getAbsolutePath()};
        Assert.assertEquals(new CollectDuplicateMetrics().instanceMain(args), 0);
        return readDuplicationMetrics(metrics);
    }

    private static void flag(final List<SAMRecord> records) {
        records.forEach(r -> r.setDuplicateReadFlag(true));
    }

    /**
     * Two libraries with duplicate flags already set. Flags on unmapped, secondary and supplementary records are ignored.
     * The pairs in libA are named as optical duplicates of each other.
     */
    private static SAMRecordSetBuilder markedInput() {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        addPair(builder, "RUN:1:1101:1000:1000", 0, 1000, 1200);
        flag(addPair(builder, "RUN:1:1101:1001:1000", 0, 1000, 1200));
        addFragment(builder, "f1", 0, 3000, false, "50M");
        addFragment(builder, "f2", 0, 3000, false, "50M").setDuplicateReadFlag(true);
        flag(addPairWithUnmappedMate(builder, "u1", 0, 5000));
        addSecondaryFragment(builder, "s1", 0, 7000).setDuplicateReadFlag(true);
        addSupplementaryFragment(builder, "s2", 0, 7000);
        addUnmappedPair(builder, "z1").get(0).setDuplicateReadFlag(true);

        useReadGroup(builder, "rgB", "libB");
        addPair(builder, "b1", 1, 1000, 1200);
        flag(addPair(builder, "b2", 1, 1000, 1200));
        flag(addPair(builder, "b3", 1, 1000, 1200));
        addFragment(builder, "b4", 1, 3000, true, "50M");
        return builder;
    }

    @DataProvider(name = "sortOrders")
    public Object[][] sortOrders() {
        return new Object[][]{{true}, {false}};
    }

    @Test(dataProvider = "sortOrders")
    public void testMetricsFromMarkedFile(final boolean coordinateSorted) throws IOException {
        final SAMRecordSetBuilder builder = markedInput();
        final File input = coordinateSorted ? writeCoordinateSorted(builder) : writeQueryNameSorted(builder.getHeader(), builder.getRecords());

        final Map<String, DuplicationMetrics> metrics = run(input);
        Assert.assertEquals(metrics.size(), 2);
        // libA: f1, f2 and u1's mapped end unpaired; u1's unmapped end and both ends of z1 unmapped.
        assertDuplicationMetrics(metrics.get("libA"), "libA", 3, 2, 2, 3, 2, 1, 0);
        assertDuplicationMetrics(metrics.get("libB"), "libB", 1, 3, 0, 0, 0, 2, 0);
    }

    @Test
    public void testReadGroupsWithoutLibraryOrReads() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgEmpty", "libEmpty");
        useReadGroup(builder, "rgNoLibrary", null);
        addPair(builder, "p1", 0, 1000, 1200);
        flag(addPair(builder, "p2", 0, 1000, 1200));

        final Map<String, DuplicationMetrics> metrics = run(writeCoordinateSorted(builder));
        Assert.assertEquals(metrics.size(), 2);
        assertDuplicationMetrics(metrics.get("Unknown Library"), "Unknown Library", 0, 2, 0, 0, 0, 1, 0);
        assertDuplicationMetrics(metrics.get("libEmpty"), "libEmpty", 0, 0, 0, 0, 0, 0, 0);
    }

    @Test
    public void testRecordsWithoutReadGroup() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        final List<SAMRecord> records = new ArrayList<>();
        records.add(addFragment(builder, "f1", 0, 1000, false, "50M"));
        records.add(addFragment(builder, "f2", 0, 1000, false, "50M"));
        records.add(addFragment(builder, "f3", 0, 1000, false, "50M"));
        records.forEach(r -> r.setAttribute(SAMTag.RG.name(), null));
        records.get(1).setDuplicateReadFlag(true);
        records.get(2).setDuplicateReadFlag(true);

        final Map<String, DuplicationMetrics> metrics = run(writeCoordinateSorted(builder.getHeader(), records));
        Assert.assertEquals(metrics.size(), 1);
        assertDuplicationMetrics(metrics.get("Unknown Library"), "Unknown Library", 3, 0, 0, 0, 2, 0, 0);
    }

    @Test
    public void testEmptyInput() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        final Map<String, DuplicationMetrics> metrics = run(writeCoordinateSorted(builder));
        Assert.assertEquals(metrics.size(), 1);
        assertDuplicationMetrics(metrics.get("libA"), "libA", 0, 0, 0, 0, 0, 0, 0);
    }
}
