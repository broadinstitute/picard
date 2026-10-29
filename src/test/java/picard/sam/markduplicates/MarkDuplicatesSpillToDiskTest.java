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
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.DuplicationMetrics;
import picard.sam.markduplicates.MarkDuplicatesEndToEndRunner.Result;
import picard.sam.markduplicates.util.ReadEndsForMarkDuplicates;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import static picard.sam.markduplicates.MarkDuplicatesEndToEndRunner.assertMetrics;

/**
 * Tests that {@link MarkDuplicates} gives the same answer on coordinate-sorted input when its read ends spill to
 * disk: the pair and fragment sorting collections (a tiny SORTING_COLLECTION_SIZE_RATIO), the map of read ends
 * waiting for their mates on other contigs (MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1). The input has duplicate sets of
 * one to four pairs over three contigs and three read groups in two libraries, mates on other contigs, optical
 * duplicates, fragments alone and at pair ends, and unmapped reads; the expected flags, tags and metrics are worked
 * out as the input is built.
 */
public class MarkDuplicatesSpillToDiskTest {
    private static final boolean F = false, R = true;
    private static final String[] READ_GROUPS = {"rgA", "rgB", "rgC"};
    private static final Map<String, String> LIBRARY = new HashMap<>();

    static {
        LIBRARY.put("rgA", "libA");
        LIBRARY.put("rgB", "libB");
        LIBRARY.put("rgC", "libA");
    }

    /** The input and its expected outcome. */
    private static final class Expected {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner();
        final Set<String> duplicates = new TreeSet<>();
        final Set<String> optical = new TreeSet<>();
        // per library: unpaired, pairs, unmapped, unpaired duplicates, pair duplicates, optical pair duplicates
        final Map<String, long[]> counts = new TreeMap<>();
        final Map<Double, Double> allSets = new TreeMap<>();
        final Map<Double, Double> opticalSets = new TreeMap<>();
        final Map<Double, Double> nonOpticalSets = new TreeMap<>();

        long[] counts(final String readGroup) {
            return counts.computeIfAbsent(LIBRARY.get(readGroup), l -> new long[6]);
        }

        void addBothEnds(final Set<String> keys, final String name) {
            keys.add(name + "/1");
            keys.add(name + "/2");
        }
    }

    private static String illuminaName(final String prefix, final int x) {
        return prefix + ":1:1101:" + x + ":1000";
    }

    private static Expected buildInput() {
        final Expected expected = new Expected();
        final MarkDuplicatesEndToEndRunner input = expected.input;
        for (final String readGroup : READ_GROUPS) {
            input.readGroup(readGroup, LIBRARY.get(readGroup));
        }

        // duplicate sets of 1 to 4 pairs; the mate of each set is on the same or a later or an earlier contig
        for (int k = 0; k < 12; k++) {
            final String readGroup = READ_GROUPS[k % 3];
            final int contig1 = k % 3;
            final int contig2 = (k / 3) % 3;
            final int start1 = 1000 + 500 * k;
            final int start2 = 20000 + 300 * (11 - k);
            final int size = 1 + k % 4;
            int optical = 0;
            for (int j = 0; j < size; j++) {
                // odd members are within the optical distance of the keeper, even ones are far from everything
                final String name = illuminaName("s" + k + "p" + j, j % 2 == 1 ? 1000 + 10 * j : 1000 + 1000 * j);
                input.pair(name, readGroup, contig1, start1, F, "10M", contig2, start2, R, "10M", j == 0 ? 30 : 20);
                if (j > 0) {
                    expected.addBothEnds(expected.duplicates, name);
                    if (j % 2 == 1) {
                        expected.addBothEnds(expected.optical, name);
                        optical++;
                    }
                }
            }
            final long[] counts = expected.counts(readGroup);
            counts[1] += size;
            counts[4] += size - 1;
            counts[5] += optical;
            expected.allSets.merge((double) size, 1.0, Double::sum);
            if (optical > 0) {
                expected.opticalSets.merge((double) optical + 1, 1.0, Double::sum);
            }
            expected.nonOpticalSets.merge((double) size - optical, 1.0, Double::sum);

            // a fragment at the first end of every other set
            if (k % 2 == 1) {
                final String name = "atPair" + k;
                input.fragment(name, readGroup, contig1, start1, F, "10M", 40);
                expected.duplicates.add(name);
                counts[0]++;
                counts[3]++;
            }
        }

        // sets of 1 to 3 fragments on either strand
        for (int k = 0; k < 6; k++) {
            final String readGroup = READ_GROUPS[(k + 1) % 3];
            final long[] counts = expected.counts(readGroup);
            for (int j = 0; j <= k % 3; j++) {
                final String name = "f" + k + "m" + j;
                input.fragment(name, readGroup, k % 3, 50000 + 100 * k, k % 2 == 1, "10M", j == 0 ? 30 : 20);
                counts[0]++;
                if (j > 0) {
                    expected.duplicates.add(name);
                    counts[3]++;
                }
            }
        }

        // reads with unmapped mates, each alone at its position, and reads with no position
        for (int k = 0; k < 3; k++) {
            input.pairWithUnmappedMate("mateUnmapped" + k, READ_GROUPS[k], k, 70000, R, "10M", 30);
            input.unplacedFragment("unplaced" + k, READ_GROUPS[k], 30);
            final long[] counts = expected.counts(READ_GROUPS[k]);
            counts[0]++;
            counts[2] += 2;
        }
        return expected;
    }

    private static void assertExpected(final Expected expected, final Result result) {
        Assert.assertEquals(result.duplicates(), expected.duplicates);
        final Set<String> sequencing = result.records.stream()
                .filter(r -> MarkDuplicates.DUPLICATE_TYPE_SEQUENCING.equals(r.getAttribute(MarkDuplicates.DUPLICATE_TYPE_TAG)))
                .map(MarkDuplicatesEndToEndRunner::key).collect(Collectors.toCollection(TreeSet::new));
        Assert.assertEquals(sequencing, expected.optical);
        Assert.assertEquals(result.metricsFile.getMetrics().size(), expected.counts.size());
        for (final Map.Entry<String, long[]> entry : expected.counts.entrySet()) {
            final long[] c = entry.getValue();
            assertMetrics(result.metrics(entry.getKey()), entry.getKey(), c[0], c[1], 0, c[2], c[3], c[4], c[5]);
        }
        Assert.assertEquals(result.nonZeroBins("all_sets"), expected.allSets);
        Assert.assertEquals(result.nonZeroBins("optical_sets"), expected.opticalSets);
        Assert.assertEquals(result.nonZeroBins("non_optical_sets"), expected.nonOpticalSets);
    }

    /** A SORTING_COLLECTION_SIZE_RATIO that leaves room for two read ends in memory in this JVM. */
    private static String tinySortingCollectionRatio() {
        return Double.toString(2.5 * ReadEndsForMarkDuplicates.getSizeOf() / Runtime.getRuntime().maxMemory());
    }

    @DataProvider(name = "spillArguments")
    public Object[][] spillArguments() {
        return new Object[][]{
                {new String[]{"MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1"}},
                {new String[]{"SORTING_COLLECTION_SIZE_RATIO=" + tinySortingCollectionRatio()}},
                {new String[]{"SORTING_COLLECTION_SIZE_RATIO=" + tinySortingCollectionRatio(),
                        "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1"}},
        };
    }

    @Test(dataProvider = "spillArguments")
    public void testSpillingGivesTheSameAnswer(final String[] spillArguments) {
        final Expected expected = buildInput();
        final Result inMemory = expected.input.run("TAGGING_POLICY=All");
        assertExpected(expected, inMemory);

        final List<String> args = new ArrayList<>(Arrays.asList(spillArguments));
        args.add("TAGGING_POLICY=All");
        final Result spilled = expected.input.run(args.toArray(new String[0]));
        assertExpected(expected, spilled);

        Assert.assertEquals(spilled.records.stream().map(SAMRecord::getSAMString).collect(Collectors.toList()),
                inMemory.records.stream().map(SAMRecord::getSAMString).collect(Collectors.toList()));
        final List<DuplicationMetrics> spilledMetrics = spilled.metricsFile.getMetrics();
        final List<DuplicationMetrics> inMemoryMetrics = inMemory.metricsFile.getMetrics();
        for (int i = 0; i < inMemoryMetrics.size(); i++) {
            Assert.assertEquals(spilledMetrics.get(i).toString(), inMemoryMetrics.get(i).toString());
        }
    }
}
