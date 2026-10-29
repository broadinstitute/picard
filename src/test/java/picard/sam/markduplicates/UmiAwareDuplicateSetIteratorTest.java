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

import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;
import htsjdk.samtools.DuplicateSet;
import htsjdk.samtools.DuplicateSetIterator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordDuplicateComparator;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.util.CloseableIterator;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import static picard.sam.markduplicates.MateCigarAndUmiTestSupport.*;

/**
 * Tests of {@link UmiAwareDuplicateSetIterator} over an htsjdk {@link DuplicateSetIterator}: the sets it splits each
 * duplicate set into, the end of iteration, and the UMI metrics it accumulates per library and finishes on close.
 */
public class UmiAwareDuplicateSetIteratorTest {

    private static DuplicateSetIterator duplicateSetIterator(final SAMRecordSetBuilder builder) {
        final SAMFileHeader header = builder.getHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final List<SAMRecord> records = new ArrayList<>(builder.getRecords());
        records.sort(new SAMRecordCoordinateComparator());
        final Iterator<SAMRecord> iterator = records.iterator();
        final CloseableIterator<SAMRecord> closeable = new CloseableIterator<SAMRecord>() {
            @Override
            public void close() {
            }

            @Override
            public boolean hasNext() {
                return iterator.hasNext();
            }

            @Override
            public SAMRecord next() {
                return iterator.next();
            }
        };
        final SAMRecordDuplicateComparator comparator = new SAMRecordDuplicateComparator(List.of(header));
        comparator.setScoringStrategy(ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH);
        return new DuplicateSetIterator(closeable, header, false, comparator);
    }

    @Test
    public void testSplitsSetsByUmiAndFinishesMetricsOnClose() {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        final String[] umis = {"AAAA", "AAAA", "CCCC", "AAAT"};
        for (int i = 0; i < umis.length; i++) {
            final String umi = umis[i];
            addPair(builder, "u" + i, 0, 1000, 1200).forEach(r -> r.setAttribute("RX", umi));
        }
        addPair(builder, "v0", 0, 3000, 3200).forEach(r -> r.setAttribute("RX", "GGGG"));

        final Map<String, UmiMetrics> metricsByLibrary = new HashMap<>();
        final UmiAwareDuplicateSetIterator iterator = new UmiAwareDuplicateSetIterator(duplicateSetIterator(builder), 1, "RX", "MI",
                false, false, metricsByLibrary);
        final Set<Set<String>> sets = new HashSet<>();
        while (iterator.hasNext()) {
            final DuplicateSet set = iterator.next();
            sets.add(set.getRecords().stream().map(MateCigarAndUmiTestSupport::key).collect(Collectors.toSet()));
        }
        Assert.assertFalse(iterator.hasNext());
        Assert.assertEquals(sets, Set.of(
                Set.of("u0/1", "u1/1", "u3/1"), Set.of("u2/1"), Set.of("u0/2", "u1/2", "u3/2"), Set.of("u2/2"),
                Set.of("v0/1"), Set.of("v0/2")));

        iterator.close();
        Assert.assertFalse(iterator.hasNext());
        Assert.assertEquals(metricsByLibrary.keySet(), Set.of("libA"));
        final UmiMetrics metrics = metricsByLibrary.get("libA");
        Assert.assertEquals(metrics.LIBRARY, "libA");
        Assert.assertEquals(metrics.MEAN_UMI_LENGTH, 4.0);
        Assert.assertEquals(metrics.DUPLICATE_SETS_IGNORING_UMI, 4);
        Assert.assertEquals(metrics.DUPLICATE_SETS_WITH_UMI, 6);
        // Ten reads: AAAA x4, CCCC x2, AAAT x2 (both corrected to AAAA) and GGGG x2; 2 errors in 40 bases is Q13.
        Assert.assertEquals(metrics.OBSERVED_UNIQUE_UMIS, 4);
        Assert.assertEquals(metrics.INFERRED_UNIQUE_UMIS, 3);
        Assert.assertEquals(metrics.OBSERVED_BASE_ERRORS, 2);
        Assert.assertEquals(metrics.UMI_BASE_QUALITIES, 13.0);
        Assert.assertEquals(metrics.PCT_UMI_WITH_N, 0.0);
        final double log4 = Math.log(4);
        Assert.assertEquals(metrics.OBSERVED_UMI_ENTROPY,
                -(0.4 * Math.log(0.4) + 3 * 0.2 * Math.log(0.2)) / log4, 1e-9);
        Assert.assertEquals(metrics.INFERRED_UMI_ENTROPY,
                -(0.6 * Math.log(0.6) + 2 * 0.2 * Math.log(0.2)) / log4, 1e-9);
    }

    @Test
    public void testEmptyInput() {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        final Map<String, UmiMetrics> metricsByLibrary = new HashMap<>();
        final UmiAwareDuplicateSetIterator iterator = new UmiAwareDuplicateSetIterator(duplicateSetIterator(builder), 1, "RX", "MI",
                false, false, metricsByLibrary);
        Assert.assertFalse(iterator.hasNext());
        iterator.close();
        Assert.assertTrue(metricsByLibrary.isEmpty());
    }
}
