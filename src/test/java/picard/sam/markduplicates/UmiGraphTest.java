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

import htsjdk.samtools.DuplicateSet;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import static picard.sam.markduplicates.MateCigarAndUmiTestSupport.*;

/**
 * Tests of {@link UmiGraph}: the connected components of UMIs within an edit distance, the UMI assigned to each
 * component (most common, avoiding Ns, fewest Ns, ties), the representative of each new duplicate set, the molecular
 * identifier and inferred-UMI attributes it sets, duplex UMIs, and missing UMIs.
 */
public class UmiGraphTest {

    /**
     * The first ends of FR pairs u0, u1, ... at the same position (first end at 1000, mate at 1200), with the given
     * UMIs in RX; a null UMI leaves the tag off. The names make u0 the best record, then u1, and so on.
     */
    private static DuplicateSet duplicateSet(final List<String> umis) {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        final DuplicateSet set = new DuplicateSet();
        for (int i = 0; i < umis.size(); i++) {
            final SAMRecord first = addPair(builder, "u" + i, 0, 1000, 1200).get(0);
            first.setAttribute("RX", umis.get(i));
            Assert.assertEquals(set.add(first), 0);
        }
        return set;
    }

    private static Set<Set<String>> namesBySet(final List<DuplicateSet> sets) {
        return sets.stream()
                .map(set -> set.getRecords().stream().map(SAMRecord::getReadName).collect(Collectors.toSet()))
                .collect(Collectors.toSet());
    }

    private static Set<Set<String>> partition(final String... groups) {
        final Set<Set<String>> partition = new HashSet<>();
        for (final String group : groups) {
            final Set<String> names = new HashSet<>();
            for (final String index : group.split(",")) {
                names.add("u" + index);
            }
            partition.add(names);
        }
        return partition;
    }

    private static String inferredUmi(final SAMRecord record) {
        return (String) record.getTransientAttribute(UmiUtil.INFERRED_UMI_TRANSIENT_TAG);
    }

    @DataProvider(name = "components")
    public Object[][] components() {
        return new Object[][]{
                // UMIs, maximum edit distance, the sets of pair indices
                {List.of("ACGT"), 0, partition("0")},
                {List.of("ACGT", "ACGT"), 0, partition("0,1")},
                {List.of("AAAA", "CCCC", "GGGG"), 0, partition("0", "1", "2")},
                {List.of("AAAA", "AAAT", "AATT"), 0, partition("0", "1", "2")},
                {List.of("AAAA", "AAAT", "AATT"), 1, partition("0,1,2")},
                {List.of("AAAA", "AATT"), 1, partition("0", "1")},
                {List.of("AAAA", "AATT"), 2, partition("0,1")},
                {List.of("AAAA", "ATTT"), 2, partition("0", "1")},
                {List.of("AAAA", "TTTT"), 4, partition("0,1")},
                // The example in the class Javadoc: ATCC-AACC and AACC-AACG are joined, so ATCC and AACG are too.
                {List.of("ATCC", "AACC", "AACG"), 1, partition("0,1,2")},
                {List.of("AAAA", "AAAT", "CCCC", "CCCG", "CCGG", "TTTT"), 1, partition("0,1", "2,3,4", "5")},
                {List.of("ANAA", "ANAA"), 0, partition("0,1")},
                {List.of("ANAA", "AAAA"), 0, partition("0", "1")},
                {List.of("ANAA", "AANA"), 1, partition("0", "1")},
                {List.of("ANAA", "AANA"), 2, partition("0,1")},
        };
    }

    @Test(dataProvider = "components")
    public void testConnectedComponents(final List<String> umis, final int maxEditDistance, final Set<Set<String>> expected) {
        final List<DuplicateSet> sets = new UmiGraph(duplicateSet(umis), "RX", "MI", false, false)
                .joinUmisIntoDuplicateSets(maxEditDistance);
        Assert.assertEquals(namesBySet(sets), expected);
    }

    @DataProvider(name = "assignedUmis")
    public Object[][] assignedUmis() {
        return new Object[][]{
                // UMIs of one component, maximum edit distance, the UMI assigned to all of them
                {List.of("ATCC", "AACC", "AACC", "AACG"), 1, "AACC"},
                {List.of("AAAT", "AAAA", "AAAA"), 1, "AAAA"},
                // A tie goes to the UMI of the best record.
                {List.of("AAAT", "AAAA"), 1, "AAAT"},
                {List.of("AAAA", "AAAT"), 1, "AAAA"},
                // A UMI without N is preferred even when it is the least common.
                {List.of("AANA", "AANA", "AANA", "AAAA"), 1, "AAAA"},
                // With an N in every UMI, the one with the fewest Ns, then the earliest of those.
                {List.of("NNAA", "NAAA", "NNAA"), 1, "NAAA"},
                {List.of("NNAA", "ANNA"), 4, "NNAA"},
                {List.of("NNNA", "NNAA", "NAAN"), 4, "NNAA"},
        };
    }

    /**
     * A tie between equally common UMIs, or between UMIs with equally few Ns, goes to the UMI of the record that comes
     * first in the duplicate set. The documentation does not specify this; the expected value follows the current
     * implementation.
     */
    @Test(dataProvider = "assignedUmis")
    public void testAssignedUmi(final List<String> umis, final int maxEditDistance, final String expectedUmi) {
        final List<DuplicateSet> sets = new UmiGraph(duplicateSet(umis), "RX", "MI", false, false)
                .joinUmisIntoDuplicateSets(maxEditDistance);
        Assert.assertEquals(sets.size(), 1);
        final List<SAMRecord> records = sets.get(0).getRecords();
        Assert.assertEquals(records.size(), umis.size());
        for (final SAMRecord record : records) {
            Assert.assertEquals(inferredUmi(record), expectedUmi);
            Assert.assertEquals(record.getAttribute("MI"), "chr1:1200/" + expectedUmi);
            Assert.assertEquals(record.getAttribute("RX"), umis.get(Integer.parseInt(record.getReadName().substring(1))));
        }
    }

    @Test
    public void testEachSetHasItsOwnUmiAndRepresentative() {
        final List<DuplicateSet> sets = new UmiGraph(duplicateSet(List.of("CCCC", "AAAA", "CCCC", "AAAT", "AAAT", "AAAT")), "RX", "MI", false, false)
                .joinUmisIntoDuplicateSets(1);
        Assert.assertEquals(sets.size(), 2);
        for (final DuplicateSet set : sets) {
            final String representative = set.getRepresentative().getReadName();
            final Set<String> names = set.getRecords().stream().map(SAMRecord::getReadName).collect(Collectors.toSet());
            if (names.contains("u0")) {
                Assert.assertEquals(names, Set.of("u0", "u2"));
                Assert.assertEquals(representative, "u0");
                set.getRecords().forEach(r -> Assert.assertEquals(inferredUmi(r), "CCCC"));
            } else {
                Assert.assertEquals(names, Set.of("u1", "u3", "u4", "u5"));
                Assert.assertEquals(representative, "u1");
                set.getRecords().forEach(r -> Assert.assertEquals(inferredUmi(r), "AAAT"));
            }
            // The representative of each set is the only record not flagged as a duplicate.
            set.getRecords().forEach(r -> Assert.assertEquals(r.getDuplicateReadFlag(), !r.getReadName().equals(representative)));
        }
    }

    @Test
    public void testEmptySet() {
        Assert.assertTrue(new UmiGraph(new DuplicateSet(), "RX", "MI", false, false).joinUmisIntoDuplicateSets(1).isEmpty());
    }

    /**
     * Without a molecular identifier tag no MI is written. AAAA and AAAT, seen once each, both get AAAA, the UMI of
     * the first record. The documentation does not specify this; the expected value follows the current
     * implementation.
     */
    @Test
    public void testNoMolecularIdentifierTag() {
        final List<DuplicateSet> sets = new UmiGraph(duplicateSet(List.of("AAAA", "AAAT")), "RX", null, false, false)
                .joinUmisIntoDuplicateSets(1);
        Assert.assertEquals(sets.size(), 1);
        for (final SAMRecord record : sets.get(0).getRecords()) {
            Assert.assertNull(record.getAttribute("MI"));
            Assert.assertEquals(inferredUmi(record), "AAAA");
        }
    }

    @Test
    public void testOtherUmiTag() {
        final DuplicateSet set = duplicateSet(List.of("AAAA", "CCCC"));
        set.getRecords().forEach(r -> {
            r.setAttribute("ZU", r.getStringAttribute("RX"));
            r.setAttribute("RX", null);
        });
        final List<DuplicateSet> sets = new UmiGraph(set, "ZU", "XM", false, false).joinUmisIntoDuplicateSets(0);
        Assert.assertEquals(namesBySet(sets), partition("0", "1"));
        for (final DuplicateSet result : sets) {
            final SAMRecord record = result.getRecords().get(0);
            Assert.assertEquals(record.getAttribute("XM"), "chr1:1200/" + record.getStringAttribute("ZU"));
            Assert.assertNull(record.getAttribute("MI"));
        }
    }

    @Test
    public void testMissingUmiRejected() {
        final List<String> umis = new ArrayList<>();
        umis.add("AAAA");
        umis.add(null);
        final PicardException exception = Assert.expectThrows(PicardException.class,
                () -> new UmiGraph(duplicateSet(umis), "RX", "MI", false, false));
        Assert.assertTrue(exception.getMessage().contains("u1"), exception.getMessage());
    }

    /**
     * When missing UMIs are allowed, records without one share the empty UMI and get no molecular identifier tag. The
     * documentation does not specify this; the expected value follows the current implementation.
     */
    @Test
    public void testMissingUmisAllowed() {
        final List<String> umis = new ArrayList<>();
        umis.add(null);
        umis.add(null);
        final List<DuplicateSet> sets = new UmiGraph(duplicateSet(umis), "RX", "MI", true, false).joinUmisIntoDuplicateSets(0);
        Assert.assertEquals(namesBySet(sets), partition("0,1"));
        for (final SAMRecord record : sets.get(0).getRecords()) {
            Assert.assertNull(record.getAttribute("RX"));
            Assert.assertNull(record.getAttribute("MI"));
            Assert.assertEquals(inferredUmi(record), "");
        }
    }

    @Test
    public void testDuplexUmisAreJoinedAcrossStrands() {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        // Both records start at 1000 with the mate's 5' end at 1249: d1's is the first end of a top-strand pair, d2's
        // the second end of a bottom-strand pair, whose UMI halves are swapped.
        final SAMRecord top = addPair(builder, "d1", 0, 1000, false, "50M", 1200, true, "50M").get(0);
        final SAMRecord bottom = addPair(builder, "d2", 0, 1200, true, "50M", 1000, false, "50M").get(1);
        top.setAttribute("RX", "AAA-CCC");
        bottom.setAttribute("RX", "CCC-AAA");
        final DuplicateSet set = new DuplicateSet();
        Assert.assertEquals(set.add(top), 0);
        Assert.assertEquals(set.add(bottom), 0);

        final List<DuplicateSet> duplex = new UmiGraph(set, "RX", "MI", false, true).joinUmisIntoDuplicateSets(0);
        Assert.assertEquals(duplex.size(), 1);
        Assert.assertEquals(inferredUmi(top), "AAA-CCC");
        Assert.assertEquals(inferredUmi(bottom), "AAA-CCC");
        Assert.assertEquals(top.getAttribute("MI"), "chr1:1200/AAA-CCC/A");
        Assert.assertEquals(bottom.getAttribute("MI"), "chr1:1200/AAA-CCC/B");
        Assert.assertEquals(bottom.getAttribute("RX"), "CCC-AAA");

        final List<DuplicateSet> single = new UmiGraph(set, "RX", "MI", false, false).joinUmisIntoDuplicateSets(5);
        Assert.assertEquals(single.size(), 2);
        Assert.assertEquals(top.getAttribute("MI"), "chr1:1200/AAA-CCC");
        Assert.assertEquals(bottom.getAttribute("MI"), "chr1:1200/CCC-AAA");
    }
}
