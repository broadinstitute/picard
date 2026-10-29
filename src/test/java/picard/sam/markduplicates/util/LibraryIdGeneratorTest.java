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

package picard.sam.markduplicates.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Histogram;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.DuplicationMetrics;
import picard.sam.FlowBasedDuplicationMetrics;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;

/**
 * Tests for {@link LibraryIdGenerator}: library names for read groups and records (including records with no read
 * group, a read group missing from the header, and a read group with no LB), library ids shared by the read groups of
 * one library and distinct between libraries, the metrics-by-library map built from the header, and the histograms
 * it keeps.
 */
public class LibraryIdGeneratorTest {

    private static final String UNKNOWN = "Unknown Library";

    private static SAMFileHeader header(final String[]... readGroupsAndLibraries) {
        final SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr1", 1000));
        for (final String[] rgAndLibrary : readGroupsAndLibraries) {
            final SAMReadGroupRecord rg = new SAMReadGroupRecord(rgAndLibrary[0]);
            if (rgAndLibrary[1] != null) rg.setLibrary(rgAndLibrary[1]);
            rg.setSample("sample");
            header.addReadGroup(rg);
        }
        return header;
    }

    /** Read groups: rgA and rgB in libA, rgC in libC, rgNoLib without LB. */
    private static SAMFileHeader multiLibraryHeader() {
        return header(new String[]{"rgA", "libA"}, new String[]{"rgB", "libA"}, new String[]{"rgC", "libC"}, new String[]{"rgNoLib", null});
    }

    private static SAMRecord record(final SAMFileHeader header, final String readGroup) {
        final SAMRecord rec = new SAMRecord(header);
        rec.setReadName("r");
        if (readGroup != null) rec.setAttribute("RG", readGroup);
        return rec;
    }

    @DataProvider(name = "libraryNames")
    public Object[][] libraryNames() {
        return new Object[][]{
                {"rgA", "libA"},
                {"rgB", "libA"},
                {"rgC", "libC"},
                {"rgNoLib", UNKNOWN},
                {"rgNotInHeader", UNKNOWN},
                {null, UNKNOWN}
        };
    }

    @Test(dataProvider = "libraryNames")
    public void testGetLibraryName(final String readGroup, final String expected) {
        final SAMFileHeader header = multiLibraryHeader();
        Assert.assertEquals(LibraryIdGenerator.getLibraryName(header, record(header, readGroup)), expected);
    }

    @Test
    public void testGetReadGroupLibraryName() {
        final SAMReadGroupRecord withLibrary = new SAMReadGroupRecord("rg1");
        withLibrary.setLibrary("lib1");
        Assert.assertEquals(LibraryIdGenerator.getReadGroupLibraryName(withLibrary), "lib1");
        Assert.assertEquals(LibraryIdGenerator.getReadGroupLibraryName(new SAMReadGroupRecord("rg2")), UNKNOWN);
    }

    /**
     * Records of one library share an id whatever their read group, records of different libraries get different ids,
     * and records without a known library share the id of the unknown library.
     */
    @Test
    public void testLibraryIdsFollowTheLibrary() {
        final SAMFileHeader header = multiLibraryHeader();
        final LibraryIdGenerator generator = new LibraryIdGenerator(header);

        final short libC = generator.getLibraryId(record(header, "rgC"));
        final short unknown = generator.getLibraryId(record(header, null));
        final short libA = generator.getLibraryId(record(header, "rgB"));
        Assert.assertEquals(generator.getLibraryId(record(header, "rgA")), libA);
        Assert.assertEquals(generator.getLibraryId(record(header, "rgNoLib")), unknown);
        Assert.assertEquals(generator.getLibraryId(record(header, "rgNotInHeader")), unknown);
        Assert.assertEquals(generator.getLibraryId(record(header, "rgC")), libC);
        Assert.assertNotEquals(libA, libC);
        Assert.assertNotEquals(libA, unknown);
        Assert.assertNotEquals(libC, unknown);

        final Map<String, Short> ids = generator.getLibraryIdsMap();
        Assert.assertEquals(ids.size(), 3);
        Assert.assertEquals(ids.get("libC"), Short.valueOf(libC));
        Assert.assertEquals(ids.get(UNKNOWN), Short.valueOf(unknown));
        Assert.assertEquals(ids.get("libA"), Short.valueOf(libA));
    }

    @Test
    public void testLibraryIdsWithoutReadGroups() {
        final SAMFileHeader header = header();
        final LibraryIdGenerator generator = new LibraryIdGenerator(header);
        Assert.assertTrue(generator.getMetricsByLibraryMap().isEmpty());
        final short unknown = generator.getLibraryId(record(header, null));
        Assert.assertEquals(generator.getLibraryId(record(header, "rgX")), unknown);
        Assert.assertEquals(generator.getLibraryIdsMap().size(), 1);
        Assert.assertNull(generator.getMetricsByLibrary(UNKNOWN));
    }

    @DataProvider(name = "flowMetrics")
    public Object[][] flowMetrics() {
        return new Object[][]{{false}, {true}};
    }

    /** One metrics object per distinct library in the header, keyed and iterated in library-name order. */
    @Test(dataProvider = "flowMetrics")
    public void testMetricsByLibraryFromHeader(final boolean flowMetrics) {
        final SAMFileHeader header = header(new String[]{"rg1", "zeta"}, new String[]{"rg2", "alpha"},
                new String[]{"rg3", "zeta"}, new String[]{"rg4", null}, new String[]{"rg5", "Mu"});
        final LibraryIdGenerator generator = new LibraryIdGenerator(header, flowMetrics);

        final Map<String, DuplicationMetrics> metrics = generator.getMetricsByLibraryMap();
        Assert.assertEquals(new ArrayList<>(metrics.keySet()), Arrays.asList("Mu", UNKNOWN, "alpha", "zeta"));
        for (final Map.Entry<String, DuplicationMetrics> entry : metrics.entrySet()) {
            final DuplicationMetrics m = entry.getValue();
            Assert.assertEquals(m.LIBRARY, entry.getKey());
            Assert.assertEquals(m instanceof FlowBasedDuplicationMetrics, flowMetrics);
            Assert.assertEquals(m.READ_PAIRS_EXAMINED, 0);
            Assert.assertEquals(m.UNPAIRED_READ_DUPLICATES, 0);
            Assert.assertSame(generator.getMetricsByLibrary(entry.getKey()), m);
        }
        Assert.assertNull(generator.getMetricsByLibrary("absent"));
        // creating the metrics does not assign library ids
        Assert.assertTrue(generator.getLibraryIdsMap().isEmpty());
    }

    @Test
    public void testDefaultConstructorMakesPlainMetrics() {
        final LibraryIdGenerator generator = new LibraryIdGenerator(header(new String[]{"rg1", "lib1"}));
        Assert.assertEquals(generator.getMetricsByLibrary("lib1").getClass(), DuplicationMetrics.class);
    }

    @Test
    public void testAddMetricsByLibrary() {
        final LibraryIdGenerator generator = new LibraryIdGenerator(header(new String[]{"rg1", "lib1"}));
        final DuplicationMetrics added = new DuplicationMetrics();
        added.LIBRARY = "lib0";
        generator.addMetricsByLibrary("lib0", added);
        Assert.assertSame(generator.getMetricsByLibrary("lib0"), added);
        Assert.assertEquals(new ArrayList<>(generator.getMetricsByLibraryMap().keySet()), Arrays.asList("lib0", "lib1"));

        final DuplicationMetrics replacement = new DuplicationMetrics();
        generator.addMetricsByLibrary("lib1", replacement);
        Assert.assertSame(generator.getMetricsByLibrary("lib1"), replacement);
        Assert.assertEquals(generator.getMetricsByLibraryMap().size(), 2);
    }

    @Test
    public void testHistograms() {
        final LibraryIdGenerator generator = new LibraryIdGenerator(header(new String[]{"rg1", "lib1"}));

        assertEmptyHistogram(generator.getDuplicateCountHist(), "set_size", "all_sets");
        assertEmptyHistogram(generator.getNonOpticalDuplicateCountHist(), "set_size", "non_optical_sets");
        assertEmptyHistogram(generator.getOpticalDuplicateCountHist(), "set_size", "optical_sets");
        Assert.assertTrue(generator.getOpticalDuplicatesByLibraryIdMap().isEmpty());
        Assert.assertEquals(generator.getNumberOfOpticalDuplicateClusters(), 0L);

        // the getters return the live histograms
        generator.getDuplicateCountHist().increment(2.0);
        Assert.assertEquals(generator.getDuplicateCountHist().get(2.0).getValue(), 1.0);
        Assert.assertTrue(generator.getNonOpticalDuplicateCountHist().isEmpty());
        Assert.assertTrue(generator.getOpticalDuplicateCountHist().isEmpty());
    }

    private static void assertEmptyHistogram(final Histogram<Double> histogram, final String binLabel, final String valueLabel) {
        Assert.assertTrue(histogram.isEmpty());
        Assert.assertEquals(histogram.getBinLabel(), binLabel);
        Assert.assertEquals(histogram.getValueLabel(), valueLabel);
    }

    /** The number of optical duplicate clusters is the sum over all libraries. */
    @Test
    public void testNumberOfOpticalDuplicateClusters() {
        final LibraryIdGenerator generator = new LibraryIdGenerator(header(new String[]{"rg1", "lib1"}, new String[]{"rg2", "lib2"}));
        final Histogram<Short> byLibrary = generator.getOpticalDuplicatesByLibraryIdMap();
        byLibrary.increment((short) 1, 3.0);
        byLibrary.increment((short) 2);
        byLibrary.increment((short) 2);
        Assert.assertEquals(generator.getNumberOfOpticalDuplicateClusters(), 5L);
        Assert.assertSame(generator.getOpticalDuplicatesByLibraryIdMap(), byLibrary);
    }
}
