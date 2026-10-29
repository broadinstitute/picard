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
import picard.PicardException;

import java.util.List;

import static picard.sam.markduplicates.MateCigarAndUmiTestSupport.*;

/**
 * Tests of {@link UmiUtil}: top-strand normalization of single and duplex UMIs and the UMIs it rejects, the strand of
 * a read from the unclipped 5' ends of it and its mate, the molecular identifier it writes, and UMI lengths.
 */
public class UmiUtilBehaviorTest {

    private static SAMRecordSetBuilder builder() {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        return builder;
    }

    /** The first end of a pair on the top strand: forward at 1000, mate reverse at 1200. */
    private static SAMRecord topRead(final String umi) {
        final SAMRecord record = addPair(builder(), "top", 0, 1000, false, "50M", 1200, true, "50M").get(0);
        record.setAttribute("RX", umi);
        return record;
    }

    /** The first end of a pair on the bottom strand: reverse at 1200, mate forward at 1000. */
    private static SAMRecord bottomRead(final String umi) {
        final SAMRecord record = addPair(builder(), "bottom", 0, 1200, true, "50M", 1000, false, "50M").get(0);
        record.setAttribute("RX", umi);
        return record;
    }

    /** A mapped first end whose mate is unmapped, so its strand is unknown. */
    private static SAMRecord unknownStrandRead(final String umi) {
        final SAMRecord record = addPairWithUnmappedMate(builder(), "unknown", 0, 1000).get(0);
        record.setAttribute("RX", umi);
        return record;
    }

    @DataProvider(name = "normalizedUmis")
    public Object[][] normalizedUmis() {
        return new Object[][]{
                {topRead("AC-GT"), false, "AC-GT"},
                {bottomRead("AC-GT"), false, "AC-GT"},
                {bottomRead("ACGT"), false, "ACGT"},
                {topRead("AC-GT"), true, "AC-GT"},
                {bottomRead("AC-GT"), true, "GT-AC"},
                {bottomRead("ACN-GTN"), true, "GTN-ACN"},
                {unknownStrandRead("AC-GT"), true, "AC-GT"},
                {topRead(""), false, ""},
                {topRead(null), false, null},
                {topRead(null), true, null},
        };
    }

    @Test(dataProvider = "normalizedUmis")
    public void testTopStrandNormalizedUmi(final SAMRecord record, final boolean duplex, final String expected) {
        Assert.assertEquals(UmiUtil.getTopStrandNormalizedUmi(record, "RX", duplex), expected);
    }

    @Test
    public void testNullTagGivesNull() {
        Assert.assertNull(UmiUtil.getTopStrandNormalizedUmi(topRead("AC-GT"), null, true));
        Assert.assertNull(UmiUtil.getTopStrandNormalizedUmi(topRead("AC-GT"), null, false));
    }

    @DataProvider(name = "rejectedUmis")
    public Object[][] rejectedUmis() {
        final String illegal = "illegal characters";
        final String notTwoParts = "X-Y";
        return new Object[][]{
                {"ACGU", false, illegal},
                {"AC_GT", false, illegal},
                {"AC GT", false, illegal},
                {"AC*T", true, illegal},
                {"AC+GT", true, illegal},
                {"ACGT", true, notTwoParts},
                {"A-C-G", true, notTwoParts},
                {"ACG-", true, notTwoParts},
        };
    }

    @Test(dataProvider = "rejectedUmis")
    public void testRejectedUmis(final String umi, final boolean duplex, final String messageFragment) {
        final PicardException exception = Assert.expectThrows(PicardException.class,
                () -> UmiUtil.getTopStrandNormalizedUmi(topRead(umi), "RX", duplex));
        Assert.assertTrue(exception.getMessage().contains(messageFragment), exception.getMessage());
    }

    @DataProvider(name = "strands")
    public Object[][] strands() {
        return new Object[][]{
                // contig, start, negative, CIGAR of each end; which end; expected strand
                {0, 1000, false, "50M", 0, 1200, true, "50M", 0, UmiUtil.ReadStrand.TOP},
                {0, 1000, false, "50M", 0, 1200, true, "50M", 1, UmiUtil.ReadStrand.TOP},
                {0, 1200, true, "50M", 0, 1000, false, "50M", 0, UmiUtil.ReadStrand.BOTTOM},
                {0, 1200, true, "50M", 0, 1000, false, "50M", 1, UmiUtil.ReadStrand.BOTTOM},
                // The read's own soft clip moves its 5' end from 105 to 95, before the mate's 5' end at 99.
                {0, 105, false, "10S40M", 0, 50, true, "50M", 0, UmiUtil.ReadStrand.TOP},
                {0, 105, false, "10S40M", 0, 50, true, "50M", 1, UmiUtil.ReadStrand.TOP},
                // The mate's soft clip moves its 5' end from 96 to 101, after the read's 5' end at 100.
                {0, 100, false, "50M", 0, 57, true, "40M5S", 0, UmiUtil.ReadStrand.TOP},
                {0, 100, false, "50M", 0, 57, true, "40M5S", 1, UmiUtil.ReadStrand.TOP},
                {0, 100, false, "50M", 0, 52, true, "45M", 0, UmiUtil.ReadStrand.BOTTOM},
                // On different contigs the contig order decides.
                {0, 5000, true, "50M", 1, 100, false, "50M", 0, UmiUtil.ReadStrand.TOP},
                {0, 5000, true, "50M", 1, 100, false, "50M", 1, UmiUtil.ReadStrand.TOP},
                {2, 100, false, "50M", 1, 5000, true, "50M", 0, UmiUtil.ReadStrand.BOTTOM},
                {2, 100, false, "50M", 1, 5000, true, "50M", 1, UmiUtil.ReadStrand.BOTTOM},
        };
    }

    @Test(dataProvider = "strands")
    public void testGetStrand(final int contig1, final int start1, final boolean negative1, final String cigar1,
                              final int contig2, final int start2, final boolean negative2, final String cigar2,
                              final int end, final UmiUtil.ReadStrand expected) {
        final List<SAMRecord> pair = builder().addPair("pair", contig1, contig2, start1, start2, false, false,
                cigar1, cigar2, negative1, negative2, false, false, QUALITY);
        Assert.assertEquals(UmiUtil.getStrand(pair.get(end)), expected);
    }

    @Test
    public void testStrandUnknownWithAnUnmappedEnd() {
        final List<SAMRecord> pair = addPairWithUnmappedMate(builder(), "pair", 0, 1000);
        Assert.assertEquals(UmiUtil.getStrand(pair.get(0)), UmiUtil.ReadStrand.UNKNOWN);
        Assert.assertEquals(UmiUtil.getStrand(pair.get(1)), UmiUtil.ReadStrand.UNKNOWN);
    }

    @Test
    public void testSetMolecularIdentifier() {
        final SAMRecordSetBuilder builder = builder();
        final List<SAMRecord> fr = addPair(builder, "fr", 2, 1000, 1200);
        final List<SAMRecord> rf = addPair(builder, "rf", 0, 1200, true, "50M", 1000, false, "50M");
        final SAMRecord mateUnmapped = addPairWithUnmappedMate(builder, "mu", 0, 3000).get(0);

        // Forward reads use the mate's start and reverse reads their own, so both ends of a pair agree.
        UmiUtil.setMolecularIdentifier(fr.get(0), "ACGT", "MI", false);
        UmiUtil.setMolecularIdentifier(fr.get(1), "ACGT", "MI", false);
        Assert.assertEquals(fr.get(0).getAttribute("MI"), "chr3:1200/ACGT");
        Assert.assertEquals(fr.get(1).getAttribute("MI"), "chr3:1200/ACGT");

        UmiUtil.setMolecularIdentifier(fr.get(0), "AC-GT", "MI", true);
        UmiUtil.setMolecularIdentifier(fr.get(1), "AC-GT", "MI", true);
        UmiUtil.setMolecularIdentifier(rf.get(0), "AC-GT", "MI", true);
        UmiUtil.setMolecularIdentifier(rf.get(1), "AC-GT", "MI", true);
        UmiUtil.setMolecularIdentifier(mateUnmapped, "AC-GT", "MI", true);
        Assert.assertEquals(fr.get(0).getAttribute("MI"), "chr3:1200/AC-GT/A");
        Assert.assertEquals(fr.get(1).getAttribute("MI"), "chr3:1200/AC-GT/A");
        Assert.assertEquals(rf.get(0).getAttribute("MI"), "chr1:1200/AC-GT/B");
        Assert.assertEquals(rf.get(1).getAttribute("MI"), "chr1:1200/AC-GT/B");
        Assert.assertEquals(mateUnmapped.getAttribute("MI"), "chr1:3000/AC-GT");

        UmiUtil.setMolecularIdentifier(rf.get(0), "TTTT", "XM", false);
        Assert.assertEquals(rf.get(0).getAttribute("XM"), "chr1:1200/TTTT");

        final SAMRecord untagged = addPair(builder, "untagged", 0, 1000, 1200).get(0);
        final int attributes = untagged.getAttributes().size();
        UmiUtil.setMolecularIdentifier(untagged, "ACGT", null, true);
        Assert.assertEquals(untagged.getAttributes().size(), attributes);
        Assert.assertNull(untagged.getAttribute("MI"));
    }

    @DataProvider(name = "umiLengths")
    public Object[][] umiLengths() {
        return new Object[][]{
                {"", 0},
                {"ACGT", 4},
                {"ACGTN", 5},
                {"AC-GT", 4},
                {"AAA-CC", 5},
                {"A-C-G-T", 4},
                {"--", 0},
        };
    }

    @Test(dataProvider = "umiLengths")
    public void testGetUmiLength(final String umi, final int expected) {
        Assert.assertEquals(UmiUtil.getUmiLength(umi), expected);
    }
}
