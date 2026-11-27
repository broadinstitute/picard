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
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.markduplicates.MarkDuplicatesForFlowArgumentCollection.FLOW_DUPLICATE_SELECTION_STRATEGY;
import picard.sam.markduplicates.util.LibraryIdGenerator;
import picard.sam.markduplicates.util.ReadEnds;
import picard.sam.markduplicates.util.ReadEndsForMarkDuplicates;
import picard.sam.markduplicates.util.ReadEndsForMarkDuplicatesWithBarcodes;

import java.util.Arrays;

/**
 * Tests of {@link MarkDuplicatesForFlowHelper} and {@link MarkDuplicatesForFlowArgumentCollection} below the
 * command line: the start and end coordinates of a read under each flow argument and tm clipping tag, the read ends
 * and scores built for a read, the score of a pair end, the comparability test, the clipping tag predicates, both
 * quality strategies, the argument defaults, and the argument combination the helper rejects.
 */
public class MarkDuplicatesForFlowHelperBehaviorTest {

    private static final String FLOW_ORDER = "TGCA";

    private static String acgt(final int length) {
        final StringBuilder sb = new StringBuilder(length);
        for (int i = 0; i < length; i++) {
            sb.append("ACGT".charAt(i % 4));
        }
        return sb.toString();
    }

    private static byte[] quals(final int length, final int q) {
        final byte[] out = new byte[length];
        Arrays.fill(out, (byte) q);
        return out;
    }

    private static byte[] qualities(final int... values) {
        final byte[] out = new byte[values.length];
        for (int i = 0; i < values.length; i++) {
            out[i] = (byte) values[i];
        }
        return out;
    }

    /** "own": the record's read group has the flow order; "other": only another read group has it; "none": no flow order. */
    private static SAMFileHeader header(final String flowOrderSource) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        header.addSequence(new SAMSequenceRecord("chr1", 100_000));
        final SAMReadGroupRecord rgA = new SAMReadGroupRecord("rgA");
        rgA.setSample("sample");
        rgA.setLibrary("libA");
        if (flowOrderSource.equals("own") || flowOrderSource.equals("noReadGroupOnRecord")) {
            rgA.setFlowOrder(FLOW_ORDER);
        }
        header.addReadGroup(rgA);
        if (flowOrderSource.equals("other")) {
            final SAMReadGroupRecord rgB = new SAMReadGroupRecord("rgB");
            rgB.setSample("sample");
            rgB.setLibrary("libB");
            rgB.setFlowOrder(FLOW_ORDER);
            header.addReadGroup(rgB);
        }
        return header;
    }

    private static SAMRecord record(final SAMFileHeader header, final int start, final String cigar, final boolean negative,
                                    final String bases, final byte[] quals, final String tm) {
        final SAMRecord rec = new SAMRecord(header);
        rec.setReadName("read");
        rec.setReferenceIndex(0);
        rec.setAlignmentStart(start);
        rec.setCigarString(cigar);
        rec.setReadNegativeStrandFlag(negative);
        rec.setMappingQuality(60);
        rec.setReadBases(bases.getBytes());
        rec.setBaseQualities(quals);
        rec.setAttribute("RG", "rgA");
        if (tm != null) {
            rec.setAttribute(MarkDuplicatesForFlowHelper.CLIPPING_TAG_NAME, tm);
        }
        return rec;
    }

    /** Comma-separated: CLIPPED, QKNOWN, USEEND, SKIP=n, ENDQ (end quality strategy), THRESHOLD=n. */
    private static MarkDuplicatesForFlowArgumentCollection arguments(final String spec) {
        final MarkDuplicatesForFlowArgumentCollection args = new MarkDuplicatesForFlowArgumentCollection();
        args.FLOW_MODE = true;
        for (final String token : spec.split(",")) {
            if (token.isEmpty()) {
                continue;
            }
            if (token.equals("CLIPPED")) {
                args.FLOW_USE_UNPAIRED_CLIPPED_END = true;
            } else if (token.equals("QKNOWN")) {
                args.FLOW_Q_IS_KNOWN_END = true;
            } else if (token.equals("USEEND")) {
                args.FLOW_USE_END_IN_UNPAIRED_READS = true;
            } else if (token.equals("ENDQ")) {
                args.FLOW_DUP_STRATEGY = FLOW_DUPLICATE_SELECTION_STRATEGY.FLOW_END_QUALITY_STRATEGY;
            } else if (token.startsWith("SKIP=")) {
                args.FLOW_SKIP_FIRST_N_FLOWS = Integer.parseInt(token.substring(5));
            } else if (token.startsWith("THRESHOLD=")) {
                args.FLOW_EFFECTIVE_QUALITY_THRESHOLD = Integer.parseInt(token.substring(10));
            } else {
                throw new IllegalArgumentException(token);
            }
        }
        return args;
    }

    private static final String TTGCA = ("TTGCA" + acgt(40).substring(1)).substring(0, 40);
    private static final String TCAGT = ("TCAGT" + acgt(40).substring(1)).substring(0, 40);

    @DataProvider(name = "coordinates")
    public Object[][] coordinates() {
        // "2S40M3S" at 100 aligns 100-139; its unclipped span is 98-142
        final String clipped = "2S40M3S";
        final String b45 = acgt(45);
        return new Object[][]{
                {"own", clipped, false, b45, null, true, true, "", 98},
                {"own", clipped, false, b45, null, false, false, "", 142},
                {"own", clipped, false, b45, null, false, false, "CLIPPED", 139},
                {"own", clipped, false, b45, "Z", false, false, "CLIPPED", 139},
                // a trimmed end that is not adapter clipped is unknown
                {"own", clipped, false, b45, "Z", false, false, "", 0},
                {"own", clipped, false, b45, "Z", true, true, "", 98},
                {"own", clipped, false, b45, "Z", false, true, "", 142},
                {"own", clipped, false, b45, "A", false, false, "", 142},
                {"own", clipped, false, b45, "A", false, false, "QKNOWN", 142},
                {"own", clipped, false, b45, "AZ", false, false, "", 142},
                {"own", clipped, false, b45, "AZ", false, false, "QKNOWN", 142},
                {"own", clipped, false, b45, "X", false, false, "", 142},
                // reverse strand: the 5' end is the unclipped end, the 3' end the unclipped start
                {"own", clipped, true, b45, null, false, true, "", 142},
                {"own", clipped, true, b45, null, true, false, "", 98},
                {"own", clipped, true, b45, "Z", true, false, "", 0},
                {"own", clipped, true, b45, null, true, false, "CLIPPED", 100},
                // without a flow order every flow argument and the tm tag are ignored
                {"none", clipped, false, b45, "Z", false, false, "", 142},
                {"none", clipped, false, b45, null, true, true, "CLIPPED", 98},
                {"none", clipped, false, b45, null, false, false, "CLIPPED", 142},
                {"none", "40M", false, TTGCA, null, true, true, "SKIP=2", 100},
                {"other", clipped, false, b45, "Z", false, false, "", 0},
                {"other", "40M", false, TTGCA, null, true, true, "SKIP=2", 103},
                {"noReadGroupOnRecord", clipped, false, b45, "Z", false, false, "", 0},
                // skipping flows of TGCA: TT is the first flow, G the second, C the third
                {"own", "40M", false, TTGCA, null, true, true, "SKIP=1", 102},
                {"own", "40M", false, TTGCA, null, true, true, "SKIP=2", 103},
                {"own", "40M", false, TTGCA, null, true, true, "SKIP=3", 104},
                // TCAG: the G flow between T and C is empty
                {"own", "40M", false, TCAGT, null, true, true, "SKIP=1", 101},
                {"own", "40M", false, TCAGT, null, true, true, "SKIP=2", 101},
                {"own", "40M", false, TCAGT, null, true, true, "SKIP=3", 102},
                {"own", "2S38M", false, TTGCA, null, true, true, "SKIP=1", 100},
                // skipping applies to the 5' end only
                {"own", "40M", false, TTGCA, null, false, false, "SKIP=2", 139},
                {"own", "40M", false, TTGCA, "Z", false, false, "SKIP=2", 0},
        };
    }

    /**
     * A record whose read group has no flow order, or that has no read group, uses the flow order of another read group
     * (the "other" and "noReadGroupOnRecord" rows). The documentation does not specify this; the expected value follows
     * the current implementation.
     */
    @Test(dataProvider = "coordinates")
    public void testGetReadEndCoordinate(final String flowOrderSource, final String cigar, final boolean negative,
                                         final String bases, final String tm, final boolean startEnd, final boolean certain,
                                         final String args, final int expected) {
        final SAMRecord rec = record(header(flowOrderSource), 100, cigar, negative, bases, quals(bases.length(), 30), tm);
        if (flowOrderSource.equals("noReadGroupOnRecord")) {
            rec.setAttribute("RG", null);
        }
        Assert.assertEquals(MarkDuplicatesForFlowHelper.getReadEndCoordinate(rec, startEnd, certain, arguments(args)), expected);
    }

    private static MarkDuplicates markDuplicates(final SAMFileHeader header, final String spec) {
        final MarkDuplicates md = new MarkDuplicates();
        md.flowBasedArguments = arguments(spec);
        md.setupOpticalDuplicateFinder();
        md.libraryIdGenerator = new LibraryIdGenerator(header, true);
        return md;
    }

    private static ReadEndsForMarkDuplicates buildReadEnds(final SAMRecord rec, final String spec) {
        final SAMFileHeader header = rec.getHeader();
        return new MarkDuplicatesForFlowHelper(markDuplicates(header, spec)).buildReadEnds(header, 7, rec, false);
    }

    @DataProvider(name = "readEnds")
    public Object[][] readEnds() {
        final String b40 = acgt(40);
        final String b45 = acgt(45);
        return new Object[][]{
                // {cigar, negative, bases, tm, arguments, read1Coordinate, read2Coordinate, orientation, score}
                {"40M", false, b40, null, "", 100, -1, ReadEnds.F, 1200},
                {"40M", false, b40, null, "USEEND", 100, 139, ReadEnds.F, 1200},
                {"40M", true, b40, null, "USEEND", 139, 100, ReadEnds.R, 1200},
                {"2S40M3S", true, b45, null, "USEEND", 142, 98, ReadEnds.R, 1350},
                {"40M5S", false, b45, null, "USEEND,CLIPPED", 100, 139, ReadEnds.F, 1350},
                // an unknown end scores -1 so that a read with a known end is kept over it
                {"40M", false, b40, "Z", "USEEND", 100, 0, ReadEnds.F, -1},
                {"40M", false, b40, "Z", "", 100, -1, ReadEnds.F, 1200},
                {"40M", false, b40, "A", "USEEND", 100, 139, ReadEnds.F, 1200},
                {"40M", false, TTGCA, null, "SKIP=1", 102, -1, ReadEnds.F, 1200},
                {"40M", false, b40, null, "THRESHOLD=31", 100, -1, ReadEnds.F, 0},
                {"40M", false, b40, null, "THRESHOLD=29", 100, -1, ReadEnds.F, 1200},
                {"40M", false, b40, null, "ENDQ", 100, -1, ReadEnds.F, 30},
        };
    }

    @Test(dataProvider = "readEnds")
    public void testBuildReadEnds(final String cigar, final boolean negative, final String bases, final String tm,
                                  final String spec, final int read1Coordinate, final int read2Coordinate,
                                  final byte orientation, final int score) {
        final SAMRecord rec = record(header("own"), 100, cigar, negative, bases, quals(bases.length(), 30), tm);
        final ReadEndsForMarkDuplicates ends = buildReadEnds(rec, spec);
        Assert.assertEquals(ends.read1ReferenceIndex, 0);
        Assert.assertEquals(ends.read1Coordinate, read1Coordinate);
        Assert.assertEquals(ends.read2Coordinate, read2Coordinate);
        Assert.assertEquals(ends.orientation, orientation);
        Assert.assertEquals(ends.score, score);
        Assert.assertEquals(ends.read1IndexInFile, 7);
        Assert.assertEquals(ends.read2ReferenceIndex, -1);
        Assert.assertFalse(ends.isPaired());
        Assert.assertEquals(ends.getLibraryId(), 1);
    }

    @DataProvider(name = "endQualityScores")
    public Object[][] endQualityScores() {
        // 40 bases without homopolymers: the end strategy looks at bases 0-9 and 30-39
        return new Object[][]{{9, 5}, {10, 30}, {29, 30}, {30, 5}, {0, 5}, {39, 5}};
    }

    @Test(dataProvider = "endQualityScores")
    public void testEndQualityStrategyLooksTenBasesFromEachEnd(final int lowBase, final int score) {
        final byte[] q = quals(40, 30);
        q[lowBase] = 5;
        final SAMRecord rec = record(header("own"), 100, "40M", false, acgt(40), q, null);
        Assert.assertEquals(buildReadEnds(rec, "ENDQ").score, score);
    }

    @Test
    public void testScoreAdjustments() {
        final SAMFileHeader header = header("own");
        final SAMRecord failed = record(header, 100, "40M", false, acgt(40), quals(40, 30), null);
        failed.setReadFailsVendorQualityCheckFlag(true);
        Assert.assertEquals(buildReadEnds(failed, "").score, 1200 + Short.MIN_VALUE / 2);

        // capped at half the short range
        final SAMRecord longRead = record(header, 100, "600M", false, acgt(600), quals(600, 40), null);
        Assert.assertEquals(buildReadEnds(longRead, "").score, Short.MAX_VALUE / 2);
    }

    @Test
    public void testPairedReadsWithMappedMatesAreRejected() {
        final SAMRecord rec = record(header("own"), 100, "40M", false, acgt(40), quals(40, 30), null);
        rec.setReadPairedFlag(true);
        rec.setFirstOfPairFlag(true);
        rec.setMateReferenceIndex(0);
        rec.setMateAlignmentStart(300);
        try {
            buildReadEnds(rec, "");
            Assert.fail("expected an IllegalArgumentException");
        } catch (final IllegalArgumentException e) {
            Assert.assertTrue(e.getMessage().contains("paired reads"), e.getMessage());
        }

        rec.setMateUnmappedFlag(true);
        final ReadEndsForMarkDuplicates ends = buildReadEnds(rec, "USEEND");
        Assert.assertFalse(ends.isPaired());
        Assert.assertEquals(ends.read1Coordinate, 100);
        Assert.assertEquals(ends.read2Coordinate, 139);
    }

    /** A homopolymer read whose flow score (1060) differs from its base-quality sum (700). */
    private static SAMRecord homopolymerRead() {
        final byte[] q = quals(40, 22);
        q[0] = 40;
        Arrays.fill(q, 1, 10, (byte) 2);
        return record(header("own"), 100, "40M", false, "AAAAAAAAAA" + acgt(31).substring(1), q, null);
    }

    @Test
    public void testReadDuplicateScore() {
        final SAMFileHeader header = header("own");
        final MarkDuplicatesForFlowHelper flow = new MarkDuplicatesForFlowHelper(markDuplicates(header, ""));
        final ReadEndsForMarkDuplicates ends = new ReadEndsForMarkDuplicates();

        ends.read2Coordinate = 139;
        Assert.assertEquals(flow.getReadDuplicateScore(homopolymerRead(), ends), 1060);
        ends.read2Coordinate = -1;
        Assert.assertEquals(flow.getReadDuplicateScore(homopolymerRead(), ends), 1060);
        ends.read2Coordinate = 0;
        Assert.assertEquals(flow.getReadDuplicateScore(homopolymerRead(), ends), -1);

        final MarkDuplicates notFlow = markDuplicates(header, "");
        notFlow.flowBasedArguments.FLOW_MODE = false;
        ends.read2Coordinate = 139;
        Assert.assertEquals(new MarkDuplicatesForFlowHelper(notFlow).getReadDuplicateScore(homopolymerRead(), ends), 700);
    }

    private static ReadEndsForMarkDuplicatesWithBarcodes ends(final short library, final int reference, final byte orientation,
                                                              final int read2Reference, final int read2Coordinate,
                                                              final int barcode, final int readOneBarcode, final int readTwoBarcode) {
        final ReadEndsForMarkDuplicatesWithBarcodes ends = new ReadEndsForMarkDuplicatesWithBarcodes();
        ends.libraryId = library;
        ends.read1ReferenceIndex = reference;
        ends.read1Coordinate = 100;
        ends.orientation = orientation;
        ends.read2ReferenceIndex = read2Reference;
        ends.read2Coordinate = read2Coordinate;
        ends.barcode = barcode;
        ends.readOneBarcode = readOneBarcode;
        ends.readTwoBarcode = readTwoBarcode;
        return ends;
    }

    @DataProvider(name = "comparable")
    public Object[][] comparable() {
        final ReadEndsForMarkDuplicatesWithBarcodes base = ends((short) 1, 0, ReadEnds.F, -1, 500, 7, 8, 9);
        final ReadEndsForMarkDuplicatesWithBarcodes otherStart = ends((short) 1, 0, ReadEnds.F, -1, 500, 7, 8, 9);
        otherStart.read1Coordinate = 250;
        return new Object[][]{
                // {rhs, compareRead2, useBarcodes, comparable}
                {ends((short) 1, 0, ReadEnds.F, -1, 500, 7, 8, 9), true, true, true},
                {ends((short) 2, 0, ReadEnds.F, -1, 500, 7, 8, 9), false, false, false},
                {ends((short) 1, 1, ReadEnds.F, -1, 500, 7, 8, 9), false, false, false},
                {ends((short) 1, 0, ReadEnds.R, -1, 500, 7, 8, 9), false, false, false},
                // the start is not compared here
                {otherStart, true, true, true},
                {ends((short) 1, 0, ReadEnds.F, -1, 501, 7, 8, 9), false, false, true},
                {ends((short) 1, 0, ReadEnds.F, -1, 501, 7, 8, 9), true, false, false},
                {ends((short) 1, 0, ReadEnds.F, 0, 500, 7, 8, 9), false, false, true},
                {ends((short) 1, 0, ReadEnds.F, 0, 500, 7, 8, 9), true, false, false},
                {ends((short) 1, 0, ReadEnds.F, -1, 500, 6, 8, 9), false, false, true},
                {ends((short) 1, 0, ReadEnds.F, -1, 500, 6, 8, 9), false, true, false},
                {ends((short) 1, 0, ReadEnds.F, -1, 500, 7, 0, 9), false, true, false},
                {ends((short) 1, 0, ReadEnds.F, -1, 500, 7, 8, 0), false, true, false},
                {ends((short) 1, 0, ReadEnds.F, -1, 500, 7, 8, 9), false, true, true},
        };
    }

    @Test(dataProvider = "comparable")
    public void testAreComparableForDuplicates(final ReadEndsForMarkDuplicatesWithBarcodes rhs, final boolean compareRead2,
                                               final boolean useBarcodes, final boolean comparable) {
        final MarkDuplicatesForFlowHelper helper = new MarkDuplicatesForFlowHelper(markDuplicates(header("own"), ""));
        final ReadEndsForMarkDuplicatesWithBarcodes lhs = ends((short) 1, 0, ReadEnds.F, -1, 500, 7, 8, 9);
        Assert.assertEquals(helper.areComparableForDuplicates(lhs, rhs, compareRead2, useBarcodes), comparable);
        Assert.assertEquals(helper.areComparableForDuplicates(rhs, lhs, compareRead2, useBarcodes), comparable);
    }

    @DataProvider(name = "clippingTags")
    public Object[][] clippingTags() {
        return new Object[][]{
                // {tm, adapter clipped, adapter clipped with Q, quality clipped}
                {null, false, false, false},
                {"A", true, true, false},
                {"Q", false, true, true},
                {"Z", false, false, true},
                {"AQ", true, true, true},
                {"QZ", false, true, true},
                {"AZ", true, true, true},
                {"X", false, false, false},
                {"aqz", false, false, false},
        };
    }

    @Test(dataProvider = "clippingTags")
    public void testClippingTagPredicates(final String tm, final boolean adapter, final boolean adapterWithQ, final boolean quality) {
        final SAMRecord rec = record(header("own"), 100, "40M", false, acgt(40), quals(40, 30), tm);
        Assert.assertEquals(MarkDuplicatesForFlowHelper.isAdapterClipped(rec), adapter);
        Assert.assertEquals(MarkDuplicatesForFlowHelper.isAdapterClippedWithQ(rec), adapterWithQ);
        Assert.assertEquals(MarkDuplicatesForFlowHelper.isQualityClipped(rec), quality);
    }

    @DataProvider(name = "qualitySum")
    public Object[][] qualitySum() {
        return new Object[][]{
                // {bases, quals, negative strand, threshold, score}
                {"AAC", qualities(10, 20, 30), false, 0, 50},
                // read in synthesis order the homopolymer starts at its right end
                {"AAC", qualities(10, 20, 30), true, 0, 70},
                {"AC", qualities(16, 14), false, 15, 16},
                {"ACCC", qualities(20, 40, 10, 10), false, 30, 120},
                {"ACCC", qualities(20, 40, 10, 10), true, 30, 0},
                {"ACCC", qualities(20, 40, 10, 10), true, 9, 50},
                {"A", qualities(7), false, 0, 7},
        };
    }

    @Test(dataProvider = "qualitySum")
    public void testFlowSumOfBaseQualities(final String bases, final byte[] quals, final boolean negative,
                                           final int threshold, final int score) {
        final SAMRecord rec = record(header("own"), 100, bases.length() + "M", negative, bases, quals, null);
        Assert.assertEquals(MarkDuplicatesForFlowHelper.getFlowSumOfBaseQualities(rec, threshold), score);
    }

    @DataProvider(name = "nearEnds")
    public Object[][] nearEnds() {
        return new Object[][]{
                // {bases, quals, distance, score}
                {"ACGTAC", qualities(30, 5, 30, 30, 30, 30), 1, 30},
                // the homopolymer at the start extends the window to base 1
                {"AACGTT", qualities(30, 5, 30, 30, 30, 30), 1, 5},
                // the homopolymer at the end extends the window to base 4
                {"ACGTAA", qualities(30, 30, 30, 30, 5, 30), 1, 5},
                {"ACGTAC", qualities(30, 30, 30, 30, 5, 30), 1, 30},
                {"ACG", qualities(20, 10, 30), 10, 10},
                {"AAAA", qualities(40, 30, 20, 10), 1, 10},
        };
    }

    @Test(dataProvider = "nearEnds")
    public void testFlowSumOfBaseQualitiesNearEnds(final String bases, final byte[] quals, final int distance, final int score) {
        final SAMRecord rec = record(header("own"), 100, bases.length() + "M", false, bases, quals, null);
        Assert.assertEquals(MarkDuplicatesForFlowHelper.getFlowSumOfBaseQualitiesNearEnds(rec, distance), score);
    }

    @DataProvider(name = "uncertaintyValidation")
    public Object[][] uncertaintyValidation() {
        return new Object[][]{
                {false, 0, false},
                {false, 1, true},
                {false, -1, true},
                {true, 0, false},
                {true, 5, false},
        };
    }

    @Test(dataProvider = "uncertaintyValidation")
    public void testEndUncertaintyRequiresUseEnd(final boolean useEnd, final int uncertainty, final boolean rejected) {
        final MarkDuplicates md = markDuplicates(header("own"), useEnd ? "USEEND" : "");
        md.flowBasedArguments.FLOW_UNPAIRED_END_UNCERTAINTY = uncertainty;
        try {
            new MarkDuplicatesForFlowHelper(md);
            Assert.assertFalse(rejected, "expected an IllegalArgumentException");
        } catch (final IllegalArgumentException e) {
            Assert.assertTrue(rejected, e.getMessage());
            Assert.assertTrue(e.getMessage().contains("UNPAIRED_END_UNCERTAINTY"), e.getMessage());
        }
    }

    @Test
    public void testArgumentDefaults() {
        final MarkDuplicatesForFlowArgumentCollection args = new MarkDuplicatesForFlowArgumentCollection();
        Assert.assertFalse(args.FLOW_MODE);
        Assert.assertEquals(args.FLOW_DUP_STRATEGY, FLOW_DUPLICATE_SELECTION_STRATEGY.FLOW_QUALITY_SUM_STRATEGY);
        Assert.assertFalse(args.FLOW_USE_END_IN_UNPAIRED_READS);
        Assert.assertFalse(args.FLOW_USE_UNPAIRED_CLIPPED_END);
        Assert.assertEquals(args.FLOW_UNPAIRED_END_UNCERTAINTY, 0);
        Assert.assertEquals(args.FLOW_UNPAIRED_START_UNCERTAINTY, 0);
        Assert.assertEquals(args.FLOW_SKIP_FIRST_N_FLOWS, 0);
        Assert.assertFalse(args.FLOW_Q_IS_KNOWN_END);
        Assert.assertEquals(args.FLOW_EFFECTIVE_QUALITY_THRESHOLD, 15);
        Assert.assertEquals(FLOW_DUPLICATE_SELECTION_STRATEGY.values(), new FLOW_DUPLICATE_SELECTION_STRATEGY[]{
                FLOW_DUPLICATE_SELECTION_STRATEGY.FLOW_QUALITY_SUM_STRATEGY, FLOW_DUPLICATE_SELECTION_STRATEGY.FLOW_END_QUALITY_STRATEGY});
        Assert.assertFalse(new MarkDuplicates().flowBasedArguments.FLOW_MODE);
    }
}
