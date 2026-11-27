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
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.sam.DuplicationMetrics;
import picard.sam.FlowBasedDuplicationMetrics;
import picard.sam.markduplicates.util.ReadEndsForMarkDuplicates;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/**
 * End-to-end tests of {@link MarkDuplicates} with FLOW_MODE=true on unpaired reads carrying flow-specific
 * information (a read group flow order, the tm clipping tag): which reads are flagged as duplicates and the
 * flow-based metrics under FLOW_USE_END_IN_UNPAIRED_READS, FLOW_USE_UNPAIRED_CLIPPED_END,
 * FLOW_UNPAIRED_END_UNCERTAINTY, FLOW_UNPAIRED_START_UNCERTAINTY, FLOW_SKIP_FIRST_N_FLOWS, FLOW_Q_IS_KNOWN_END,
 * FLOW_DUP_STRATEGY and FLOW_EFFECTIVE_QUALITY_THRESHOLD, with barcodes, libraries, spilled sorts, and the
 * argument and input errors flow mode raises.
 */
public class MarkDuplicatesFlowModeTest extends CommandLineProgramTest {

    private static final String FLOW_ORDER = "TGCA";
    private static final String RG_A = "rgA";

    @Override
    public String getCommandLineProgramName() {
        return MarkDuplicates.class.getSimpleName();
    }

    /** Bases with no homopolymer runs. */
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

    /** Builds a header from (read group, library, flow order) triples; a null flow order leaves FO unset. */
    private static SAMFileHeader header(final String... readGroupLibraryFlowOrder) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        header.addSequence(new SAMSequenceRecord("chr1", 100_000));
        header.addSequence(new SAMSequenceRecord("chr2", 100_000));
        for (int i = 0; i < readGroupLibraryFlowOrder.length; i += 3) {
            final SAMReadGroupRecord rg = new SAMReadGroupRecord(readGroupLibraryFlowOrder[i]);
            rg.setSample("sample");
            rg.setLibrary(readGroupLibraryFlowOrder[i + 1]);
            if (readGroupLibraryFlowOrder[i + 2] != null) {
                rg.setFlowOrder(readGroupLibraryFlowOrder[i + 2]);
            }
            header.addReadGroup(rg);
        }
        return header;
    }

    private static SAMFileHeader defaultHeader() {
        return header(RG_A, "libA", FLOW_ORDER);
    }

    private static SAMRecord read(final SAMFileHeader header, final String name, final int contig, final int start,
                                  final String cigar, final boolean negative, final String bases, final byte[] quals) {
        final SAMRecord rec = new SAMRecord(header);
        rec.setReadName(name);
        rec.setReferenceIndex(contig);
        rec.setAlignmentStart(start);
        rec.setCigarString(cigar);
        rec.setReadNegativeStrandFlag(negative);
        rec.setMappingQuality(60);
        rec.setReadBases(bases.getBytes());
        rec.setBaseQualities(quals);
        rec.setAttribute("RG", RG_A);
        return rec;
    }

    /** A forward read on chr1 of ACGT bases with flat qualities and a plain match CIGAR. */
    private static SAMRecord read(final SAMFileHeader header, final String name, final int start, final int length, final int q) {
        return read(header, name, 0, start, length + "M", false, acgt(length), quals(length, q));
    }

    private static final class FlowResult {
        final Set<String> duplicates = new TreeSet<>();
        final Set<String> names = new HashSet<>();
        final Map<String, SAMRecord> records = new HashMap<>();
        final Map<String, DuplicationMetrics> metricsByLibrary = new HashMap<>();

        FlowBasedDuplicationMetrics flow(final String library) {
            final DuplicationMetrics m = metricsByLibrary.get(library);
            Assert.assertTrue(m instanceof FlowBasedDuplicationMetrics, "metrics for " + library + ": " + m);
            return (FlowBasedDuplicationMetrics) m;
        }
    }

    private static File writeSorted(final SAMFileHeader header, final List<SAMRecord> records) throws IOException {
        final File input = File.createTempFile("mdFlow.", ".sam");
        input.deleteOnExit();
        final List<SAMRecord> sorted = new ArrayList<>(records);
        sorted.sort(Comparator.comparingInt(SAMRecord::getReferenceIndex).thenComparingInt(SAMRecord::getAlignmentStart));
        try (SAMFileWriter writer = new SAMFileWriterFactory().makeSAMWriter(header, true, input)) {
            sorted.forEach(writer::addAlignment);
        }
        return input;
    }

    private FlowResult run(final SAMFileHeader header, final List<SAMRecord> records, final String... extraArgs) throws IOException {
        final File input = writeSorted(header, records);
        final File output = File.createTempFile("mdFlow.", ".sam");
        output.deleteOnExit();
        final File metrics = File.createTempFile("mdFlow.", ".duplicate_metrics");
        metrics.deleteOnExit();
        final List<String> args = new ArrayList<>(Arrays.asList(
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + output.getAbsolutePath(),
                "METRICS_FILE=" + metrics.getAbsolutePath()));
        args.addAll(Arrays.asList(extraArgs));
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final FlowResult result = new FlowResult();
        try (SamReader reader = SamReaderFactory.makeDefault().open(output)) {
            for (final SAMRecord rec : reader) {
                result.names.add(rec.getReadName());
                result.records.put(rec.getReadName(), rec);
                if (rec.getDuplicateReadFlag()) {
                    result.duplicates.add(rec.getReadName());
                }
            }
        }
        for (final DuplicationMetrics m : MetricsFile.<DuplicationMetrics>readBeans(metrics)) {
            result.metricsByLibrary.put(m.LIBRARY, m);
        }
        return result;
    }

    private FlowResult runFlow(final SAMFileHeader header, final List<SAMRecord> records, final String... extraArgs) throws IOException {
        final List<String> args = new ArrayList<>();
        args.add("FLOW_MODE=true");
        args.addAll(Arrays.asList(extraArgs));
        return run(header, records, args.toArray(new String[0]));
    }

    private static Set<String> set(final String... names) {
        return new TreeSet<>(Arrays.asList(names));
    }

    private static void assertUnpairedMetrics(final FlowBasedDuplicationMetrics m, final long examined, final long duplicates,
                                              final long withTlen, final long dupsWithTlen, final long dupsWithoutTlen) {
        Assert.assertEquals(m.UNPAIRED_READS_EXAMINED, examined, "UNPAIRED_READS_EXAMINED");
        Assert.assertEquals(m.UNPAIRED_READ_DUPLICATES, duplicates, "UNPAIRED_READ_DUPLICATES");
        Assert.assertEquals(m.READ_PAIRS_EXAMINED, 0);
        Assert.assertEquals(m.READ_PAIR_DUPLICATES, 0);
        Assert.assertEquals(m.READ_PAIR_OPTICAL_DUPLICATES, 0);
        Assert.assertNull(m.ESTIMATED_LIBRARY_SIZE);
        Assert.assertEquals(m.PERCENT_DUPLICATION, examined == 0 ? 0.0 : duplicates / (double) examined, 1e-6);
        Assert.assertEquals(m.UNPAIRED_WITH_TLEN, withTlen, "UNPAIRED_WITH_TLEN");
        Assert.assertEquals(m.UNPAIRED_DUPS_WITH_TLEN, dupsWithTlen, "UNPAIRED_DUPS_WITH_TLEN");
        Assert.assertEquals(m.UNPAIRED_DUPS_WITHOUT_TLEN, dupsWithoutTlen, "UNPAIRED_DUPS_WITHOUT_TLEN");
        Assert.assertEquals(m.UNPAIRED_DUP_RATE_WITH_TLEN, withTlen == 0 ? 0.0 : dupsWithTlen / (double) withTlen, 1e-6);
        Assert.assertEquals(m.UNPAIRED_DUP_RATE_WITHOUT_TLEN,
                examined == withTlen ? 0.0 : dupsWithoutTlen / (double) (examined - withTlen), 1e-6);
    }

    @Test
    public void testStartAndOrientationDefineDuplicatesByDefault() throws IOException {
        final SAMFileHeader header = defaultHeader();
        final List<SAMRecord> records = Arrays.asList(
                read(header, "r1", 100, 40, 30),
                read(header, "r2", 100, 30, 30),
                // unclipped start 100
                read(header, "r3", 0, 105, "5S35M", false, acgt(40), quals(40, 29)),
                // reverse strand with its 5' end at 100
                read(header, "r4", 0, 61, "40M", true, acgt(40), quals(40, 30)),
                read(header, "r5", 101, 40, 30));
        final FlowResult result = runFlow(header, records);
        Assert.assertEquals(result.names, set("r1", "r2", "r3", "r4", "r5"));
        Assert.assertEquals(result.duplicates, set("r2", "r3"));
        Assert.assertEquals(result.metricsByLibrary.keySet(), Collections.singleton("libA"));
        assertUnpairedMetrics(result.flow("libA"), 5, 2, 0, 0, 2);
    }

    @Test
    public void testFlowScoreRatherThanBaseQualitySumPicksTheKeeper() throws IOException {
        // a 10-base homopolymer whose first base has quality 40 and the rest 2: the flow score counts all ten at 40
        final String bases = "AAAAAAAAAA" + acgt(32).substring(1, 31);
        final byte[] hmerQuals = quals(40, 22);
        hmerQuals[0] = 40;
        Arrays.fill(hmerQuals, 1, 10, (byte) 2);
        final SAMFileHeader header = defaultHeader();
        final List<SAMRecord> records = Arrays.asList(
                read(header, "x", 0, 200, "40M", false, bases, hmerQuals),
                read(header, "y", 0, 200, "40M", false, bases, quals(40, 25)));
        // flow sums: x 10 * 40 + 30 * 22 = 1060, y 40 * 25 = 1000; base-quality sums: x 40 + 660 = 700, y 1000
        Assert.assertEquals(runFlow(header, records).duplicates, set("y"));
        Assert.assertEquals(run(header, records).duplicates, set("x"));
    }

    @DataProvider(name = "qualityThreshold")
    public Object[][] qualityThreshold() {
        // x: 40 bases at 16 (sum 640 when counted); y: 20 bases at 30 and 20 at 3 (sum 600, or 660 when the 3s count)
        return new Object[][]{
                {null, set("y")},
                {"FLOW_EFFECTIVE_QUALITY_THRESHOLD=17", set("x")},
                {"FLOW_EFFECTIVE_QUALITY_THRESHOLD=1", set("x")},
        };
    }

    @Test(dataProvider = "qualityThreshold")
    public void testEffectiveQualityThreshold(final String arg, final Set<String> duplicates) throws IOException {
        final SAMFileHeader header = defaultHeader();
        final byte[] yQuals = quals(40, 3);
        Arrays.fill(yQuals, 0, 20, (byte) 30);
        final List<SAMRecord> records = Arrays.asList(
                read(header, "x", 300, 40, 16),
                read(header, "y", 0, 300, "40M", false, acgt(40), yQuals));
        final FlowResult result = arg == null ? runFlow(header, records) : runFlow(header, records, arg);
        Assert.assertEquals(result.duplicates, duplicates);
    }

    @DataProvider(name = "strategy")
    public Object[][] strategy() {
        return new Object[][]{
                {"FLOW_DUP_STRATEGY=FLOW_QUALITY_SUM_STRATEGY", set("x")},
                {"FLOW_DUP_STRATEGY=FLOW_END_QUALITY_STRATEGY", set("y")},
        };
    }

    @Test(dataProvider = "strategy")
    public void testDuplicateSelectionStrategy(final String arg, final Set<String> duplicates) throws IOException {
        // x: quality 2 in the middle (sum 900, lowest near the ends 30); y: one 10 near the start (sum 1170, ends 10)
        final byte[] xQuals = quals(40, 30);
        Arrays.fill(xQuals, 15, 25, (byte) 2);
        final byte[] yQuals = quals(40, 30);
        yQuals[3] = 10;
        final SAMFileHeader header = defaultHeader();
        final List<SAMRecord> records = Arrays.asList(
                read(header, "x", 0, 400, "40M", false, acgt(40), xQuals),
                read(header, "y", 0, 400, "40M", false, acgt(40), yQuals));
        Assert.assertEquals(runFlow(header, records, arg).duplicates, duplicates);
    }

    private static List<SAMRecord> endInput(final SAMFileHeader header) {
        return Arrays.asList(
                read(header, "e1", 500, 40, 30),
                read(header, "e2", 500, 40, 20),
                read(header, "e3", 500, 41, 30),
                // unclipped end equal to e1's
                read(header, "e4", 0, 500, "35M5S", false, acgt(40), quals(40, 25)),
                // reverse strand: 5' end 639 for both, 3' ends 600 and 599
                read(header, "v1", 0, 600, "40M", true, acgt(40), quals(40, 30)),
                read(header, "v2", 0, 599, "41M", true, acgt(41), quals(41, 20)));
    }

    @Test
    public void testUseEndInUnpairedReads() throws IOException {
        final SAMFileHeader header = defaultHeader();

        final FlowResult startOnly = runFlow(header, endInput(header));
        Assert.assertEquals(startOnly.duplicates, set("e1", "e2", "e4", "v2"));
        // e4 is soft clipped at its 3' end, which makes its fragment length known
        assertUnpairedMetrics(startOnly.flow("libA"), 6, 4, 1, 1, 3);

        final FlowResult withEnd = runFlow(header, endInput(header), "FLOW_USE_END_IN_UNPAIRED_READS=true");
        Assert.assertEquals(withEnd.duplicates, set("e2", "e4"));
        assertUnpairedMetrics(withEnd.flow("libA"), 6, 2, 1, 1, 1);
    }

    @Test
    public void testSpilledSortGivesTheSameDuplicates() throws IOException {
        // room for about two read ends before the fragment sort spills to disk
        final double ratio = 2.5 * ReadEndsForMarkDuplicates.getSizeOf() / Runtime.getRuntime().maxMemory();
        final SAMFileHeader header = defaultHeader();
        final FlowResult startOnly = runFlow(header, endInput(header), "SORTING_COLLECTION_SIZE_RATIO=" + ratio);
        Assert.assertEquals(startOnly.duplicates, set("e1", "e2", "e4", "v2"));
        final FlowResult withEnd = runFlow(header, endInput(header), "SORTING_COLLECTION_SIZE_RATIO=" + ratio,
                "FLOW_USE_END_IN_UNPAIRED_READS=true");
        Assert.assertEquals(withEnd.duplicates, set("e2", "e4"));
    }

    @Test
    public void testDuplicateTypeTagsAndRemoval() throws IOException {
        final SAMFileHeader header = defaultHeader();
        final Set<String> all = set("e1", "e2", "e3", "e4", "v1", "v2");
        final Set<String> duplicates = set("e1", "e2", "e4", "v2");

        final FlowResult tagged = runFlow(header, endInput(header), "TAGGING_POLICY=All");
        Assert.assertEquals(tagged.duplicates, duplicates);
        for (final String name : all) {
            Assert.assertEquals(tagged.records.get(name).getAttribute("DT"), duplicates.contains(name) ? "LB" : null, name);
        }

        final FlowResult removed = runFlow(header, endInput(header), "REMOVE_DUPLICATES=true");
        Assert.assertEquals(removed.names, set("e3", "v1"));
        assertUnpairedMetrics(removed.flow("libA"), 6, 4, 1, 1, 3);
    }

    @Test
    public void testInputsWithoutDuplicateCandidates() throws IOException {
        final SAMFileHeader header = defaultHeader();
        final List<SAMRecord> unmapped = new ArrayList<>();
        for (final String name : new String[]{"u1", "u2"}) {
            final SAMRecord rec = new SAMRecord(header);
            rec.setReadName(name);
            rec.setReadUnmappedFlag(true);
            rec.setReadBases(acgt(40).getBytes());
            rec.setBaseQualities(quals(40, 30));
            rec.setAttribute("RG", RG_A);
            unmapped.add(rec);
        }
        final FlowResult none = runFlow(header, unmapped, "FLOW_USE_END_IN_UNPAIRED_READS=true");
        Assert.assertEquals(none.names, set("u1", "u2"));
        Assert.assertTrue(none.duplicates.isEmpty());
        Assert.assertEquals(none.flow("libA").UNMAPPED_READS, 2);
        assertUnpairedMetrics(none.flow("libA"), 0, 0, 0, 0, 0);

        final FlowResult single = runFlow(header, Collections.singletonList(read(header, "only", 2200, 40, 30)),
                "FLOW_USE_END_IN_UNPAIRED_READS=true");
        Assert.assertTrue(single.duplicates.isEmpty());
        assertUnpairedMetrics(single.flow("libA"), 1, 0, 0, 0, 0);
    }

    @Test
    public void testStartAndEndUncertaintyTogether() throws IOException {
        final SAMFileHeader header = defaultHeader();
        // (start, end): n0 (2300, 2339), n1 (2302, 2341) within 2 of n0 at both ends, n2 (2303, 2350) only at the start
        final List<SAMRecord> records = Arrays.asList(
                read(header, "n0", 2300, 40, 40),
                read(header, "n1", 2302, 40, 30),
                read(header, "n2", 2303, 48, 30));
        final FlowResult result = runFlow(header, records, "FLOW_USE_END_IN_UNPAIRED_READS=true",
                "FLOW_UNPAIRED_START_UNCERTAINTY=2", "FLOW_UNPAIRED_END_UNCERTAINTY=2");
        Assert.assertEquals(result.duplicates, set("n1"));
    }

    @DataProvider(name = "endUncertainty")
    public Object[][] endUncertainty() {
        // reads start together; the first has the best score; ends are start + length - 1
        return new Object[][]{
                {2, new int[]{40, 42, 45}, set("n1")},
                {2, new int[]{40, 42, 44}, set("n1", "n2")},
                {2, new int[]{40, 44, 42}, set("n1", "n2")},
                {1, new int[]{40, 42}, set()},
                {0, new int[]{40, 40, 41}, set("n1")},
                // two chunks, 739-744 and 750-755; n3 outscores n2
                {5, new int[]{40, 45, 51, 56}, set("n1", "n2")},
        };
    }

    @Test(dataProvider = "endUncertainty")
    public void testUnpairedEndUncertainty(final int uncertainty, final int[] lengths, final Set<String> duplicates) throws IOException {
        final SAMFileHeader header = defaultHeader();
        final List<SAMRecord> records = new ArrayList<>();
        for (int i = 0; i < lengths.length; i++) {
            records.add(read(header, "n" + i, 700, lengths[i], i == 0 ? 40 : 30));
        }
        final FlowResult result = runFlow(header, records, "FLOW_USE_END_IN_UNPAIRED_READS=true",
                "FLOW_UNPAIRED_END_UNCERTAINTY=" + uncertainty);
        Assert.assertEquals(result.duplicates, duplicates);
    }

    @DataProvider(name = "startUncertainty")
    public Object[][] startUncertainty() {
        // the first read has the best score
        return new Object[][]{
                {2, new int[]{800, 802, 804, 807}, set("n1", "n2")},
                {2, new int[]{800, 803}, set()},
                {0, new int[]{800, 800, 801}, set("n1")},
                {3, new int[]{800, 803}, set("n1")},
        };
    }

    @Test(dataProvider = "startUncertainty")
    public void testUnpairedStartUncertainty(final int uncertainty, final int[] starts, final Set<String> duplicates) throws IOException {
        final SAMFileHeader header = defaultHeader();
        final List<SAMRecord> records = new ArrayList<>();
        for (int i = 0; i < starts.length; i++) {
            records.add(read(header, "n" + i, starts[i], 40, i == 0 ? 40 : 30));
        }
        final FlowResult result = runFlow(header, records, "FLOW_UNPAIRED_START_UNCERTAINTY=" + uncertainty);
        Assert.assertEquals(result.duplicates, duplicates);
    }

    @DataProvider(name = "clippingTag")
    public Object[][] clippingTag() {
        return new Object[][]{
                {new String[]{"FLOW_USE_END_IN_UNPAIRED_READS=true"}, set("z1")},
                {new String[]{"FLOW_USE_END_IN_UNPAIRED_READS=true", "FLOW_Q_IS_KNOWN_END=false"}, set("z1")},
                // without the end, z1 has the best score of the three reads sharing a start
                {new String[]{}, set("k1", "k2")},
        };
    }

    @Test(dataProvider = "clippingTag")
    public void testClippingTagDecidesWhetherTheEndIsKnown(final String[] args, final Set<String> duplicates) throws IOException {
        final SAMFileHeader header = defaultHeader();
        final SAMRecord k1 = read(header, "k1", 900, 40, 30);
        // adapter clipped: the end is known
        final SAMRecord k2 = read(header, "k2", 900, 45, 30);
        k2.setAttribute("tm", "A");
        // trimmed but not adapter clipped: the end is unknown, so it joins the first read with a known end
        final SAMRecord z1 = read(header, "z1", 900, 50, 40);
        z1.setAttribute("tm", "Z");
        final FlowResult result = runFlow(header, Arrays.asList(k1, k2, z1), args);
        Assert.assertEquals(result.duplicates, duplicates);
        final boolean useEnd = args.length > 0;
        assertUnpairedMetrics(result.flow("libA"), 3, duplicates.size(), 1, useEnd ? 0 : 1, 1);
    }

    /** Reads soft clipped at their 3' ends only; the argument is documented for the 3' end. */
    private static List<SAMRecord> clippedInput(final SAMFileHeader header) {
        return Arrays.asList(
                // aligned 1000-1039
                read(header, "c1", 0, 1000, "40M", false, acgt(40), quals(40, 40)),
                read(header, "c2", 0, 1000, "40M3S", false, acgt(43), quals(43, 30)),
                // aligned 1100-1139 on the reverse strand, whose 3' end is the alignment start
                read(header, "d1", 0, 1100, "40M", true, acgt(40), quals(40, 40)),
                read(header, "d2", 0, 1100, "2S40M", true, acgt(42), quals(42, 30)));
    }

    @DataProvider(name = "clippedEnd")
    public Object[][] clippedEnd() {
        return new Object[][]{
                {new String[]{"FLOW_USE_END_IN_UNPAIRED_READS=true", "FLOW_USE_UNPAIRED_CLIPPED_END=true"}, set("c2", "d2")},
                {new String[]{"FLOW_USE_END_IN_UNPAIRED_READS=true", "FLOW_USE_UNPAIRED_CLIPPED_END=false"}, set()},
                {new String[]{"FLOW_USE_UNPAIRED_CLIPPED_END=true"}, set("c2", "d2")},
                {new String[]{"FLOW_USE_UNPAIRED_CLIPPED_END=false"}, set("c2", "d2")},
        };
    }

    @Test(dataProvider = "clippedEnd")
    public void testUseUnpairedClippedEnd(final String[] args, final Set<String> duplicates) throws IOException {
        final SAMFileHeader header = defaultHeader();
        final FlowResult result = runFlow(header, clippedInput(header), args);
        Assert.assertEquals(result.duplicates, duplicates);
        // c2 and d2 are the reads with a known fragment length
        assertUnpairedMetrics(result.flow("libA"), 4, duplicates.size(), 2, duplicates.size(), 0);
    }

    @DataProvider(name = "skipFlows")
    public Object[][] skipFlows() {
        // s1 TTGCA... at 1200 and s2 TGCA... at 1201 are the same molecule with one extra T on s1;
        // s3 TCAG... and s4 TGCA... both start at 1300 and only share their first flow
        return new Object[][]{
                {0, "own", set("s4")},
                {1, "own", set("s2", "s4")},
                {2, "own", set("s2")},
                {3, "own", set("s2")},
                {2, "other", set("s2")},
                {2, "none", set("s4")},
                {1, "none", set("s4")},
        };
    }

    /**
     * When the reads' read group has no flow order, that of another read group is used (the "other" row). The
     * documentation does not specify this; the expected value follows the current implementation.
     */
    @Test(dataProvider = "skipFlows")
    public void testSkipFirstNFlows(final int skip, final String flowOrderSource, final Set<String> duplicates) throws IOException {
        final SAMFileHeader header;
        switch (flowOrderSource) {
            case "own":
                header = defaultHeader();
                break;
            case "other":
                // the reads' read group has no flow order; another read group's is used
                header = header(RG_A, "libA", null, "rgB", "libB", FLOW_ORDER);
                break;
            case "none":
                header = header(RG_A, "libA", null);
                break;
            default:
                throw new IllegalArgumentException(flowOrderSource);
        }
        final String filler = acgt(40).substring(1);
        final List<SAMRecord> records = Arrays.asList(
                read(header, "s1", 0, 1200, "40M", false, ("TTGCA" + filler).substring(0, 40), quals(40, 40)),
                read(header, "s2", 0, 1201, "39M", false, ("TGCA" + filler).substring(0, 39), quals(39, 30)),
                read(header, "s3", 0, 1300, "40M", false, ("TCAGT" + filler).substring(0, 40), quals(40, 40)),
                read(header, "s4", 0, 1300, "40M", false, ("TGCAG" + filler).substring(0, 40), quals(40, 30)));
        final FlowResult result = runFlow(header, records, "FLOW_SKIP_FIRST_N_FLOWS=" + skip);
        Assert.assertEquals(result.duplicates, duplicates);
    }

    @Test
    public void testVendorFailedReadIsNotKept() throws IOException {
        final SAMFileHeader header = defaultHeader();
        final SAMRecord failed = read(header, "failed", 1400, 40, 40);
        failed.setReadFailsVendorQualityCheckFlag(true);
        final SAMRecord passed = read(header, "passed", 1400, 40, 20);
        Assert.assertEquals(runFlow(header, Arrays.asList(failed, passed)).duplicates, set("failed"));
        Assert.assertEquals(runFlow(header, Arrays.asList(failed, passed), "FLOW_USE_END_IN_UNPAIRED_READS=true").duplicates,
                set("failed"));
    }

    @DataProvider(name = "barcodes")
    public Object[][] barcodes() {
        return new Object[][]{
                {null, "BC", new String[]{"AA", "AA", "CC"}, set("b2", "b3")},
                {"BARCODE_TAG=BC", "BC", new String[]{"AA", "AA", "CC"}, set("b2")},
                {"READ_ONE_BARCODE_TAG=RX", "RX", new String[]{"AA", "CC", "AA"}, set("b3")},
        };
    }

    @Test(dataProvider = "barcodes")
    public void testBarcodes(final String arg, final String tag, final String[] values, final Set<String> duplicates) throws IOException {
        final SAMFileHeader header = defaultHeader();
        final List<SAMRecord> records = new ArrayList<>();
        for (int i = 0; i < values.length; i++) {
            final SAMRecord rec = read(header, "b" + (i + 1), 1500, 40, i == 0 ? 40 : 30);
            rec.setAttribute(tag, values[i]);
            records.add(rec);
        }
        for (final String end : new String[]{"FLOW_USE_END_IN_UNPAIRED_READS=false", "FLOW_USE_END_IN_UNPAIRED_READS=true"}) {
            final FlowResult result = arg == null ? runFlow(header, records, end) : runFlow(header, records, end, arg);
            Assert.assertEquals(result.duplicates, duplicates, end);
        }
    }

    @Test
    public void testLibrariesContigsAndStrandsAreKeptApart() throws IOException {
        final SAMFileHeader header = header(RG_A, "libA", FLOW_ORDER, "rgB", "libB", FLOW_ORDER);
        final SAMRecord other = read(header, "l3", 1600, 40, 30);
        other.setAttribute("RG", "rgB");
        final List<SAMRecord> records = Arrays.asList(
                read(header, "l1", 1600, 40, 40),
                read(header, "l2", 1600, 40, 30),
                other,
                read(header, "l4", 1, 1600, "40M", false, acgt(40), quals(40, 30)),
                read(header, "l5", 0, 1561, "40M", true, acgt(40), quals(40, 30)));
        final FlowResult result = runFlow(header, records, "FLOW_USE_END_IN_UNPAIRED_READS=true", "FLOW_UNPAIRED_START_UNCERTAINTY=1");
        Assert.assertEquals(result.duplicates, set("l2"));
        assertUnpairedMetrics(result.flow("libA"), 4, 1, 0, 0, 1);
        assertUnpairedMetrics(result.flow("libB"), 1, 0, 0, 0, 0);
    }

    @Test
    public void testPairedReadsWithUnmappedMatesAreTreatedAsFragments() throws IOException {
        final SAMFileHeader header = defaultHeader();
        final List<SAMRecord> records = new ArrayList<>();
        for (final String name : new String[]{"p1", "p2"}) {
            final SAMRecord mapped = read(header, name, 1700, 40, name.equals("p1") ? 40 : 30);
            mapped.setReadPairedFlag(true);
            mapped.setFirstOfPairFlag(true);
            mapped.setMateUnmappedFlag(true);
            mapped.setMateReferenceIndex(0);
            mapped.setMateAlignmentStart(1700);
            final SAMRecord unmapped = new SAMRecord(header);
            unmapped.setReadName(name);
            unmapped.setReadPairedFlag(true);
            unmapped.setSecondOfPairFlag(true);
            unmapped.setReadUnmappedFlag(true);
            unmapped.setReferenceIndex(0);
            unmapped.setAlignmentStart(1700);
            unmapped.setMateReferenceIndex(0);
            unmapped.setMateAlignmentStart(1700);
            unmapped.setReadBases(acgt(40).getBytes());
            unmapped.setBaseQualities(quals(40, 30));
            unmapped.setAttribute("RG", RG_A);
            records.add(mapped);
            records.add(unmapped);
        }
        final FlowResult result = runFlow(header, records);
        Assert.assertEquals(result.duplicates, set("p2"));
        final FlowBasedDuplicationMetrics m = result.flow("libA");
        Assert.assertEquals(m.UNMAPPED_READS, 2);
        assertUnpairedMetrics(m, 2, 1, 0, 0, 1);
    }

    private static File writeOne(final SAMFileHeader header, final SAMRecord rec) throws IOException {
        return writeSorted(header, Collections.singletonList(rec));
    }

    private void assertFlowRunThrows(final File input, final String message, final String... extraArgs) throws IOException {
        final File output = File.createTempFile("mdFlow.", ".sam");
        output.deleteOnExit();
        final File metrics = File.createTempFile("mdFlow.", ".duplicate_metrics");
        metrics.deleteOnExit();
        final List<String> args = new ArrayList<>(Arrays.asList("INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + output.getAbsolutePath(), "METRICS_FILE=" + metrics.getAbsolutePath(), "FLOW_MODE=true"));
        args.addAll(Arrays.asList(extraArgs));
        try {
            runPicardCommandLine(args);
            Assert.fail("expected an IllegalArgumentException");
        } catch (final IllegalArgumentException e) {
            Assert.assertTrue(e.getMessage().contains(message), e.getMessage());
        }
    }

    @Test
    public void testPairedReadsWithMappedMatesAreRejected() throws IOException {
        final SAMFileHeader header = defaultHeader();
        final SAMRecord rec = read(header, "paired", 1800, 40, 30);
        rec.setReadPairedFlag(true);
        rec.setFirstOfPairFlag(true);
        rec.setMateReferenceIndex(0);
        rec.setMateAlignmentStart(1900);
        rec.setMateNegativeStrandFlag(true);
        assertFlowRunThrows(writeOne(header, rec), "paired reads");
    }

    @DataProvider(name = "endUncertaintyWithoutEnd")
    public Object[][] endUncertaintyWithoutEnd() {
        return new Object[][]{{"FLOW_UNPAIRED_END_UNCERTAINTY=1"}, {"FLOW_UNPAIRED_END_UNCERTAINTY=-1"}};
    }

    @Test(dataProvider = "endUncertaintyWithoutEnd")
    public void testEndUncertaintyRequiresUseEnd(final String arg) throws IOException {
        final SAMFileHeader header = defaultHeader();
        assertFlowRunThrows(writeOne(header, read(header, "a", 1900, 40, 30)),
                "UNPAIRED_END_UNCERTAINTY", arg);
    }

    @Test
    public void testFlowArgumentsWithoutFlowModeAreIgnored() throws IOException {
        final SAMFileHeader header = defaultHeader();
        final List<SAMRecord> records = Arrays.asList(read(header, "a", 2000, 40, 30), read(header, "b", 2000, 41, 20));
        final FlowResult result = run(header, records, "FLOW_UNPAIRED_END_UNCERTAINTY=1", "FLOW_USE_UNPAIRED_CLIPPED_END=true");
        Assert.assertEquals(result.duplicates, set("b"));
        final DuplicationMetrics m = result.metricsByLibrary.get("libA");
        Assert.assertEquals(m.getClass(), DuplicationMetrics.class);
        Assert.assertEquals(m.UNPAIRED_READS_EXAMINED, 2);
        Assert.assertEquals(m.UNPAIRED_READ_DUPLICATES, 1);
    }

    @Test
    public void testUnknownStrategyFailsTheCommandLine() throws IOException {
        final SAMFileHeader header = defaultHeader();
        final File input = writeOne(header, read(header, "a", 2100, 40, 30));
        final File output = File.createTempFile("mdFlow.", ".sam");
        output.deleteOnExit();
        final File metrics = File.createTempFile("mdFlow.", ".duplicate_metrics");
        metrics.deleteOnExit();
        Assert.assertEquals(runPicardCommandLine(Arrays.asList("INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + output.getAbsolutePath(), "METRICS_FILE=" + metrics.getAbsolutePath(),
                "FLOW_MODE=true", "FLOW_DUP_STRATEGY=FLOW_RANDOM_STRATEGY")), 1);
    }
}
