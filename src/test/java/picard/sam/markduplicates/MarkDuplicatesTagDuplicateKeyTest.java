/*
 * The MIT License
 *
 * Copyright (c) 2025 The Broad Institute
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
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Tests for the opt-in {@code TAG_DUPLICATE_KEY} argument, which writes the canonical fragment duplicate
 * key to the {@code kf} attribute on every primary, mapped record.
 * <p>
 * The load-bearing property is that the emitted keys reflect the key MarkDuplicates actually marks on: the
 * partition of templates induced by their {@code kf} key(s) must be identical to the partition induced by
 * the {@code DI} (duplicate-set-index) tag. If that holds, an external tool can recover Picard's exact
 * duplicate sets from the tags without re-deriving Picard's clipping/coordinate conventions.
 */
public class MarkDuplicatesTagDuplicateKeyTest {

    private static final int READ_LENGTH = 36;
    private static final int BASE_QUALITY = 30;
    /** The SAM tag used to carry UMIs in the barcode-aware tests. */
    private static final String UMI_TAG = "RX";

    /** Builds a coordinate-sorted SAMRecordSetBuilder ready for adding reads. */
    private SAMRecordSetBuilder newBuilder() {
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
        builder.setReadLength(READ_LENGTH);
        return builder;
    }

    /** Materializes the builder's (sorted) records into a temp BAM file. */
    private File writeToBam(final SAMRecordSetBuilder builder) throws IOException {
        return writeToBam(builder, null);
    }

    /**
     * Materializes the builder's (sorted) records into a temp BAM file, optionally stamping each record with a
     * UMI in the {@link #UMI_TAG} attribute looked up by read name (both ends of a template get the same UMI).
     */
    private File writeToBam(final SAMRecordSetBuilder builder, final Map<String, String> umiByName) throws IOException {
        final File bam = File.createTempFile("tag_dup_key.input.", ".bam");
        bam.deleteOnExit();
        try (final SamReader reader = builder.getSamReader();
             final SAMFileWriter writer = new SAMFileWriterFactory().makeWriter(reader.getFileHeader(), true, bam, null)) {
            for (final SAMRecord rec : reader) {
                if (umiByName != null) {
                    final String umi = umiByName.get(rec.getReadName());
                    if (umi != null) {
                        rec.setAttribute(UMI_TAG, umi);
                    }
                }
                writer.addAlignment(rec);
            }
        }
        return bam;
    }

    /** Runs MarkDuplicates on the input, returning the output BAM. DI tagging is always on; kf tagging is optional. */
    private File runMarkDuplicates(final File input, final boolean tagDuplicateKey) throws IOException {
        return runMarkDuplicates(input, tagDuplicateKey, null);
    }

    /**
     * Runs MarkDuplicates on the input, returning the output BAM. DI tagging is always on; kf tagging is
     * optional. When {@code barcodeTag} is non-null it is passed as BARCODE_TAG so duplicate marking (and the
     * kf key) become UMI-aware.
     */
    private File runMarkDuplicates(final File input, final boolean tagDuplicateKey, final String barcodeTag) throws IOException {
        final File output = File.createTempFile("tag_dup_key.output.", ".bam");
        final File metrics = File.createTempFile("tag_dup_key.metrics.", ".txt");
        output.deleteOnExit();
        metrics.deleteOnExit();

        final List<String> args = new ArrayList<>();
        args.add("INPUT=" + input.getAbsolutePath());
        args.add("OUTPUT=" + output.getAbsolutePath());
        args.add("METRICS_FILE=" + metrics.getAbsolutePath());
        // PROGRAM_RECORD_ID=null avoids the jar-manifest version lookup that fails when running outside a jar.
        args.add("PROGRAM_RECORD_ID=null");
        args.add("TAG_DUPLICATE_SET_MEMBERS=true");
        if (tagDuplicateKey) {
            args.add("TAG_DUPLICATE_KEY=true");
        }
        if (barcodeTag != null) {
            args.add("BARCODE_TAG=" + barcodeTag);
        }
        Assert.assertEquals(new MarkDuplicates().instanceMain(args.toArray(new String[0])), 0);
        return output;
    }

    /** Reads all records from a BAM into a list. */
    private List<SAMRecord> readAll(final File bam) {
        final List<SAMRecord> records = new ArrayList<>();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(bam)) {
            for (final SAMRecord rec : reader) {
                records.add(rec);
            }
        } catch (final IOException e) {
            throw new RuntimeException(e);
        }
        return records;
    }

    /**
     * Computes a template's group key from its primary, mapped records' {@code kf} tags. For a pair this is
     * the (sorted) combination of the two ends' fragment keys, which is a faithful proxy for Picard's
     * canonical pair key; for an orphan it is the single fragment key.
     */
    private String groupKeyFor(final List<String> fragmentKeys) {
        final List<String> sorted = new ArrayList<>(fragmentKeys);
        Collections.sort(sorted);
        return String.join("|", sorted);
    }

    /**
     * The headline assertion: the partition of templates induced by their kf key(s) equals the partition
     * induced by DI. Asserts over the output of a run that has both tags enabled.
     */
    private void assertKeyPartitionEqualsDuplicateIndexPartition(final File output) {
        final List<SAMRecord> records = readAll(output);

        // Collect, per read name (template), its fragment keys and its DI (if any).
        final Map<String, List<String>> fragmentKeysByName = new HashMap<>();
        final Map<String, Integer> diByName = new HashMap<>();
        for (final SAMRecord rec : records) {
            if (rec.getReadUnmappedFlag() || rec.isSecondaryOrSupplementary()) {
                continue;
            }
            final String kf = rec.getStringAttribute(MarkDuplicates.DUPLICATE_KEY_TAG);
            Assert.assertNotNull(kf, "every primary mapped record must carry a kf tag: " + rec.getReadName());
            fragmentKeysByName.computeIfAbsent(rec.getReadName(), k -> new ArrayList<>()).add(kf);

            final Integer di = rec.getIntegerAttribute(MarkDuplicates.DUPLICATE_SET_INDEX_TAG);
            if (di != null) {
                diByName.put(rec.getReadName(), di);
            }
        }

        // Partition templates by group key (from kf) and by DI.
        final Map<String, Set<String>> namesByGroupKey = new HashMap<>();
        for (final Map.Entry<String, List<String>> e : fragmentKeysByName.entrySet()) {
            namesByGroupKey.computeIfAbsent(groupKeyFor(e.getValue()), k -> new HashSet<>()).add(e.getKey());
        }
        final Map<Integer, Set<String>> namesByDuplicateIndex = new HashMap<>();
        for (final Map.Entry<String, Integer> e : diByName.entrySet()) {
            namesByDuplicateIndex.computeIfAbsent(e.getValue(), k -> new HashSet<>()).add(e.getKey());
        }

        // DI only appears on members of a duplicate set (size >= 2); the group-key partition additionally
        // contains the singletons. So compare the DI partition against the multi-member group-key classes.
        final Set<Set<String>> duplicateIndexPartition = new HashSet<>(namesByDuplicateIndex.values());
        final Set<Set<String>> multiMemberGroupKeyPartition = new HashSet<>();
        for (final Set<String> names : namesByGroupKey.values()) {
            if (names.size() >= 2) {
                multiMemberGroupKeyPartition.add(names);
            }
        }
        Assert.assertEquals(multiMemberGroupKeyPartition, duplicateIndexPartition,
                "kf-induced duplicate sets must equal DI-induced duplicate sets");
    }

    /** Maps each template's read name to the kf key of its first-of-pair (read one) primary, mapped record. */
    private Map<String, String> read1FragmentKeyByName(final File output) {
        final Map<String, String> keys = new HashMap<>();
        for (final SAMRecord rec : readAll(output)) {
            if (rec.getReadUnmappedFlag() || rec.isSecondaryOrSupplementary()) {
                continue;
            }
            if (rec.getReadPairedFlag() && !rec.getFirstOfPairFlag()) {
                continue;
            }
            keys.put(rec.getReadName(), rec.getStringAttribute(MarkDuplicates.DUPLICATE_KEY_TAG));
        }
        return keys;
    }

    @Test
    public void kfKeysInduceTheSameDuplicateSetsAsDuplicateIndex() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        // Duplicate set of 2 on chr1.
        builder.addPair("dupA_rep", 0, 1000, 1200, false, false, "36M", "36M", false, true, BASE_QUALITY);
        builder.addPair("dupA_dup", 0, 1000, 1200, false, false, "36M", "36M", false, true, BASE_QUALITY);
        // Duplicate set of 3 on chr2.
        builder.addPair("dupB_rep", 1, 5000, 5300, false, false, "36M", "36M", false, true, BASE_QUALITY);
        builder.addPair("dupB_dup1", 1, 5000, 5300, false, false, "36M", "36M", false, true, BASE_QUALITY);
        builder.addPair("dupB_dup2", 1, 5000, 5300, false, false, "36M", "36M", false, true, BASE_QUALITY);
        // A unique pair (singleton): must still carry kf, but no DI.
        builder.addPair("unique", 0, 8000, 8200, false, false, "36M", "36M", false, true, BASE_QUALITY);

        final File output = runMarkDuplicates(writeToBam(builder), true);
        assertKeyPartitionEqualsDuplicateIndexPartition(output);

        // The unique pair has no DI but must still carry kf on both ends.
        int uniqueKfCount = 0;
        for (final SAMRecord rec : readAll(output)) {
            if (rec.getReadName().equals("unique")) {
                Assert.assertNull(rec.getIntegerAttribute(MarkDuplicates.DUPLICATE_SET_INDEX_TAG));
                Assert.assertNotNull(rec.getStringAttribute(MarkDuplicates.DUPLICATE_KEY_TAG));
                uniqueKfCount++;
            }
        }
        Assert.assertEquals(uniqueKfCount, 2, "singleton pair should have kf on both ends");
    }

    @Test
    public void kfPartitionMatchesDuplicateIndexAcrossOrientationsAndChromosomes() throws IOException {
        // The pair-key proxy (sorted kf strings) must reproduce Picard's duplicate
        // sets even where the canonicalization is non-trivial: differing per-end
        // strand/orientation, FF/RR orientations, and inter-chromosomal pairs. The
        // assertion verifies the kf-induced partition equals the DI partition, so a
        // wrong merge (e.g. ignoring orientation) or split would fail it.
        final SAMRecordSetBuilder builder = newBuilder();

        // FR dup set (read1 forward, read2 reverse) on chr1.
        builder.addPair("fr_rep", 0, 1000, 1500, false, false, "36M", "36M", false, true, BASE_QUALITY);
        builder.addPair("fr_dup", 0, 1000, 1500, false, false, "36M", "36M", false, true, BASE_QUALITY);

        // RF dup set (read1 reverse, read2 forward) at the SAME starts: a different
        // orientation, so Picard keeps it as a separate set — kf must NOT merge it
        // with the FR set above.
        builder.addPair("rf_rep", 0, 1000, 1500, false, false, "36M", "36M", true, false, BASE_QUALITY);
        builder.addPair("rf_dup", 0, 1000, 1500, false, false, "36M", "36M", true, false, BASE_QUALITY);

        // FF dup set (both forward) on chr2.
        builder.addPair("ff_rep", 1, 2000, 2400, false, false, "36M", "36M", false, false, BASE_QUALITY);
        builder.addPair("ff_dup", 1, 2000, 2400, false, false, "36M", "36M", false, false, BASE_QUALITY);

        // Inter-chromosomal dup set: read1 on chr1, read2 on chr2 (14-arg addPair
        // with two reference indexes).
        builder.addPair("ic_rep", 0, 1, 7000, 8000, false, false, "36M", "36M", false, true, false, false, BASE_QUALITY);
        builder.addPair("ic_dup", 0, 1, 7000, 8000, false, false, "36M", "36M", false, true, false, false, BASE_QUALITY);

        // Two genuinely different keys (read1 at different positions) must NOT merge.
        builder.addPair("uniqA", 0, 9000, 9300, false, false, "36M", "36M", false, true, BASE_QUALITY);
        builder.addPair("uniqB", 0, 9100, 9300, false, false, "36M", "36M", false, true, BASE_QUALITY);

        final File output = runMarkDuplicates(writeToBam(builder), true);
        assertKeyPartitionEqualsDuplicateIndexPartition(output);

        // Explicit over-merge guard: the FR and RF sets at the same starts must
        // carry distinct kf-derived pair keys (different orientation -> different key).
        final java.util.Map<String, String> pairKeyByName = new HashMap<>();
        for (final SAMRecord rec : readAll(output)) {
            if (rec.getReadUnmappedFlag() || rec.isSecondaryOrSupplementary()) {
                continue;
            }
            final String kf = rec.getStringAttribute(MarkDuplicates.DUPLICATE_KEY_TAG);
            pairKeyByName.merge(rec.getReadName(), kf, (a, b) -> {
                final String[] ends = {a, b};
                java.util.Arrays.sort(ends);
                return ends[0] + "|" + ends[1];
            });
        }
        Assert.assertNotEquals(
                pairKeyByName.get("fr_rep"), pairKeyByName.get("rf_rep"),
                "FR and RF pairs at the same starts must not share a kf pair key");
    }

    @Test
    public void orphanKeyMatchesPairEndKeyAndOrphanIsMarkedDuplicate() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        // A pair whose forward end starts at chr1:2000.
        builder.addPair("pair", 0, 2000, 2300, false, false, "36M", "36M", false, true, BASE_QUALITY);
        // A mapped orphan (mate unmapped) at the same forward position: its fragment key must equal the
        // pair's forward-end key, and Picard must mark the orphan as a duplicate ("fragments don't beat pairs").
        builder.addPair("orphan", 0, 2000, 2000, false, true, "36M", "36M", false, false, BASE_QUALITY);

        final File output = runMarkDuplicates(writeToBam(builder), true);
        final List<SAMRecord> records = readAll(output);

        final Set<String> pairKeys = new HashSet<>();
        String orphanKey = null;
        boolean orphanMarkedDuplicate = false;
        for (final SAMRecord rec : records) {
            if (rec.getReadUnmappedFlag() || rec.isSecondaryOrSupplementary()) {
                continue;
            }
            final String kf = rec.getStringAttribute(MarkDuplicates.DUPLICATE_KEY_TAG);
            if (rec.getReadName().equals("pair")) {
                pairKeys.add(kf);
            } else if (rec.getReadName().equals("orphan")) {
                orphanKey = kf;
                orphanMarkedDuplicate = rec.getDuplicateReadFlag();
            }
        }
        Assert.assertNotNull(orphanKey);
        Assert.assertTrue(pairKeys.contains(orphanKey),
                "orphan fragment key " + orphanKey + " should match a pair end key " + pairKeys);
        Assert.assertTrue(orphanMarkedDuplicate, "orphan should be marked duplicate as a fragment dup of the pair");
    }

    @Test
    public void umiDuplicateSetsMatchDuplicateIndexWhenBarcodeTagSet() throws IOException {
        // Three pairs at identical coordinates. Under UMI-aware marking, only the two sharing a UMI are
        // duplicates; the third (different UMI) stays a singleton despite the identical position. This is the
        // exact case a barcode-less key got wrong, so the partition-equality assertion is the regression guard.
        final SAMRecordSetBuilder builder = newBuilder();
        builder.addPair("umiA_rep", 0, 1000, 1200, false, false, "36M", "36M", false, true, BASE_QUALITY);
        builder.addPair("umiA_dup", 0, 1000, 1200, false, false, "36M", "36M", false, true, BASE_QUALITY);
        builder.addPair("umiB",     0, 1000, 1200, false, false, "36M", "36M", false, true, BASE_QUALITY);

        final Map<String, String> umiByName = new HashMap<>();
        umiByName.put("umiA_rep", "AAAAAAAA");
        umiByName.put("umiA_dup", "AAAAAAAA");
        umiByName.put("umiB", "CCCCCCCC");

        final File output = runMarkDuplicates(writeToBam(builder, umiByName), true, UMI_TAG);
        assertKeyPartitionEqualsDuplicateIndexPartition(output);

        // umiB shares its position with the umiA set but, having a distinct UMI, must remain a singleton (no DI).
        for (final SAMRecord rec : readAll(output)) {
            if (rec.getReadName().equals("umiB")) {
                Assert.assertNull(rec.getIntegerAttribute(MarkDuplicates.DUPLICATE_SET_INDEX_TAG),
                        "umiB has a distinct UMI and must not join the umiA duplicate set");
            }
        }
    }

    @Test
    public void kfIncludesBarcodeFieldsOnlyWhenBarcodeTagSet() throws IOException {
        // Two pairs at identical coordinates but with different UMIs.
        final SAMRecordSetBuilder builder = newBuilder();
        builder.addPair("p1", 0, 1000, 1200, false, false, "36M", "36M", false, true, BASE_QUALITY);
        builder.addPair("p2", 0, 1000, 1200, false, false, "36M", "36M", false, true, BASE_QUALITY);
        final Map<String, String> umiByName = new HashMap<>();
        umiByName.put("p1", "AAAAAAAA");
        umiByName.put("p2", "CCCCCCCC");
        final File input = writeToBam(builder, umiByName);

        // With BARCODE_TAG: the readable UMI is appended as a single positional field (lib=..;..:..:strand:<umi>),
        // so the key carries the verbatim UMI and the two different UMIs produce different keys.
        final Map<String, String> withBarcode = read1FragmentKeyByName(runMarkDuplicates(input, true, UMI_TAG));
        Assert.assertTrue(withBarcode.get("p1").endsWith(":AAAAAAAA"),
                "kf must append the verbatim UMI when BARCODE_TAG is set: " + withBarcode.get("p1"));
        Assert.assertEquals(fieldsAfterLibrary(withBarcode.get("p1")), 4,
                "barcode-aware key must have referenceIndex:coord:strand plus the joined barcode field");
        Assert.assertNotEquals(withBarcode.get("p1"), withBarcode.get("p2"),
                "reads at the same position with different UMIs must get different kf keys");

        // Without BARCODE_TAG: no barcode fields (backward compatible) and the position-identical ends share a key.
        final Map<String, String> noBarcode = read1FragmentKeyByName(runMarkDuplicates(input, true));
        Assert.assertEquals(fieldsAfterLibrary(noBarcode.get("p1")), 3,
                "barcode-free key must be just referenceIndex:coord:strand: " + noBarcode.get("p1"));
        Assert.assertFalse(noBarcode.get("p1").contains("AAAAAAAA"),
                "kf must not include the UMI when no BARCODE_TAG is set: " + noBarcode.get("p1"));
        Assert.assertEquals(noBarcode.get("p1"), noBarcode.get("p2"),
                "without UMIs, identical-position read-one ends must share a kf key");
    }

    /** Counts the colon-separated fields after the {@code lib=...;} prefix of a kf key. */
    private int fieldsAfterLibrary(final String key) {
        final String afterLibrary = key.substring(key.indexOf(';') + 1);
        return afterLibrary.split(":", -1).length;
    }

    @Test
    public void kfJoinsMultipleBarcodeSequencesWithHyphens() throws IOException {
        // With both a UMI (RX) and a read-one barcode (BC) present, the two non-empty sequences are joined in
        // key order with a single hyphen ('<umi>-<readOneBarcode>'), with no extra hyphens or colons.
        final SAMRecordSetBuilder builder = newBuilder();
        builder.addPair("p", 0, 1000, 1200, false, false, "36M", "36M", false, true, BASE_QUALITY);

        final File input = File.createTempFile("tag_dup_key.input.", ".bam");
        input.deleteOnExit();
        try (final SamReader reader = builder.getSamReader();
             final SAMFileWriter writer = new SAMFileWriterFactory().makeWriter(reader.getFileHeader(), true, input, null)) {
            for (final SAMRecord rec : reader) {
                rec.setAttribute("RX", "ACTG");
                rec.setAttribute("BC", "TTGA");
                writer.addAlignment(rec);
            }
        }

        final File output = File.createTempFile("tag_dup_key.output.", ".bam");
        final File metrics = File.createTempFile("tag_dup_key.metrics.", ".txt");
        output.deleteOnExit();
        metrics.deleteOnExit();
        Assert.assertEquals(new MarkDuplicates().instanceMain(new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + output.getAbsolutePath(),
                "METRICS_FILE=" + metrics.getAbsolutePath(),
                "PROGRAM_RECORD_ID=null",
                "TAG_DUPLICATE_KEY=true",
                "BARCODE_TAG=RX",
                "READ_ONE_BARCODE_TAG=BC"}), 0);

        // Read one carries both the UMI and the read-one barcode -> 'ACTG-TTGA'.
        final String read1Key = read1FragmentKeyByName(output).get("p");
        Assert.assertTrue(read1Key.endsWith(":ACTG-TTGA"),
                "read-one key must join UMI and read-one barcode with a single hyphen: " + read1Key);
        Assert.assertEquals(fieldsAfterLibrary(read1Key), 4,
                "the joined barcodes must be a single field, not split across colons: " + read1Key);
    }

    @Test
    public void readTwoBarcodeIsRenderedAndDistinguishesDuplicates() throws IOException {
        // Three pairs at identical coordinates, distinguished only by a read-two barcode (READ_TWO_BARCODE_TAG),
        // which populates the readTwoBarcode field on second-of-pair records. Two pairs share the barcode (so
        // they are duplicates), the third differs (singleton); the kf-induced partition must still match DI.
        final SAMRecordSetBuilder builder = newBuilder();
        builder.addPair("pA", 0, 1000, 1200, false, false, "36M", "36M", false, true, BASE_QUALITY);
        builder.addPair("pB", 0, 1000, 1200, false, false, "36M", "36M", false, true, BASE_QUALITY);
        builder.addPair("pC", 0, 1000, 1200, false, false, "36M", "36M", false, true, BASE_QUALITY);

        final Map<String, String> readTwoBarcodeByName = new HashMap<>();
        readTwoBarcodeByName.put("pA", "ACAC");
        readTwoBarcodeByName.put("pB", "ACAC");
        readTwoBarcodeByName.put("pC", "GTGT");

        final File input = File.createTempFile("tag_dup_key.input.", ".bam");
        input.deleteOnExit();
        try (final SamReader reader = builder.getSamReader();
             final SAMFileWriter writer = new SAMFileWriterFactory().makeWriter(reader.getFileHeader(), true, input, null)) {
            for (final SAMRecord rec : reader) {
                // READ_TWO_BARCODE_TAG applies only to the second-of-pair read.
                if (rec.getReadPairedFlag() && !rec.getFirstOfPairFlag()) {
                    rec.setAttribute("BC", readTwoBarcodeByName.get(rec.getReadName()));
                }
                writer.addAlignment(rec);
            }
        }

        final File output = File.createTempFile("tag_dup_key.output.", ".bam");
        final File metrics = File.createTempFile("tag_dup_key.metrics.", ".txt");
        output.deleteOnExit();
        metrics.deleteOnExit();
        Assert.assertEquals(new MarkDuplicates().instanceMain(new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + output.getAbsolutePath(),
                "METRICS_FILE=" + metrics.getAbsolutePath(),
                "PROGRAM_RECORD_ID=null",
                "TAG_DUPLICATE_SET_MEMBERS=true",
                "TAG_DUPLICATE_KEY=true",
                "READ_TWO_BARCODE_TAG=BC"}), 0);

        assertKeyPartitionEqualsDuplicateIndexPartition(output);

        // The read-two end carries its barcode in the key, and pC's distinct barcode keeps it out of the pA/pB set.
        for (final SAMRecord rec : readAll(output)) {
            if (rec.getReadUnmappedFlag() || rec.isSecondaryOrSupplementary()) {
                continue;
            }
            final boolean isReadTwo = rec.getReadPairedFlag() && !rec.getFirstOfPairFlag();
            if (isReadTwo && rec.getReadName().equals("pA")) {
                Assert.assertTrue(rec.getStringAttribute(MarkDuplicates.DUPLICATE_KEY_TAG).endsWith(":ACAC"),
                        "read-two end must carry its barcode in the key: " + rec.getStringAttribute(MarkDuplicates.DUPLICATE_KEY_TAG));
            }
            if (rec.getReadName().equals("pC")) {
                Assert.assertNull(rec.getIntegerAttribute(MarkDuplicates.DUPLICATE_SET_INDEX_TAG),
                        "pC has a distinct read-two barcode and must not join the pA/pB duplicate set");
            }
        }
    }

    @Test
    public void tagDuplicateKeyIsRejectedWithFlowMode() throws IOException {
        // The fragment key cannot faithfully represent flow mode's position-uncertain matching, so combining
        // the two options must fail command-line validation (non-zero exit) rather than emit a misleading key.
        final SAMRecordSetBuilder builder = newBuilder();
        builder.addPair("a", 0, 1000, 1200, false, false, "36M", "36M", false, true, BASE_QUALITY);
        final File input = writeToBam(builder);

        final File output = File.createTempFile("tag_dup_key.output.", ".bam");
        final File metrics = File.createTempFile("tag_dup_key.metrics.", ".txt");
        output.deleteOnExit();
        metrics.deleteOnExit();

        final int returnCode = new MarkDuplicates().instanceMain(new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + output.getAbsolutePath(),
                "METRICS_FILE=" + metrics.getAbsolutePath(),
                "PROGRAM_RECORD_ID=null",
                "TAG_DUPLICATE_KEY=true",
                "FLOW_MODE=true"});
        Assert.assertNotEquals(returnCode, 0, "TAG_DUPLICATE_KEY combined with FLOW_MODE must fail validation");
    }

    @Test
    public void noKeyTagWhenArgumentIsOff() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        builder.addPair("a", 0, 1000, 1200, false, false, "36M", "36M", false, true, BASE_QUALITY);
        builder.addPair("b", 0, 1000, 1200, false, false, "36M", "36M", false, true, BASE_QUALITY);

        final File output = runMarkDuplicates(writeToBam(builder), false);
        for (final SAMRecord rec : readAll(output)) {
            Assert.assertNull(rec.getStringAttribute(MarkDuplicates.DUPLICATE_KEY_TAG),
                    "kf must not be written when TAG_DUPLICATE_KEY is false");
        }
    }
}
