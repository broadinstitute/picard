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

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.HashMap;
import java.util.Map;

/**
 * Tests for the two {@link ReadEndsForMarkDuplicatesMap} implementations, {@link MemoryBasedReadEndsForMarkDuplicatesMap}
 * and {@link DiskBasedReadEndsForMarkDuplicatesMap}: put, remove, size and sizeInRam against a reference model,
 * with records stored for more mate sequences than the disk-based map may keep files open for, each codec the
 * map is built with in MarkDuplicates, sequences revisited after being spilled, and keys that are absent.
 */
public class ReadEndsForMarkDuplicatesMapTest {

    private enum Kind {MEMORY, DISK, DISK_FORCED, DISK_BARCODES}

    @DataProvider(name = "maps")
    public Object[][] maps() {
        return new Object[][]{
                {Kind.MEMORY, 0},
                {Kind.DISK, 1},
                {Kind.DISK, 2},
                {Kind.DISK, 100},
                {Kind.DISK_FORCED, 2},
                {Kind.DISK_BARCODES, 2}
        };
    }

    private static ReadEndsForMarkDuplicatesMap newMap(final Kind kind, final int maxOpenFiles) {
        switch (kind) {
            case MEMORY:
                return new MemoryBasedReadEndsForMarkDuplicatesMap();
            case DISK:
                return new DiskBasedReadEndsForMarkDuplicatesMap(maxOpenFiles, new ReadEndsForMarkDuplicatesCodec());
            case DISK_FORCED:
                final ReadEndsForMarkDuplicatesCodec forced = new ReadEndsForMarkDuplicatesCodec();
                forced.setForceSerialization();
                return new DiskBasedReadEndsForMarkDuplicatesMap(maxOpenFiles, forced);
            case DISK_BARCODES:
                return new DiskBasedReadEndsForMarkDuplicatesMap(maxOpenFiles, new ReadEndsForMarkDuplicatesWithBarcodesCodec());
            default:
                throw new IllegalArgumentException(kind.toString());
        }
    }

    /**
     * A record as MarkDuplicates stores it while waiting for the mate: a fragment orientation, with the mate's
     * coordinate and file index unset unless the codec keeps them (flow mode).
     */
    private static ReadEndsForMarkDuplicates ends(final Kind kind, final int id, final int mateSequence) {
        final ReadEndsForMarkDuplicates ends = kind == Kind.DISK_BARCODES ? new ReadEndsForMarkDuplicatesWithBarcodes() : new ReadEndsForMarkDuplicates();
        ends.libraryId = (short) (id % 3 + 1);
        ends.orientation = id % 2 == 0 ? ReadEnds.F : ReadEnds.R;
        ends.read1ReferenceIndex = 0;
        ends.read1Coordinate = 100 * id + 1;
        ends.read1IndexInFile = 1000L * id;
        ends.read2ReferenceIndex = mateSequence;
        if (kind == Kind.DISK_FORCED || kind == Kind.MEMORY) {
            ends.read2Coordinate = 100 * id + 51;
            ends.read2IndexInFile = 1000L * id + 1;
        }
        ends.readGroup = (short) (id % 2);
        ends.tile = (short) (1100 + id);
        ends.x = 10 * id;
        ends.y = Integer.MAX_VALUE - id;
        ends.orientationForOpticalDuplicates = (byte) (id % 6);
        ends.score = (short) (30 * id);
        if (ends instanceof ReadEndsForMarkDuplicatesWithBarcodes) {
            final ReadEndsForMarkDuplicatesWithBarcodes withBarcodes = (ReadEndsForMarkDuplicatesWithBarcodes) ends;
            withBarcodes.barcode = 7 * id;
            withBarcodes.readOneBarcode = -id;
            withBarcodes.readTwoBarcode = id + 13;
        }
        return ends;
    }

    /** Runs each operation on the map and on a reference model and compares the results. */
    private static final class Checker {
        private final Kind kind;
        private final ReadEndsForMarkDuplicatesMap map;
        private final Map<Integer, Map<String, ReadEndsForMarkDuplicates>> model = new HashMap<>();
        /** The disk-based map holds the sequence last passed to remove() in RAM. */
        private Integer sequenceInRam = null;

        Checker(final Kind kind, final int maxOpenFiles) {
            this.kind = kind;
            this.map = newMap(kind, maxOpenFiles);
        }

        void put(final int mateSequence, final String key, final int id) {
            final ReadEndsForMarkDuplicates ends = ends(kind, id, mateSequence);
            map.put(mateSequence, key, ends);
            model.computeIfAbsent(mateSequence, s -> new HashMap<>()).put(key, ends);
            checkSizes();
        }

        void remove(final int mateSequence, final String key) {
            final ReadEndsForMarkDuplicates actual = map.remove(mateSequence, key);
            final ReadEndsForMarkDuplicates expected = model.getOrDefault(mateSequence, new HashMap<>()).remove(key);
            sequenceInRam = mateSequence;
            if (expected == null) {
                Assert.assertNull(actual, key);
            } else {
                Assert.assertNotNull(actual, key);
                ReadEndsForMarkDuplicatesTest.assertSameFields(actual, expected);
                if (expected instanceof ReadEndsForMarkDuplicatesWithBarcodes) {
                    final ReadEndsForMarkDuplicatesWithBarcodes a = (ReadEndsForMarkDuplicatesWithBarcodes) actual;
                    final ReadEndsForMarkDuplicatesWithBarcodes e = (ReadEndsForMarkDuplicatesWithBarcodes) expected;
                    Assert.assertEquals(a.barcode, e.barcode);
                    Assert.assertEquals(a.readOneBarcode, e.readOneBarcode);
                    Assert.assertEquals(a.readTwoBarcode, e.readTwoBarcode);
                }
            }
            checkSizes();
        }

        void checkSizes() {
            final int total = model.values().stream().mapToInt(Map::size).sum();
            Assert.assertEquals(map.size(), total);
            if (kind == Kind.MEMORY) {
                Assert.assertEquals(map.sizeInRam(), total);
            } else {
                final int inRam = sequenceInRam == null ? 0 : model.getOrDefault(sequenceInRam, new HashMap<>()).size();
                Assert.assertEquals(map.sizeInRam(), inRam);
            }
        }
    }

    @Test(dataProvider = "maps")
    public void testEmptyMap(final Kind kind, final int maxOpenFiles) {
        final Checker checker = new Checker(kind, maxOpenFiles);
        checker.checkSizes();
        checker.remove(0, "absent");
        checker.remove(5, "absent");
        checker.remove(0, "absent");
    }

    @Test(dataProvider = "maps")
    public void testPutAndRemoveAcrossManySequences(final Kind kind, final int maxOpenFiles) {
        final Checker checker = new Checker(kind, maxOpenFiles);

        // before any remove nothing is held in RAM by the disk-based map
        checker.put(2, "rg1:early", 1);

        // the first read on sequence 0 looks for its mate, then stores itself
        checker.remove(0, "rg1:a");
        checker.put(0, "rg1:a", 2);
        checker.put(1, "rg1:b", 3);
        checker.put(2, "rg1:c", 4);
        checker.put(3, "rg1:d", 5);
        checker.put(4, "rg1:e", 6);
        // sequences whose output files were closed to stay under the limit are appended to again
        checker.put(2, "rg2:c", 7);
        checker.put(1, "rg2:b", 8);
        checker.put(5, "rg1:f", 9);
        checker.put(0, "rg2:a", 10);

        checker.remove(0, "rg1:a");
        checker.remove(0, "rg1:missing");
        checker.remove(0, "rg1:a");

        checker.remove(1, "rg2:b");
        checker.put(1, "rg3:b", 11);
        checker.put(3, "rg2:d", 12);
        checker.remove(1, "rg1:b");

        checker.remove(2, "rg1:c");
        checker.remove(2, "rg1:early");
        checker.remove(2, "rg2:c");
        checker.remove(3, "rg2:d");
        checker.remove(4, "rg1:e");

        // back to sequences that were spilled after being loaded
        checker.remove(1, "rg3:b");
        checker.remove(0, "rg2:a");
        checker.remove(3, "rg1:d");
        checker.remove(5, "rg1:f");
        checker.remove(9, "rg1:never");
    }

    /** Several records per sequence file, more than one file per handle, then everything back. */
    @Test(dataProvider = "maps")
    public void testManyRecordsPerSequence(final Kind kind, final int maxOpenFiles) {
        final Checker checker = new Checker(kind, maxOpenFiles);
        checker.remove(0, "start");
        int id = 0;
        for (int round = 0; round < 4; round++) {
            for (int sequence = 0; sequence < 6; sequence++) {
                checker.put(sequence, "rg:" + sequence + ":" + round, id++);
            }
        }
        for (int sequence = 0; sequence < 6; sequence++) {
            for (int round = 3; round >= 0; round--) {
                checker.remove(sequence, "rg:" + sequence + ":" + round);
            }
        }
        Assert.assertEquals(checker.map.size(), 0);
    }

    @Test
    public void testMemoryMapRemoveBeyondKnownSequences() {
        final MemoryBasedReadEndsForMarkDuplicatesMap map = new MemoryBasedReadEndsForMarkDuplicatesMap();
        Assert.assertNull(map.remove(Integer.MAX_VALUE, "x"));
        final ReadEndsForMarkDuplicates ends = ends(Kind.MEMORY, 1, 3);
        map.put(3, "x", ends);
        Assert.assertNull(map.remove(1, "x"));
        Assert.assertNull(map.remove(4, "x"));
        Assert.assertEquals(map.size(), 1);
        ReadEndsForMarkDuplicatesTest.assertSameFields(map.remove(3, "x"), ends);
        Assert.assertEquals(map.size(), 0);
        Assert.assertEquals(map.sizeInRam(), 0);
    }
}
