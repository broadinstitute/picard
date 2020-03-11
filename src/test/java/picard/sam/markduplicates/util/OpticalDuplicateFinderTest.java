package picard.sam.markduplicates.util;

import htsjdk.samtools.util.Log;
import org.testng.annotations.Test;
import org.testng.Assert;
import picard.sam.util.PhysicalLocation;
import picard.sam.util.PhysicalLocationInt;
import picard.sam.util.ReadNameParser;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 * Tests for OpticalDuplicateFinder
 *
 * @author Nils Homer
 */
public class OpticalDuplicateFinderTest {
    @Test
    public void testDefaultRegex() {
        final String readName1 = "000000000-ZZZZZ:1:1105:17981:23325";
        final String readName2 = "000000000-ZZZZZ:1:1109:22981:17995";

        final int[] tokens = new int[3];
        Assert.assertEquals(ReadNameParser.getLastThreeFields(readName1, ':', tokens), 5);
        Assert.assertEquals(ReadNameParser.getLastThreeFields(readName2, ':', tokens), 5);

        final OpticalDuplicateFinder opticalDuplicateFinder = new OpticalDuplicateFinder();
        final PhysicalLocation loc1 = new ReadEndsForMarkDuplicates();
        final PhysicalLocation loc2 = new ReadEndsForMarkDuplicates();

        Assert.assertTrue(opticalDuplicateFinder.addLocationInformation(readName1, loc1));
        Assert.assertTrue(opticalDuplicateFinder.addLocationInformation(readName2, loc2));

        final boolean[] opticalDuplicateFlags = opticalDuplicateFinder.findOpticalDuplicates(Arrays.asList(loc1, loc2), null);
        for (final boolean opticalDuplicateFlag : opticalDuplicateFlags) {
            Assert.assertFalse(opticalDuplicateFlag);
        }
    }

    @Test
    public void testVeryLongReadNames() {
        final String readName1 = "M01234:123:000000000-ZZZZZ:1:1105:17981:23325";
        final String readName2 = "M01234:123:000000000-ZZZZZ:1:1109:22981:17995";

        final int[] tokens = new int[3];
        Assert.assertEquals(ReadNameParser.getLastThreeFields(readName1, ':', tokens), 7);
        Assert.assertEquals(ReadNameParser.getLastThreeFields(readName2, ':', tokens), 7);

        final OpticalDuplicateFinder opticalDuplicateFinder = new OpticalDuplicateFinder();
        final PhysicalLocation loc1 = new ReadEndsForMarkDuplicates();
        final PhysicalLocation loc2 = new ReadEndsForMarkDuplicates();

        Assert.assertTrue(opticalDuplicateFinder.addLocationInformation(readName1, loc1));
        Assert.assertTrue(opticalDuplicateFinder.addLocationInformation(readName2, loc2));

        final boolean[] opticalDuplicateFlags = opticalDuplicateFinder.findOpticalDuplicates(Arrays.asList(loc1, loc2), null);
        for (final boolean opticalDuplicateFlag : opticalDuplicateFlags) {
            Assert.assertFalse(opticalDuplicateFlag);
        }
    }

    @Test
    public void testKeeper() {
        final Log log = Log.getInstance(OpticalDuplicateFinderTest.class);
        final OpticalDuplicateFinder finder = new OpticalDuplicateFinder(OpticalDuplicateFinder.DEFAULT_READ_NAME_REGEX, 100, log);
        final List<PhysicalLocation> locs = Arrays.asList(
                loc(7, 1500, 1500),
                loc(7, 1501, 1501),
                loc(5, 1500, 1500),
                loc(7, 1490, 1502),
                loc(7, 2500, 2500),
                loc(7,   10,   10)
        );

        assertEquals(finder.findOpticalDuplicates(locs, null       ), new boolean[] {true, true,  false, false, false, false});
        assertEquals(finder.findOpticalDuplicates(locs, locs.get(0)), new boolean[] {false, true,  false, true, false, false});
        assertEquals(finder.findOpticalDuplicates(locs, locs.get(1)), new boolean[] {true,  false, false, true, false, false});
        assertEquals(finder.findOpticalDuplicates(locs, locs.get(3)), new boolean[] {true,  true,  false, false, false, false});

        for (int i=0; i<100; ++i) {
            final Random random = new Random(i);
            final List<PhysicalLocation> shuffled = new ArrayList<>(locs);
            final List<PhysicalLocation> keepers  = Arrays.asList(locs.get(0), locs.get(1), locs.get(3));
            final PhysicalLocation keeper = keepers.get(random.nextInt(keepers.size()));
            Collections.shuffle(shuffled);

            int opticalDupeCount = countTrue(finder.findOpticalDuplicates(shuffled, keeper));
            Assert.assertEquals(opticalDupeCount, 2);
        }
    }

    /**
     * Tests the case where the "keeper" record is not in the list that is passed to the OpticalDuplicateFinder. This can happen
     * when there are, e.g. FR and RF reads, which can all be molecular duplicates of one another, but cannot be duplicates of one
     * another and are thus partitioned into two sets for optical duplicate checking.
     */
    @Test
    public void testKeeperNotInList() {
        final Log log = Log.getInstance(OpticalDuplicateFinderTest.class);
        final OpticalDuplicateFinder finder = new OpticalDuplicateFinder(OpticalDuplicateFinder.DEFAULT_READ_NAME_REGEX, 100, log);
        final List<PhysicalLocation> locs = Arrays.asList(
                loc(1, 100, 100),
                loc(1, 101, 101),
                loc(1,  99, 99),
                loc(1,  99, 102)
        );

        Assert.assertEquals(countTrue(finder.findOpticalDuplicates(locs, loc(7, 5000, 5000))), 3);
    }

    @Test
    public void testKeeperAtEndWithinCliqueOfAllOpticalDuplicates() {
        final Log log = Log.getInstance(OpticalDuplicateFinderTest.class);
        final OpticalDuplicateFinder finder = new OpticalDuplicateFinder(OpticalDuplicateFinder.DEFAULT_READ_NAME_REGEX, 15, log);
        final List<PhysicalLocation> locs = Arrays.asList(
                loc(1, 10, 0),
                loc(1, 20, 0),
                loc(1, 30, 0)
        );

        assertEquals(finder.findOpticalDuplicates(locs, locs.get(2)), new boolean[] {true, true, false});
    }

    @Test
    public void testKeeperOrderDependence() {
        final Log log = Log.getInstance(OpticalDuplicateFinderTest.class);
        final OpticalDuplicateFinder finder = new OpticalDuplicateFinder(OpticalDuplicateFinder.DEFAULT_READ_NAME_REGEX, 100, log);
        final List<PhysicalLocation> locs = Arrays.asList(
                loc(1, 100, 190),
                loc(1, 100, 280),
                loc(1,  100, 370),
                loc(1,  100, 460)
        );
        final List<PhysicalLocation> locsReordered = Arrays.asList(
                loc(1, 100, 190),
                loc(1, 100, 460),
                loc(1,  100, 370),
                loc(1,  100, 280)
        );

        Assert.assertEquals(countTrue(finder.findOpticalDuplicates(locs, locs.get(0))), 3);
        Assert.assertEquals(countTrue(finder.findOpticalDuplicates(locsReordered, locsReordered.get(0))), 3);
    }

    @Test
    public void testKeeperOrderDependenceTwoDistinctClusters() {
        final OpticalDuplicateFinder finder = new OpticalDuplicateFinder(OpticalDuplicateFinder.DEFAULT_READ_NAME_REGEX, 100, null);
        final List<PhysicalLocation> locs = Arrays.asList(
                loc(1, 100, 190),
                loc(1, 100, 280),
                // A second cluster of reads that should result in 3 extra optical duplicates
                loc(1, 300, 500),
                loc(1, 300, 590),
                loc(1, 300, 670),
                loc(1, 300, 760)
        );
        // Reordered the second cluster from the keeper cluster so they will not be grouped together
        final List<PhysicalLocation> locsReordered = Arrays.asList(
                loc(1, 100, 190),
                loc(1, 100, 280),
                loc(1, 300, 500),
                loc(1, 300, 760),
                loc(1, 300, 670),
                loc(1, 300, 590)
        );

        Assert.assertEquals(countTrue(finder.findOpticalDuplicates(locs, locs.get(0))), 4);
        Assert.assertEquals(countTrue(finder.findOpticalDuplicates(locsReordered, locsReordered.get(0))), 4);
    }

    @Test
    public void testKeeperOrderIndependenceTwoDistinctClusters() {
        final OpticalDuplicateFinder finder = new OpticalDuplicateFinder(OpticalDuplicateFinder.DEFAULT_READ_NAME_REGEX, 100, null);
        final List<PhysicalLocation> locs = Arrays.asList(
                loc(1, 100, 190),
                loc(1, 100, 280),
                // A second cluster of reads that should result in 3 extra optical duplicates
                loc(1, 300, 500),
                loc(1, 300, 590),
                loc(1, 300, 670),
                loc(1, 300, 760)
        );
        // Reordered the second cluster from the keeper cluster so they will not be grouped together
        final List<PhysicalLocation> locsReordered = Arrays.asList(
                loc(1, 100, 190),
                loc(1, 100, 280),
                loc(1, 300, 590),
                loc(1, 300, 760),
                loc(1, 300, 670),
                loc(1, 300, 500)
        );

        // In both cases we expect "loc(1, 300, 500)" to be the non-optical-duplicate
        Assert.assertEquals(countTrue(finder.findOpticalDuplicates(locs, locs.get(0))), 4);
        Assert.assertFalse(finder.findOpticalDuplicates(locs, locs.get(0))[2]);
        Assert.assertEquals(countTrue(finder.findOpticalDuplicates(locsReordered, locsReordered.get(0))), 4);
        Assert.assertFalse(finder.findOpticalDuplicates(locsReordered, locs.get(0))[5]);
    }

    @Test
    public void testMaxSetSize() {
        final List<PhysicalLocation> locs = Arrays.asList(
                loc(7, 1500, 1500),
                loc(7, 1501, 1501),
                loc(7, 1490, 1502));

        // there should be at least 1 optical duplicate
        final OpticalDuplicateFinder normalFinder = new OpticalDuplicateFinder(OpticalDuplicateFinder.DEFAULT_READ_NAME_REGEX, 100, null);
        Assert.assertTrue(countTrue(normalFinder.findOpticalDuplicates(locs, null)) > 0);

        // there should be zero optical duplicates
        final OpticalDuplicateFinder constrainedFinder = new OpticalDuplicateFinder(OpticalDuplicateFinder.DEFAULT_READ_NAME_REGEX, 100, 1, null);
        Assert.assertEquals(countTrue(constrainedFinder.findOpticalDuplicates(locs, null)), 0);
    }

    /** Helper method to create a physical location. */
    private PhysicalLocation loc(final int tile, final int x, final int y) {
        final PhysicalLocation l = new PhysicalLocationInt() {
            @Override
            public short getReadGroup() { return 1; }
        };
        l.setTile((short) tile);
        l.setX(x);
        l.setY(y);
        return l;
    }

    void assertEquals(final boolean[] actual, final boolean[] expected) {
        if (!Arrays.equals(actual, expected)) {
            throw new AssertionError("expected: " + Arrays.toString(expected) + " but was: " + Arrays.toString(actual));
        }
    }

    /** Simply counts the true values in a boolean array. */
    int countTrue(final boolean[] bs) {
        int count = 0;
        for (final boolean b : bs) if (b) ++count;
        return count;
    }
}
