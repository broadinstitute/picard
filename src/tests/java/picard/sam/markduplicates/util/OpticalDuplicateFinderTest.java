package picard.sam.markduplicates.util;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.Assert;
import picard.sam.util.PhysicalLocation;
import picard.sam.util.PhysicalLocationInt;
import picard.sam.util.PhysicalLocationShort;
import picard.sam.util.ReadNameParser;

import java.util.Arrays;

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

        final boolean[] opticalDuplicateFlags = opticalDuplicateFinder.findOpticalDuplicates(Arrays.asList(loc1, loc2));
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

        final boolean[] opticalDuplicateFlags = opticalDuplicateFinder.findOpticalDuplicates(Arrays.asList(loc1, loc2));
        for (final boolean opticalDuplicateFlag : opticalDuplicateFlags) {
            Assert.assertFalse(opticalDuplicateFlag);
        }
    }

}
