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

package picard.sam.util;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

import java.util.function.Supplier;

/**
 * Tests for {@link PhysicalLocation}, {@link PhysicalLocationInt} and {@link PhysicalLocationShort}: the unset
 * defaults, {@link PhysicalLocation#hasLocation()}, the tile, x and y getters and setters, the range in which each
 * variant stores x and y, and the read group and library accessors that the base classes do not implement.
 */
public class PhysicalLocationTest {

    @DataProvider(name = "locationFactories")
    public Object[][] locationFactories() {
        return new Object[][]{
                {(Supplier<PhysicalLocation>) PhysicalLocationInt::new},
                {(Supplier<PhysicalLocation>) PhysicalLocationShort::new}
        };
    }

    @Test
    public void testNoValue() {
        Assert.assertEquals(PhysicalLocation.NO_VALUE, -1);
    }

    @Test(dataProvider = "locationFactories")
    public void testDefaultsAreUnset(final Supplier<PhysicalLocation> factory) {
        final PhysicalLocation loc = factory.get();
        Assert.assertEquals(loc.getTile(), PhysicalLocation.NO_VALUE);
        Assert.assertEquals(loc.getX(), PhysicalLocation.NO_VALUE);
        Assert.assertEquals(loc.getY(), PhysicalLocation.NO_VALUE);
        Assert.assertFalse(loc.hasLocation());
    }

    @DataProvider(name = "hasLocationByTile")
    public Object[][] hasLocationByTile() {
        final short[] tiles = {-1, 0, 1, -2, 1101, Short.MAX_VALUE, Short.MIN_VALUE};
        final Object[][] rows = new Object[tiles.length * 2][];
        int i = 0;
        for (final short tile : tiles) {
            // Only the tile decides; x and y are left unset.
            rows[i++] = new Object[]{new PhysicalLocationInt(), tile, tile != -1};
            rows[i++] = new Object[]{new PhysicalLocationShort(), tile, tile != -1};
        }
        return rows;
    }

    @Test(dataProvider = "hasLocationByTile")
    public void testHasLocationDependsOnlyOnTile(final PhysicalLocation loc, final short tile, final boolean expected) {
        loc.setTile(tile);
        Assert.assertEquals(loc.getTile(), tile);
        Assert.assertEquals(loc.getX(), -1);
        Assert.assertEquals(loc.getY(), -1);
        Assert.assertEquals(loc.hasLocation(), expected);
    }

    @DataProvider(name = "intCoordinates")
    public Object[][] intCoordinates() {
        return new Object[][]{
                {0}, {1}, {-1}, {Short.MAX_VALUE}, {Short.MAX_VALUE + 1}, {65535}, {65536}, {100000}, {-32769},
                {Integer.MAX_VALUE}, {Integer.MIN_VALUE}
        };
    }

    @Test(dataProvider = "intCoordinates")
    public void testIntVariantStoresFullIntRange(final int value) {
        final PhysicalLocationInt loc = new PhysicalLocationInt();
        loc.setX(value);
        Assert.assertEquals(loc.getX(), value);
        Assert.assertEquals(loc.x, value);
        Assert.assertEquals(loc.getY(), -1);

        loc.setY(value);
        Assert.assertEquals(loc.getY(), value);
        Assert.assertEquals(loc.y, value);

        loc.setX(7);
        Assert.assertEquals(loc.getX(), 7);
        Assert.assertEquals(loc.getY(), value);
    }

    /** x and y in the short variant keep only the low 16 bits, as a signed short. */
    @DataProvider(name = "shortCoordinates")
    public Object[][] shortCoordinates() {
        return new Object[][]{
                {0, 0},
                {1, 1},
                {-1, -1},
                {32767, 32767},
                {32768, -32768},
                {-32768, -32768},
                {-32769, 32767},
                {65535, -1},
                {65536, 0},
                {100000, -31072}, // 100000 = 0x186A0, low 16 bits 0x86A0
                {Integer.MAX_VALUE, -1},
                {Integer.MIN_VALUE, 0}
        };
    }

    @Test(dataProvider = "shortCoordinates")
    public void testShortVariantTruncatesToShort(final int value, final int expected) {
        final PhysicalLocationShort loc = new PhysicalLocationShort();
        loc.setX(value);
        Assert.assertEquals(loc.getX(), expected);
        Assert.assertEquals(loc.getY(), -1);

        loc.setY(value);
        Assert.assertEquals(loc.getY(), expected);
        Assert.assertEquals(loc.getX(), expected);
    }

    @DataProvider(name = "tiles")
    public Object[][] tiles() {
        return new Object[][]{{(short) 0}, {(short) 1}, {(short) 2228}, {Short.MAX_VALUE}, {Short.MIN_VALUE}};
    }

    @Test(dataProvider = "tiles")
    public void testTileRoundTrip(final short tile) {
        final PhysicalLocationInt intLoc = new PhysicalLocationInt();
        intLoc.setTile(tile);
        Assert.assertEquals(intLoc.getTile(), tile);
        Assert.assertEquals(intLoc.tile, tile);

        final PhysicalLocationShort shortLoc = new PhysicalLocationShort();
        shortLoc.setTile(tile);
        Assert.assertEquals(shortLoc.getTile(), tile);
    }

    @Test(dataProvider = "locationFactories")
    public void testReadGroupAndLibraryAreNotImplemented(final Supplier<PhysicalLocation> factory) {
        final PhysicalLocation loc = factory.get();
        assertNotImplemented(() -> loc.getReadGroup());
        assertNotImplemented(() -> loc.setReadGroup((short) 1));
        assertNotImplemented(() -> loc.getLibraryId());
        assertNotImplemented(() -> loc.setLibraryId((short) 1));
    }

    private static void assertNotImplemented(final Assert.ThrowingRunnable runnable) {
        final PicardException e = Assert.expectThrows(PicardException.class, runnable);
        Assert.assertTrue(e.getMessage().contains("Not Implemented"), e.getMessage());
    }
}
