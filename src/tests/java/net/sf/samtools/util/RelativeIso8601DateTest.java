package net.sf.samtools.util;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Date;
import java.util.List;

/** @author mccowan */

public class RelativeIso8601DateTest {
    @Test
    public void testLazyInstance() {
        final RelativeIso8601Date lazy = RelativeIso8601Date.generateLazyNowInstance();
        Assert.assertEquals(lazy.toString(), RelativeIso8601Date.LAZY_NOW_LABEL);
        Assert.assertEquals(lazy.toString(), RelativeIso8601Date.LAZY_NOW_LABEL);
        Assert.assertEquals(lazy.toString(), RelativeIso8601Date.LAZY_NOW_LABEL);
        Assert.assertEquals((double) lazy.getTime(), (double) System.currentTimeMillis(), 2000d); // Up to 2 seconds off; iso truncates milliseconds.
        // Assert no exception thrown; this should be valid, because toString should now return an iso-looking date.
        new RelativeIso8601Date(lazy.toString());
    }

    @Test
    public void testNonLazyInstance() {
        // Test both constructor methods
        final List<RelativeIso8601Date> testDates = Arrays.<RelativeIso8601Date>asList(
                new RelativeIso8601Date(new Date()),
                new RelativeIso8601Date(new Iso8601Date(new Date()).toString())
        );

        for (final RelativeIso8601Date nonLazy : testDates) {
            Assert.assertFalse(nonLazy.toString().equals(RelativeIso8601Date.LAZY_NOW_LABEL));
            Assert.assertEquals((double) nonLazy.getTime(), (double) System.currentTimeMillis(), 1d);
            // Assert no exception thrown; this should be valid, because toString return an iso-looking date.
            new RelativeIso8601Date(nonLazy.toString());
        }
    }
}
