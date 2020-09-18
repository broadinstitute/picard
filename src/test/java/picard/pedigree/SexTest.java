package picard.pedigree;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

import java.util.Arrays;

public class SexTest {
    @DataProvider(name = "testSexDataProviderInt")
    public Object[][] testSexDataProviderInt() {
        return new Object[][]{
                {Sex.Male, 1},
                {Sex.Female, 2},
                {Sex.Unknown, 0},
        };
    }

    @DataProvider(name = "testSexDataProviderStr")
    public Object[][] testSexDataProviderStr() {
        return new Object[][]{
                {Sex.Female, "F"},
                {Sex.Female, "f"},
                {Sex.Female, "Female"},
                {Sex.Female, "female"},
                {Sex.Male, "M"},
                {Sex.Male, "m"},
                {Sex.Male, "Male"},
                {Sex.Male, "male"},
                {Sex.Unknown, "Unknown"},
                {Sex.Unknown, "U"},
                {Sex.Unknown, "unknown"},
                {Sex.Unknown, "u"},
        };
    }

    /**
     * Test Methods within Sex enum
     */

    @Test(dataProvider = "testSexDataProviderInt")
    public void testFromCode(final Sex sex, final int i) {
        Assert.assertTrue(sex == Sex.fromCode(i));
    }

    @Test(dataProvider = "testSexDataProviderStr")
    public void testFromString(final Sex sex, final String s) {
        Assert.assertTrue(sex == Sex.fromString(s));
    }

    @Test(expectedExceptions = PicardException.class)
    public void testFromStringException() {
        Sex.fromString("notASex");
    }

}
