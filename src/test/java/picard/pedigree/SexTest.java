package picard.pedigree;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

import java.util.Arrays;

public class SexTest {
    @DataProvider(name="testSexDataProviderInt")
    public Object[][] testSexDataProviderInt() {
        return new Object[][] {
                new Object[]{Sex.Male, 1, true},
                new Object[]{Sex.Female, 2, true},
                new Object[]{Sex.Unknown, 0, true},
                new Object[]{Sex.NotReported, 0, false},
                new Object[]{Sex.Female, 1, false},
                new Object[]{Sex.Unknown, 2, false},
                };
    }

    @DataProvider(name="testSexDataProviderStr")
    public Object[][] testSexDataProviderStr() {
        return new Object[][] {
                new Object[]{Sex.Male, "M", true},
                new Object[]{Sex.Female, "Female", true},
                new Object[]{Sex.Unknown, "Unknown", true},
                new Object[]{Sex.NotReported, "NotReported", true},
                new Object[]{Sex.Female, "M", false},
                new Object[]{Sex.Unknown, "female", false}
        };
    }

    /**
     *Test Methods within Sex enum
     */

    @Test(dataProvider="testSexDataProviderInt")
    public void testFromCode(final Sex sex, final int i, final boolean b){
        Assert.assertEquals(sex == Sex.fromCode(i), b);
    }

    @Test(dataProvider="testSexDataProviderStr")
    public void testFromString(final Sex sex, final String s, final boolean b){
        Assert.assertEquals(sex == Sex.fromString(s), b);
    }

    @Test(expectedExceptions = PicardException.class)
    public void testFromStringException(){
        Sex.fromString("notASex");
    }

}
