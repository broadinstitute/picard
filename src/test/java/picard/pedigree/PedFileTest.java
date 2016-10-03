package picard.pedigree;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by farjoun on 9/19/14.
 */
public class PedFileTest {

    @DataProvider()
    public Object[][] testFromSexMapDataProvider() {
        return new Object[][]{
                new Object[]{Arrays.asList(), Arrays.asList()},
                new Object[]{Arrays.asList("female1", "female2", "female3"), Arrays.asList("male1", "male2", "male3")},
                new Object[]{Arrays.asList("female1", "female2", "female3"), Arrays.asList()},
                new Object[]{Arrays.asList("female1"), Arrays.asList("male1", "male2", "male3")},
        };
    }

    @Test(dataProvider = "testFromSexMapDataProvider")
    public void testFromSexMap(final Collection<String> females, final Collection<String> males) throws Exception {
        final Map<String, Sex> data = new HashMap<String, Sex>();
        for (final String sample : females) {
            data.put(sample, Sex.Female);
        }
        for (final String sample : males) {
            data.put(sample, Sex.Male);
        }

        final PedFile pedFile = PedFile.fromSexMap(data);

        // Check that sizes agree
        Assert.assertEquals(pedFile.size(), females.size() + males.size());

        // Check that every entry in the PedFile came from one of the two collections
        for (final Map.Entry<String, PedFile.PedTrio> pedTrioEntry : pedFile.entrySet()) {
            Assert.assertTrue(females.contains(pedTrioEntry.getValue().getIndividualId()) |
                    males.contains(pedTrioEntry.getValue().getIndividualId()));
        }
        // Check that all the females are there and are listed as female
        for (final String female : females) {
            Assert.assertEquals(pedFile.get(female).getSex(), Sex.Female);
        }

        // Check that all the males are there and are listed as male
        for (final String male : males) {
            Assert.assertEquals(pedFile.get(male).getSex(), Sex.Male);
        }

    }
}
