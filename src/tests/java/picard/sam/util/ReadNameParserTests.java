package picard.sam.util;

import htsjdk.samtools.util.CollectionUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;

/**
 * Tests for the ReadNameParser class.
 */
public class ReadNameParserTests {
    /** Tests rapidParseInt for positive and negative numbers, as well as non-digit suffixes */
    @Test
    public void testRapidParseInt() {
        for (int i = -100; i < 100; i++) {
            Assert.assertEquals(ReadNameParser.rapidParseInt(Integer.toString(i)), i);

            // trailing characters
            Assert.assertEquals(ReadNameParser.rapidParseInt(Integer.toString(i)+"A"), i);
            Assert.assertEquals(ReadNameParser.rapidParseInt(Integer.toString(i)+"ACGT"), i);
            Assert.assertEquals(ReadNameParser.rapidParseInt(Integer.toString(i)+".1"), i);
        }
    }

    /** Tests rapidParseInt for positive and negative numbers, as well as non-digit suffixes */
    @Test
    public void testRapidParseIntFails() {
        List<String> values = CollectionUtil.makeList("foo", "bar", "abc123", "-foo", "f00", "-f00");
        for (String s : values) {
            try {
                ReadNameParser.rapidParseInt(s);
                Assert.fail("Should have failed to rapid-parse " + s + " as an int.");
            }
            catch (NumberFormatException nfe) {
                /* expected */
            }
        }
    }

    /** Helper for testGetRapidDefaultReadNameRegexSplit */
    private void doTestGetRapidDefaultReadNameRegexSplit(int numFields) {
        final int[] inputFields = new int[3];
        final int[] expectedFields = new int[3];
        String readName = "";
        for (int i = 0; i < numFields; i++) {
            if (0 < i) readName += ":";
            readName += Integer.toString(i);
        }
        inputFields[0] = inputFields[1] = inputFields[2] = -1;
        if (numFields < 3) {
            Assert.assertEquals(ReadNameParser.getLastThreeFields(readName, ':', inputFields), -1);
        }
        else {
            Assert.assertEquals(ReadNameParser.getLastThreeFields(readName, ':', inputFields), numFields);
            expectedFields[0] = expectedFields[1] = expectedFields[2] = -1;
            if (0 < numFields) expectedFields[0] = numFields-3;
            if (1 < numFields) expectedFields[1] = numFields-2;
            if (2 < numFields) expectedFields[2] = numFields-1;
            for (int i = 0; i < inputFields.length; i++) {
                Assert.assertEquals(inputFields[i], expectedFields[i]);
            }
        }
    }

    /** Tests that we split the string with the correct # of fields, and modified values */
    @Test
    public void testGetRapidDefaultReadNameRegexSplit() {
        for (int i = 1; i < 10; i++) {
            doTestGetRapidDefaultReadNameRegexSplit(i);
        }
    }

    @DataProvider(name = "testParseReadNameDataProvider")
    public Object[][] testParseReadNameDataProvider() {
        return new Object[][]{
                {"RUNID:7:1203:2886:82292", 1203, 2886, 82292},
                {"RUNID:7:1203:2884:16834", 1203, 2884, 16834}
        };
    }

    // NB: these test fail s due to overflow in the duplicate finder test.  This has been the behavior previously, so keep it for now.
    @Test(dataProvider = "testParseReadNameDataProvider", enabled = true)
    public void testParseReadNameOverflow(final String readName, final int tile, final int x, final int y) {
        ReadNameParser parser = new ReadNameParser();
        PhysicalLocation loc = new PhysicalLocationShort();
        Assert.assertTrue(parser.addLocationInformation(readName, loc));
        Assert.assertEquals(loc.getTile(), tile);
        Assert.assertEquals(loc.getX(), (short)x); // casting to short for the overflow
        Assert.assertEquals(loc.getY(), (short)y); // casting to short for the overflow
    }

    // NB: this test the case where we do not overflow in the duplicate finder test.
    @Test(dataProvider = "testParseReadNameDataProvider", enabled = true)
    public void testParseReadNameOK(final String readName, final int tile, final int x, final int y) {
        ReadNameParser parser = new ReadNameParser();
        PhysicalLocation loc = new PhysicalLocationInt();
        Assert.assertTrue(parser.addLocationInformation(readName, loc));
        Assert.assertEquals(loc.getTile(), tile);
        Assert.assertEquals(loc.getX(), x); // we store ints, so we should not overflow
        Assert.assertEquals(loc.getY(), y); // we store ints, so we should not overflow
    }

    @DataProvider(name = "testReadNameParsing")
    public Object[][] testReadNameParsingDataProvider() {
        final String lastThreeFieldsRegex = "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$";
        return new Object[][]{
                {lastThreeFieldsRegex, "RUNID:123:000000000-ZZZZZ:1:1105:17981:23325", 1105, 17981, 23325, true},
                {lastThreeFieldsRegex, "RUNID:123:000000000-ZZZZZ:1:1109:22981:17995", 1109, 22981, 17995, true},
                {lastThreeFieldsRegex, "1109:22981:17995", 1109, 22981, 17995, true},
                {lastThreeFieldsRegex, "RUNID:7:1203:2886:82292", 1203, 2886, 82292, true},
                {lastThreeFieldsRegex, "RUNID:7:1203:2884:16834", 1203, 2884, 16834, true},
                {lastThreeFieldsRegex, "1109ABC:22981DEF:17995GHI", 1109, 22981, 17995, true},
                {ReadNameParser.DEFAULT_READ_NAME_REGEX, "RUNID:123:000000000-ZZZZZ:1:1105:17981:23325", 1105, 17981, 23325, true},
                {ReadNameParser.DEFAULT_READ_NAME_REGEX, "RUNID:123:000000000-ZZZZZ:1:1109:22981:17995", 1109, 22981, 17995, true},
                {ReadNameParser.DEFAULT_READ_NAME_REGEX, "1109:22981:17995", 1109, 22981, 17995, false},
                {ReadNameParser.DEFAULT_READ_NAME_REGEX, "RUNID:7:1203:2886:82292", 1203, 2886, 82292, true},
                {ReadNameParser.DEFAULT_READ_NAME_REGEX, "RUNID:7:1203:2884:16834", 1203, 2884, 16834, true}
        };
    }

    @Test(dataProvider = "testReadNameParsing")
    public void testReadNameParsing(final String readNameRegex, final String readName, final int tile, final int x, final int y, final boolean addLocationInformationSucceeds) {
        final ReadNameParser parser = new ReadNameParser(readNameRegex);
        final PhysicalLocationInt loc = new PhysicalLocationInt();
        Assert.assertEquals(parser.addLocationInformation(readName, loc), addLocationInformationSucceeds);
        if (addLocationInformationSucceeds) { // just check the location
            Assert.assertEquals(loc.getTile(), tile);
            Assert.assertEquals(loc.getX(), x);
            Assert.assertEquals(loc.getY(), y);
        }
        else if (readNameRegex == ReadNameParser.DEFAULT_READ_NAME_REGEX) { // additional testing on the default regex
            int[] tokens = new int[3];
            ReadNameParser.getLastThreeFields(readName, ':', tokens);
            Assert.assertEquals(tokens[0], tile);
            Assert.assertEquals(tokens[1], x);
            Assert.assertEquals(tokens[2], y);
        }
    }

}
