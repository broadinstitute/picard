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

import htsjdk.samtools.util.Log;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;

/**
 * Tests {@link ReadNameParser}: the default colon-split parsing of five- and seven-field names, custom regular
 * expressions with three and five capture groups, names that do not parse and the single warning per parser, a null
 * regex, coordinates at the limits of an int, the cache of the last parsed name, and the static helpers
 * {@link ReadNameParser#getLastThreeFields} and {@link ReadNameParser#rapidParseInt}.
 */
public class ReadNameParserBehaviorTest {

    /** The example regex from the READ_NAME_REGEX documentation. */
    private static final String LAST_THREE_FIELDS_REGEX = "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$";

    private static final String SEMICOLON_REGEX = "(?:.*;)?([0-9]+)[^;]*;([0-9]+)[^;]*;([0-9]+)[^;]*$";

    private static void assertLocation(final PhysicalLocation loc, final int tile, final int x, final int y) {
        Assert.assertEquals(loc.getTile(), tile, "tile");
        Assert.assertEquals(loc.getX(), x, "x");
        Assert.assertEquals(loc.getY(), y, "y");
    }

    private static void assertUntouched(final PhysicalLocation loc) {
        assertLocation(loc, -1, -1, -1);
        Assert.assertFalse(loc.hasLocation());
    }

    /** Runs the action with the global log stream redirected and returns what was written. */
    private static String captureLog(final Runnable action) {
        final PrintStream originalStream = Log.getGlobalPrintStream();
        final Log.LogLevel originalLevel = Log.getGlobalLogLevel();
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try (final PrintStream capture = new PrintStream(bytes, true, StandardCharsets.UTF_8)) {
            Log.setGlobalPrintStream(capture);
            Log.setGlobalLogLevel(Log.LogLevel.INFO);
            action.run();
        } finally {
            Log.setGlobalPrintStream(originalStream);
            Log.setGlobalLogLevel(originalLevel);
        }
        return bytes.toString(StandardCharsets.UTF_8);
    }

    private static int occurrences(final String text, final String fragment) {
        int count = 0;
        for (int i = text.indexOf(fragment); i >= 0; i = text.indexOf(fragment, i + 1)) {
            count++;
        }
        return count;
    }

    @DataProvider(name = "defaultRegexNames")
    public Object[][] defaultRegexNames() {
        return new Object[][]{
                // read name, parsed, tile, x, y
                {"RUN:1:1101:1000:2000", true, 1101, 1000, 2000},
                {"M01234:123:000000000-ABCDE:1:1101:1000:2000", true, 1101, 1000, 2000},
                {"M01234:123:FLOWCELL:LANE:1101:1000:2000", true, 1101, 1000, 2000},
                {"RUN:1:1101:1000:2000/1", true, 1101, 1000, 2000},
                {"RUN:1:1101X:1000Y:2000Z", true, 1101, 1000, 2000},
                {"RUN:1:1101:-1000:-2000", true, 1101, -1000, -2000},
                {"RUN:1:0:0:0", true, 0, 0, 0},
                {"RUN:1:32767:5:6", true, 32767, 5, 6},
                {"RUN:1:1101:2147483647:-2147483648", true, 1101, Integer.MAX_VALUE, Integer.MIN_VALUE},
                {"A::1101:1000:2000", true, 1101, 1000, 2000},
                {"1101:1000:2000", false, 0, 0, 0},
                {"A:1101:1000:2000", false, 0, 0, 0},
                {"A:B:C:1101:1000:2000", false, 0, 0, 0},
                {"A:B:C:D:E:1101:1000:2000", false, 0, 0, 0},
                {"1000:2000", false, 0, 0, 0},
                {"read1", false, 0, 0, 0},
                {"", false, 0, 0, 0},
                {"RUN:1:X1101:1000:2000", false, 0, 0, 0},
                {"RUN:1:1101::2000", false, 0, 0, 0},
                {"RUN:1:1101:1000:", false, 0, 0, 0},
                {"RUN:1:1101:-:2000", false, 0, 0, 0},
                {"RUN:1:1101:+5:2000", false, 0, 0, 0},
        };
    }

    @Test(dataProvider = "defaultRegexNames")
    public void testDefaultRegex(final String readName, final boolean parsed, final int tile, final int x, final int y) {
        for (final ReadNameParser parser : new ReadNameParser[]{
                new ReadNameParser(),
                new ReadNameParser(ReadNameParser.DEFAULT_READ_NAME_REGEX),
                new ReadNameParser(ReadNameParser.DEFAULT_READ_NAME_REGEX, Log.getInstance(ReadNameParserBehaviorTest.class))}) {
            final PhysicalLocationInt loc = new PhysicalLocationInt();
            Assert.assertEquals(parser.addLocationInformation(readName, loc), parsed, readName);
            if (parsed) {
                assertLocation(loc, tile, x, y);
            } else {
                assertUntouched(loc);
            }
        }
    }

    @DataProvider(name = "customRegexNames")
    public Object[][] customRegexNames() {
        return new Object[][]{
                // regex, read name, parsed, tile, x, y
                {LAST_THREE_FIELDS_REGEX, "1101:1000:2000", true, 1101, 1000, 2000},
                {LAST_THREE_FIELDS_REGEX, "A:1101:1000:2000", true, 1101, 1000, 2000},
                {LAST_THREE_FIELDS_REGEX, "RUN:1:1101:1000:2000", true, 1101, 1000, 2000},
                {LAST_THREE_FIELDS_REGEX, "A:B:C:1101:1000:2000", true, 1101, 1000, 2000},
                {LAST_THREE_FIELDS_REGEX, "M01234:123:000000000-ABCDE:1:1101:1000:2000", true, 1101, 1000, 2000},
                {LAST_THREE_FIELDS_REGEX, "RUN:1:1101:1000:2000/1", true, 1101, 1000, 2000},
                {LAST_THREE_FIELDS_REGEX, "RUN:1:007:0010:0", true, 7, 10, 0},
                {LAST_THREE_FIELDS_REGEX, "RUN:1:1101:2147483647:2147483647", true, 1101, Integer.MAX_VALUE, Integer.MAX_VALUE},
                {LAST_THREE_FIELDS_REGEX, "RUN:1:1101:-1000:2000", false, 0, 0, 0},
                {LAST_THREE_FIELDS_REGEX, "RUN:1:X:1000:2000", false, 0, 0, 0},
                {LAST_THREE_FIELDS_REGEX, "1000:2000", false, 0, 0, 0},
                {LAST_THREE_FIELDS_REGEX, "RUN:1:99999999999:1000:2000", false, 0, 0, 0},
                {SEMICOLON_REGEX, "RUN;1;1101;1000;2000", true, 1101, 1000, 2000},
                {SEMICOLON_REGEX, "RUN:1:1101:1000:2000", false, 0, 0, 0},
                // matches() requires the whole name to match
                {"([0-9]+):([0-9]+):([0-9]+)", "1101:1000:2000", true, 1101, 1000, 2000},
                {"([0-9]+):([0-9]+):([0-9]+)", "X:1101:1000:2000", false, 0, 0, 0},
                {"([0-9]+):([0-9]+):([0-9]+)", "1101:1000:2000:", false, 0, 0, 0},
                // groups are taken in order: tile, x, y
                {"x([0-9]+)_y([0-9]+)_t([0-9]+)", "x10_y20_t30", true, 10, 20, 30},
                {"[A-Z]+-(\\d+)-(\\d+)-(\\d+)", "READ-3-200-100", true, 3, 200, 100},
                {"(?:[^:]*:){2}(-?[0-9]+):(-?[0-9]+):(-?[0-9]+)", "RUN:1:12:-5:-6", true, 12, -5, -6},
                // five capture groups: only the first three are read
                {"([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)", "11_22_33_44_55", true, 11, 22, 33},
                {"([0-9]+)_([0-9]+)_([0-9]+)_([A-Z]+)_([A-Z]+)", "11_22_33_AB_CD", true, 11, 22, 33},
                {"([A-Z]+):([0-9]+):([0-9]+):([0-9]+):([0-9]+)", "LANE:1101:1000:2000:7", false, 0, 0, 0},
                {"([0-9]*):([0-9]+):([0-9]+)", ":1000:2000", false, 0, 0, 0},
        };
    }

    @Test(dataProvider = "customRegexNames")
    public void testCustomRegex(final String regex, final String readName, final boolean parsed, final int tile, final int x, final int y) {
        for (final ReadNameParser parser : new ReadNameParser[]{
                new ReadNameParser(regex),
                new ReadNameParser(regex, Log.getInstance(ReadNameParserBehaviorTest.class))}) {
            final PhysicalLocationInt loc = new PhysicalLocationInt();
            Assert.assertEquals(parser.addLocationInformation(readName, loc), parsed, readName);
            if (parsed) {
                assertLocation(loc, tile, x, y);
            }
        }
    }

    /** A name that fails to parse leaves the location as it was, for failures detected before any field is stored. */
    @Test
    public void testFailedParseLeavesLocationUntouched() {
        final String[][] cases = {
                {LAST_THREE_FIELDS_REGEX, "no_coordinates_here"},
                {LAST_THREE_FIELDS_REGEX, "RUN:1:99999999999:1000:2000"},
                {"([A-Z]+):([0-9]+):([0-9]+)", "ABC:1000:2000"},
                {ReadNameParser.DEFAULT_READ_NAME_REGEX, "RUN:1:1101:99999:"},
                {ReadNameParser.DEFAULT_READ_NAME_REGEX, "RUN:1:1101:1000"},
        };
        for (final String[] c : cases) {
            final PhysicalLocationInt loc = new PhysicalLocationInt();
            loc.setTile((short) 7);
            loc.setX(8);
            loc.setY(9);
            Assert.assertFalse(new ReadNameParser(c[0]).addLocationInformation(c[1], loc), c[1]);
            assertLocation(loc, 7, 8, 9);
        }
    }

    @Test
    public void testOverflowingCoordinatesAreRejectedByRegex() {
        final ReadNameParser parser = new ReadNameParser(LAST_THREE_FIELDS_REGEX);
        Assert.assertFalse(parser.addLocationInformation("RUN:1:1101:2147483648:2000", new PhysicalLocationInt()));
        Assert.assertFalse(parser.addLocationInformation("RUN:1:1101:1000:2147483648", new PhysicalLocationInt()));
        Assert.assertFalse(parser.addLocationInformation("RUN:1:2147483648:1000:2000", new PhysicalLocationInt()));
        final PhysicalLocationInt loc = new PhysicalLocationInt();
        Assert.assertTrue(parser.addLocationInformation("RUN:1:1101:2147483647:2000", loc));
        assertLocation(loc, 1101, Integer.MAX_VALUE, 2000);
    }

    /** The documented example regex and the default parsing agree on numeric five- and seven-field names. */
    @Test
    public void testDefaultParsingAgreesWithDocumentedRegex() {
        final String[] names = {
                "RUN:1:1101:1000:2000",
                "RUN:1:2228:0:0",
                "RUN:1:11101:25000:31000",
                "M01234:123:000000000-ABCDE:1:1105:17981:23325",
                "A00123:8:HABCDEFXX:4:2678:32741:37138",
                "RUN:1:1101:1000:2000/2",
                "RUN:1:1101:1000:2000_ACGT",
                "RUN:7:1203:2886:82292",
        };
        final ReadNameParser fast = new ReadNameParser();
        final ReadNameParser regex = new ReadNameParser(LAST_THREE_FIELDS_REGEX);
        for (final String name : names) {
            final PhysicalLocationInt fastLoc = new PhysicalLocationInt();
            final PhysicalLocationInt regexLoc = new PhysicalLocationInt();
            Assert.assertTrue(fast.addLocationInformation(name, fastLoc), name);
            Assert.assertTrue(regex.addLocationInformation(name, regexLoc), name);
            assertLocation(regexLoc, fastLoc.getTile(), fastLoc.getX(), fastLoc.getY());

            final String[] fields = name.split(":");
            final int n = fields.length;
            Assert.assertEquals(fastLoc.getTile(), Integer.parseInt(fields[n - 3]));
            Assert.assertEquals(fastLoc.getX(), Integer.parseInt(fields[n - 2]));
            Assert.assertEquals(fastLoc.getY(), Integer.parseInt(fields[n - 1].replaceAll("[^0-9].*$", "")));
        }
    }

    @Test
    public void testNullRegexNeverParses() {
        for (final ReadNameParser parser : new ReadNameParser[]{
                new ReadNameParser(null),
                new ReadNameParser(null, Log.getInstance(ReadNameParserBehaviorTest.class))}) {
            final String log = captureLog(() -> {
                final PhysicalLocationInt loc = new PhysicalLocationInt();
                Assert.assertFalse(parser.addLocationInformation("RUN:1:1101:1000:2000", loc));
                assertUntouched(loc);
                Assert.assertFalse(parser.addLocationInformation("RUN:1:1101:1000:2000", loc));
                Assert.assertFalse(parser.addLocationInformation("not a name", loc));
                assertUntouched(loc);
            });
            Assert.assertEquals(log, "");
        }
    }

    @Test
    public void testDefaultRegexWarnsOnce() {
        final ReadNameParser parser = new ReadNameParser(ReadNameParser.DEFAULT_READ_NAME_REGEX, Log.getInstance(ReadNameParserBehaviorTest.class));
        // numeric fields, but neither five nor seven of them
        final String log = captureLog(() -> {
            Assert.assertFalse(parser.addLocationInformation("bad:1:2:3", new PhysicalLocationInt()));
            Assert.assertFalse(parser.addLocationInformation("x:y:z:1:2:3", new PhysicalLocationInt()));
            Assert.assertFalse(parser.addLocationInformation("RUN:1:X1101:1000:2000", new PhysicalLocationInt()));
            Assert.assertTrue(parser.addLocationInformation("RUN:1:1101:1000:2000", new PhysicalLocationInt()));
            Assert.assertFalse(parser.addLocationInformation("7:8:9", new PhysicalLocationInt()));
        });
        Assert.assertEquals(occurrences(log, "did not match read name"), 1, log);
        Assert.assertTrue(log.contains("Default READ_NAME_REGEX '" + ReadNameParser.DEFAULT_READ_NAME_REGEX + "' did not match read name 'bad:1:2:3'"), log);
        Assert.assertFalse(log.contains("x:y:z:1:2:3"), log);
        Assert.assertFalse(log.contains("7:8:9"), log);
        Assert.assertFalse(log.contains("expected to contain an integer"), log);

        // the flag belongs to the parser instance
        final ReadNameParser another = new ReadNameParser(ReadNameParser.DEFAULT_READ_NAME_REGEX, Log.getInstance(ReadNameParserBehaviorTest.class));
        final String secondLog = captureLog(() -> Assert.assertFalse(another.addLocationInformation("7:8:9", new PhysicalLocationInt())));
        Assert.assertEquals(occurrences(secondLog, "did not match read name '7:8:9'"), 1, secondLog);
    }

    @Test
    public void testCustomRegexWarnsOnce() {
        final ReadNameParser parser = new ReadNameParser(LAST_THREE_FIELDS_REGEX, Log.getInstance(ReadNameParserBehaviorTest.class));
        final String log = captureLog(() -> {
            Assert.assertFalse(parser.addLocationInformation("firstBadName", new PhysicalLocationInt()));
            Assert.assertFalse(parser.addLocationInformation("secondBadName", new PhysicalLocationInt()));
            Assert.assertFalse(parser.addLocationInformation("RUN:1:99999999999:1000:2000", new PhysicalLocationInt()));
        });
        Assert.assertEquals(occurrences(log, "did not match read name"), 1, log);
        Assert.assertTrue(log.contains("READ_NAME_REGEX '" + LAST_THREE_FIELDS_REGEX + "' did not match read name 'firstBadName'"), log);
        Assert.assertFalse(log.contains("Default READ_NAME_REGEX"), log);
        Assert.assertFalse(log.contains("secondBadName"), log);
        Assert.assertFalse(log.contains("expected to contain an integer"), log);
    }

    /** A field that is not a number warns with its own message, and uses up the single warning. */
    @Test
    public void testNumberFormatWarningSharesTheFlag() {
        final ReadNameParser fast = new ReadNameParser(ReadNameParser.DEFAULT_READ_NAME_REGEX, Log.getInstance(ReadNameParserBehaviorTest.class));
        final String fastLog = captureLog(() -> {
            Assert.assertFalse(fast.addLocationInformation("RUN:1:X1101:1000:2000", new PhysicalLocationInt()));
            Assert.assertFalse(fast.addLocationInformation("RUN:1:1101:Y1000:2000", new PhysicalLocationInt()));
            Assert.assertFalse(fast.addLocationInformation("bad:1:2:3", new PhysicalLocationInt()));
        });
        Assert.assertEquals(occurrences(fastLog, "expected to contain an integer"), 1, fastLog);
        Assert.assertTrue(fastLog.contains("RUN:1:X1101:1000:2000"), fastLog);
        Assert.assertFalse(fastLog.contains("Y1000"), fastLog);
        Assert.assertFalse(fastLog.contains("did not match"), fastLog);

        final ReadNameParser regex = new ReadNameParser(LAST_THREE_FIELDS_REGEX, Log.getInstance(ReadNameParserBehaviorTest.class));
        final String regexLog = captureLog(() -> {
            Assert.assertFalse(regex.addLocationInformation("RUN:1:99999999999:1000:2000", new PhysicalLocationInt()));
            Assert.assertFalse(regex.addLocationInformation("firstBadName", new PhysicalLocationInt()));
        });
        Assert.assertEquals(occurrences(regexLog, "expected to contain an integer"), 1, regexLog);
        Assert.assertTrue(regexLog.contains("RUN:1:99999999999:1000:2000"), regexLog);
        Assert.assertFalse(regexLog.contains("did not match"), regexLog);
    }

    @Test
    public void testNoWarningWithoutLog() {
        final String log = captureLog(() -> {
            Assert.assertFalse(new ReadNameParser().addLocationInformation("badName", new PhysicalLocationInt()));
            Assert.assertFalse(new ReadNameParser(LAST_THREE_FIELDS_REGEX).addLocationInformation("badName", new PhysicalLocationInt()));
            Assert.assertFalse(new ReadNameParser().addLocationInformation("RUN:1:X:1:2", new PhysicalLocationInt()));
            Assert.assertFalse(new ReadNameParser(LAST_THREE_FIELDS_REGEX, null).addLocationInformation("RUN:1:99999999999:1:2", new PhysicalLocationInt()));
        });
        Assert.assertEquals(log, "");
    }

    @DataProvider(name = "parsers")
    public Object[][] parsers() {
        return new Object[][]{
                {ReadNameParser.DEFAULT_READ_NAME_REGEX},
                {LAST_THREE_FIELDS_REGEX},
        };
    }

    /** A repeated name is served from the last successful parse, which a failed parse does not replace. */
    @Test(dataProvider = "parsers")
    public void testRepeatedNames(final String regex) {
        final ReadNameParser parser = new ReadNameParser(regex);
        final String first = "RUN:1:1101:1000:2000";
        final String second = "RUN:1:1102:5:6";

        final PhysicalLocationInt loc1 = new PhysicalLocationInt();
        Assert.assertTrue(parser.addLocationInformation(first, loc1));
        assertLocation(loc1, 1101, 1000, 2000);

        final PhysicalLocationInt loc2 = new PhysicalLocationInt();
        Assert.assertTrue(parser.addLocationInformation(first, loc2));
        assertLocation(loc2, 1101, 1000, 2000);

        final PhysicalLocationInt loc3 = new PhysicalLocationInt();
        Assert.assertFalse(parser.addLocationInformation("badName", loc3));
        assertUntouched(loc3);

        final PhysicalLocationInt loc4 = new PhysicalLocationInt();
        Assert.assertTrue(parser.addLocationInformation(first, loc4));
        assertLocation(loc4, 1101, 1000, 2000);

        final PhysicalLocationInt loc5 = new PhysicalLocationInt();
        Assert.assertTrue(parser.addLocationInformation(second, loc5));
        assertLocation(loc5, 1102, 5, 6);

        final PhysicalLocationInt loc6 = new PhysicalLocationInt();
        Assert.assertTrue(parser.addLocationInformation(new String(first.toCharArray()), loc6));
        assertLocation(loc6, 1101, 1000, 2000);

        // the location object is overwritten, not merged
        final PhysicalLocationInt reused = new PhysicalLocationInt();
        Assert.assertTrue(parser.addLocationInformation(second, reused));
        Assert.assertTrue(parser.addLocationInformation(second, reused));
        assertLocation(reused, 1102, 5, 6);
        Assert.assertTrue(parser.addLocationInformation(first, reused));
        assertLocation(reused, 1101, 1000, 2000);
    }

    /** Only tile, x and y are written; read group and library are left alone. */
    @Test
    public void testOnlyTileAndCoordinatesAreWritten() {
        final PhysicalLocation loc = new PhysicalLocationInt() {
            short readGroup = 3;
            short libraryId = 4;

            @Override
            public short getReadGroup() { return readGroup; }

            @Override
            public void setReadGroup(final short rg) { Assert.fail("read group was set"); }

            @Override
            public short getLibraryId() { return libraryId; }

            @Override
            public void setLibraryId(final short id) { Assert.fail("library id was set"); }
        };
        final ReadNameParser parser = new ReadNameParser();
        Assert.assertTrue(parser.addLocationInformation("RUN:1:1101:1000:2000", loc));
        Assert.assertTrue(parser.addLocationInformation("RUN:1:1101:1000:2000", loc));
        assertLocation(loc, 1101, 1000, 2000);
        Assert.assertEquals(loc.getReadGroup(), 3);
        Assert.assertEquals(loc.getLibraryId(), 4);
    }

    @DataProvider(name = "lastThreeFieldsData")
    public Object[][] lastThreeFieldsData() {
        return new Object[][]{
                // name, delimiter, returned field count, tokens
                {"1:2:3", ':', 3, new int[]{1, 2, 3}},
                {"a:1:2:3", ':', 4, new int[]{1, 2, 3}},
                {"a:b:1:2:3", ':', 5, new int[]{1, 2, 3}},
                {"a:b:c:1:2:3", ':', 6, new int[]{1, 2, 3}},
                {"a:b:c:d:1:2:3", ':', 7, new int[]{1, 2, 3}},
                {"a:b:c:d:e:1:2:3", ':', 8, new int[]{1, 2, 3}},
                {"a::1:2:3", ':', 5, new int[]{1, 2, 3}},
                {"x_1101_-20_30abc", '_', 4, new int[]{1101, -20, 30}},
                {"1101:20:30", '_', -1, new int[]{-1, -1, -1}},
                {"1:2", ':', -1, new int[]{-1, -1, -1}},
                {"7", ':', -1, new int[]{-1, -1, -1}},
                {"", ':', -1, new int[]{-1, -1, -1}},
        };
    }

    @Test(dataProvider = "lastThreeFieldsData")
    public void testGetLastThreeFields(final String name, final char delimiter, final int fields, final int[] tokens) {
        final int[] actual = {99, 99, 99};
        Assert.assertEquals(ReadNameParser.getLastThreeFields(name, delimiter, actual), fields);
        Assert.assertEquals(actual, tokens);
    }

    @Test
    public void testGetLastThreeFieldsRejectsNonNumericFields() {
        for (final String name : new String[]{"a:b:1:2:", "a:b:1::3", "a:b:x:2:3", "a:b:1:-:3", "x:2"}) {
            Assert.assertThrows(NumberFormatException.class, () -> ReadNameParser.getLastThreeFields(name, ':', new int[3]));
        }
        // fields before the last three are never parsed
        Assert.assertEquals(ReadNameParser.getLastThreeFields("x:y:z:1:2:3", ':', new int[3]), 6);
    }

    @DataProvider(name = "rapidParseIntData")
    public Object[][] rapidParseIntData() {
        return new Object[][]{
                {"0", 0},
                {"-0", 0},
                {"7", 7},
                {"0042", 42},
                {"123abc456", 123},
                {"-17:5", -17},
                {"2147483647", Integer.MAX_VALUE},
                {"-2147483648", Integer.MIN_VALUE},
                {"2147483646x", Integer.MAX_VALUE - 1},
        };
    }

    @Test(dataProvider = "rapidParseIntData")
    public void testRapidParseInt(final String input, final int expected) {
        Assert.assertEquals(ReadNameParser.rapidParseInt(input), expected);
        Assert.assertEquals(ReadNameParser.rapidParseInt(input), Integer.parseInt(input.replaceAll("(?<=[0-9])[^0-9].*$", "")));
    }

    @Test
    public void testRapidParseIntFailures() {
        for (final String input : new String[]{"", "-", "--1", "+1", " 1", "x1", "-x1"}) {
            final NumberFormatException e = Assert.expectThrows(NumberFormatException.class, () -> ReadNameParser.rapidParseInt(input));
            Assert.assertTrue(e.getMessage().contains("'" + input + "'"), e.getMessage());
        }
    }
}
