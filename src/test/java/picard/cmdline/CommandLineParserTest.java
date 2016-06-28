/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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
package picard.cmdline;

import htsjdk.samtools.util.CollectionUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.programgroups.Testing;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class CommandLineParserTest {

    enum FrobnicationFlavor {
        FOO, BAR, BAZ
    }

    @CommandLineProgramProperties(
            usage = "Usage: frobnicate [options] input-file output-file\n\nRead input-file, frobnicate it, and write frobnicated results to output-file\n",
            usageShort = "Read input-file, frobnicate it, and write frobnicated results to output-file"
    )
    class FrobnicateOptions {

        @PositionalArguments(minElements = 2, maxElements = 2)
        public List<File> positionalArguments = new ArrayList<File>();

        @Option(shortName = "T", doc = "Frobnication threshold setting.")
        public Integer FROBNICATION_THRESHOLD = 20;

        @Option
        public FrobnicationFlavor FROBNICATION_FLAVOR;

        @Option(doc = "Allowed shmiggle types.", minElements = 1, maxElements = 3)
        public List<String> SHMIGGLE_TYPE = new ArrayList<String>();

        @Option
        public Boolean TRUTHINESS;
    }

    @CommandLineProgramProperties(
            usage = "Usage: frobnicate [options] input-file output-file\n\nRead input-file, frobnicate it, and write frobnicated results to output-file\n",
            usageShort = "Read input-file, frobnicate it, and write frobnicated results to output-file"
    )
    class FrobnicateOptionsWithNullList {

        @PositionalArguments(minElements = 2, maxElements = 2)
        public List<File> positionalArguments = new ArrayList<File>();

        @Option(shortName = "T", doc = "Frobnication threshold setting.")
        public Integer FROBNICATION_THRESHOLD = 20;

        @Option
        public FrobnicationFlavor FROBNICATION_FLAVOR;

        @Option(doc = "Allowed shmiggle types.", minElements = 0, maxElements = 3)
        public List<String> SHMIGGLE_TYPE = new ArrayList<String>();

        @Option
        public Boolean TRUTHINESS;
    }

    @CommandLineProgramProperties(
            usage = "Usage: framistat [options]\n\nCompute the plebnick of the freebozzle.\n",
            usageShort = "ompute the plebnick of the freebozzle"
    )
    class OptionsWithoutPositional {
        public static final int DEFAULT_FROBNICATION_THRESHOLD = 20;
        @Option(shortName = "T", doc = "Frobnication threshold setting.")
        public Integer FROBNICATION_THRESHOLD = DEFAULT_FROBNICATION_THRESHOLD;

        @Option
        public FrobnicationFlavor FROBNICATION_FLAVOR;

        @Option(doc = "Allowed shmiggle types.", minElements = 1, maxElements = 3)
        public List<String> SHMIGGLE_TYPE = new ArrayList<String>();

        @Option
        public Boolean TRUTHINESS;
    }

    class OptionsWithCaseClash {
        @Option
        public String FROB;
        @Option
        public String frob;
    }

    class OptionsWithSameShortName {
        @Option(shortName = "SAME_SHORT_NAME", overridable = true, optional = true)
        public String SAME_SHORT_NAME;
        @Option(shortName = "SOMETHING_ELSE", overridable = true, optional = true)
        public String DIFF_SHORT_NAME;
    }

    class MutexOptions {
        @Option(mutex = {"M", "N", "Y", "Z"})
        public String A;
        @Option(mutex = {"M", "N", "Y", "Z"})
        public String B;
        @Option(mutex = {"A", "B", "Y", "Z"})
        public String M;
        @Option(mutex = {"A", "B", "Y", "Z"})
        public String N;
        @Option(mutex = {"A", "B", "M", "N"})
        public String Y;
        @Option(mutex = {"A", "B", "M", "N"})
        public String Z;

    }


    @Test
    public void testUsage() {
        final FrobnicateOptions fo = new FrobnicateOptions();
        final CommandLineParser clp = new CommandLineParser(fo);
        clp.usage(System.out, false);
    }

    @Test
    public void testUsageWithDefault() {
        final FrobnicateOptions fo = new FrobnicateOptions();
        final CommandLineParser clp = new CommandLineParser(fo);
        clp.usage(System.out, true);
    }

    @Test
    public void testUsageWithoutPositional() {
        final OptionsWithoutPositional fo = new OptionsWithoutPositional();
        final CommandLineParser clp = new CommandLineParser(fo);
        clp.usage(System.out, false);
    }

    @Test
    public void testUsageWithoutPositionalWithDefault() {
        final OptionsWithoutPositional fo = new OptionsWithoutPositional();
        final CommandLineParser clp = new CommandLineParser(fo);
        clp.usage(System.out, true);
    }

    /**
     * If the short name is set to be the same as the long name we still want the argument to appear in the commandLine.
     */
    @Test
    public void testForIdenticalShortName() {
        final String[] args = {
                "SAME_SHORT_NAME=FOO",
                "SOMETHING_ELSE=BAR"
        };
        final OptionsWithSameShortName fo = new OptionsWithSameShortName();
        final CommandLineParser clp = new CommandLineParser(fo);
        clp.parseOptions(System.err, args);
        final String commandLine = clp.getCommandLine();
        Assert.assertTrue(commandLine.contains("DIFF_SHORT_NAME"));
        Assert.assertTrue(commandLine.contains("SAME_SHORT_NAME"));
    }


    @Test
    public void testPositive() {
        final String[] args = {
                "T=17",
                "FROBNICATION_FLAVOR=BAR",
                "TRUTHINESS=False",
                "SHMIGGLE_TYPE=shmiggle1",
                "SHMIGGLE_TYPE=shmiggle2",
                "positional1",
                "positional2",
        };
        final FrobnicateOptions fo = new FrobnicateOptions();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertTrue(clp.parseOptions(System.err, args));
        Assert.assertEquals(fo.positionalArguments.size(), 2);
        final File[] expectedPositionalArguments = {new File("positional1"), new File("positional2")};
        Assert.assertEquals(fo.positionalArguments.toArray(), expectedPositionalArguments);
        Assert.assertEquals(fo.FROBNICATION_THRESHOLD.intValue(), 17);
        Assert.assertEquals(fo.FROBNICATION_FLAVOR, FrobnicationFlavor.BAR);
        Assert.assertEquals(fo.SHMIGGLE_TYPE.size(), 2);
        final String[] expectedShmiggleTypes = {"shmiggle1", "shmiggle2"};
        Assert.assertEquals(fo.SHMIGGLE_TYPE.toArray(), expectedShmiggleTypes);
        Assert.assertFalse(fo.TRUTHINESS);
    }

    /**
     * Allow a whitespace btw equal sign and option value.
     */
    @Test
    public void testPositiveWithSpaces() {
        final String[] args = {
                "T=", "17",
                "FROBNICATION_FLAVOR=", "BAR",
                "TRUTHINESS=", "False",
                "SHMIGGLE_TYPE=", "shmiggle1",
                "SHMIGGLE_TYPE=", "shmiggle2",
                "positional1",
                "positional2",
        };
        final FrobnicateOptions fo = new FrobnicateOptions();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertTrue(clp.parseOptions(System.err, args));
        Assert.assertEquals(fo.positionalArguments.size(), 2);
        final File[] expectedPositionalArguments = {new File("positional1"), new File("positional2")};
        Assert.assertEquals(fo.positionalArguments.toArray(), expectedPositionalArguments);
        Assert.assertEquals(fo.FROBNICATION_THRESHOLD.intValue(), 17);
        Assert.assertEquals(fo.FROBNICATION_FLAVOR, FrobnicationFlavor.BAR);
        Assert.assertEquals(fo.SHMIGGLE_TYPE.size(), 2);
        final String[] expectedShmiggleTypes = {"shmiggle1", "shmiggle2"};
        Assert.assertEquals(fo.SHMIGGLE_TYPE.toArray(), expectedShmiggleTypes);
        Assert.assertFalse(fo.TRUTHINESS);
    }

    @Test
    public void testPositiveWithoutPositional() {
        final String[] args = {
                "T=17",
                "FROBNICATION_FLAVOR=BAR",
                "TRUTHINESS=False",
                "SHMIGGLE_TYPE=shmiggle1",
                "SHMIGGLE_TYPE=shmiggle2",
        };
        final OptionsWithoutPositional fo = new OptionsWithoutPositional();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertTrue(clp.parseOptions(System.err, args));
        Assert.assertEquals(fo.FROBNICATION_THRESHOLD.intValue(), 17);
        Assert.assertEquals(fo.FROBNICATION_FLAVOR, FrobnicationFlavor.BAR);
        Assert.assertEquals(fo.SHMIGGLE_TYPE.size(), 2);
        final String[] expectedShmiggleTypes = {"shmiggle1", "shmiggle2"};
        Assert.assertEquals(fo.SHMIGGLE_TYPE.toArray(), expectedShmiggleTypes);
        Assert.assertFalse(fo.TRUTHINESS);
    }

    /**
     * If last character of command line is the equal sign in an option=value pair,
     * make sure no crash, and that the value is empty string.
     */
    @Test
    public void testPositiveTerminalEqualSign() {
        final String[] args = {
                "T=17",
                "FROBNICATION_FLAVOR=BAR",
                "TRUTHINESS=False",
                "SHMIGGLE_TYPE=shmiggle1",
                "SHMIGGLE_TYPE=",
        };
        final OptionsWithoutPositional fo = new OptionsWithoutPositional();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertTrue(clp.parseOptions(System.err, args));
        Assert.assertEquals(fo.FROBNICATION_THRESHOLD.intValue(), 17);
        Assert.assertEquals(fo.FROBNICATION_FLAVOR, FrobnicationFlavor.BAR);
        Assert.assertEquals(fo.SHMIGGLE_TYPE.size(), 2);
        final String[] expectedShmiggleTypes = {"shmiggle1", ""};
        Assert.assertEquals(fo.SHMIGGLE_TYPE.toArray(), expectedShmiggleTypes);
        Assert.assertFalse(fo.TRUTHINESS);
    }

    @Test
    public void testDefault() {
        final String[] args = {
                "FROBNICATION_FLAVOR=BAR",
                "TRUTHINESS=False",
                "SHMIGGLE_TYPE=shmiggle1",
                "SHMIGGLE_TYPE=shmiggle2",
                "positional1",
                "positional2",
        };
        final FrobnicateOptions fo = new FrobnicateOptions();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertTrue(clp.parseOptions(System.err, args));
        Assert.assertEquals(fo.FROBNICATION_THRESHOLD.intValue(), 20);
    }

    @Test
    public void testMissingRequiredArgument() {
        final String[] args = {
                "TRUTHINESS=False",
                "SHMIGGLE_TYPE=shmiggle1",
                "SHMIGGLE_TYPE=shmiggle2",
                "positional1",
                "positional2",
        };
        final FrobnicateOptions fo = new FrobnicateOptions();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertFalse(clp.parseOptions(System.err, args));
    }

    @Test
    public void testBadValue() {
        final String[] args = {
                "FROBNICATION_THRESHOLD=ABC",
                "FROBNICATION_FLAVOR=BAR",
                "TRUTHINESS=False",
                "SHMIGGLE_TYPE=shmiggle1",
                "SHMIGGLE_TYPE=shmiggle2",
                "positional1",
                "positional2",
        };
        final FrobnicateOptions fo = new FrobnicateOptions();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertFalse(clp.parseOptions(System.err, args));
    }

    @Test
    public void testBadEnumValue() {
        final String[] args = {
                "FROBNICATION_FLAVOR=HiMom",
                "TRUTHINESS=False",
                "SHMIGGLE_TYPE=shmiggle1",
                "SHMIGGLE_TYPE=shmiggle2",
                "positional1",
                "positional2",
        };
        final FrobnicateOptions fo = new FrobnicateOptions();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertFalse(clp.parseOptions(System.err, args));
    }

    @Test
    public void testNotEnoughOfListOption() {
        final String[] args = {
                "FROBNICATION_FLAVOR=BAR",
                "TRUTHINESS=False",
                "positional1",
                "positional2",
        };
        final FrobnicateOptions fo = new FrobnicateOptions();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertFalse(clp.parseOptions(System.err, args));
    }

    @Test
    public void testTooManyListOption() {
        final String[] args = {
                "FROBNICATION_FLAVOR=BAR",
                "TRUTHINESS=False",
                "SHMIGGLE_TYPE=shmiggle1",
                "SHMIGGLE_TYPE=shmiggle2",
                "SHMIGGLE_TYPE=shmiggle3",
                "SHMIGGLE_TYPE=shmiggle4",
                "positional1",
                "positional2",
        };
        final FrobnicateOptions fo = new FrobnicateOptions();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertFalse(clp.parseOptions(System.err, args));
    }

    @Test
    public void testTooManyPositional() {
        final String[] args = {
                "FROBNICATION_FLAVOR=BAR",
                "TRUTHINESS=False",
                "SHMIGGLE_TYPE=shmiggle1",
                "SHMIGGLE_TYPE=shmiggle2",
                "positional1",
                "positional2",
                "positional3",
        };
        final FrobnicateOptions fo = new FrobnicateOptions();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertFalse(clp.parseOptions(System.err, args));
    }

    @Test
    public void testNotEnoughPositional() {
        final String[] args = {
                "FROBNICATION_FLAVOR=BAR",
                "TRUTHINESS=False",
                "SHMIGGLE_TYPE=shmiggle1",
                "SHMIGGLE_TYPE=shmiggle2",
        };
        final FrobnicateOptions fo = new FrobnicateOptions();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertFalse(clp.parseOptions(System.err, args));
    }

    @Test
    public void testUnexpectedPositional() {
        final String[] args = {
                "T=17",
                "FROBNICATION_FLAVOR=BAR",
                "TRUTHINESS=False",
                "SHMIGGLE_TYPE=shmiggle1",
                "SHMIGGLE_TYPE=shmiggle2",
                "positional"
        };
        final OptionsWithoutPositional fo = new OptionsWithoutPositional();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertFalse(clp.parseOptions(System.err, args));
    }

    @Test(expectedExceptions = CommandLineParserDefinitionException.class)
    public void testOptionDefinitionCaseClash() {
        final OptionsWithCaseClash options = new OptionsWithCaseClash();
        new CommandLineParser(options);
        Assert.fail("Should not be reached.");
    }

    @Test
    public void testOptionUseCaseClash() {
        final String[] args = {
                "FROBNICATION_FLAVOR=BAR",
                "FrOBNICATION_fLAVOR=BAR",
        };
        final FrobnicateOptions fo = new FrobnicateOptions();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertFalse(clp.parseOptions(System.err, args));
    }

    @Test
    public void testNullValue() {
        final String[] args = {
                "FROBNICATION_THRESHOLD=null",
                "FROBNICATION_FLAVOR=BAR",
                "TRUTHINESS=False",
                "SHMIGGLE_TYPE=null",
                "positional1",
                "positional2",
        };

        final FrobnicateOptionsWithNullList fownl = new FrobnicateOptionsWithNullList();
        fownl.SHMIGGLE_TYPE.add("shmiggle1"); //providing null value should clear this list

        final CommandLineParser clp = new CommandLineParser(fownl);
        Assert.assertTrue(clp.parseOptions(System.err, args));
        Assert.assertEquals(fownl.positionalArguments.size(), 2);
        final File[] expectedPositionalArguments = {new File("positional1"), new File("positional2")};
        Assert.assertEquals(fownl.positionalArguments.toArray(), expectedPositionalArguments);
        Assert.assertEquals(fownl.FROBNICATION_THRESHOLD, null); //test null value         
        Assert.assertEquals(fownl.SHMIGGLE_TYPE.size(), 0); //test null value for list        
        Assert.assertFalse(fownl.TRUTHINESS);

        //verify that required arg can't be set to null
        args[2] = "TRUTHINESS=null";
        final CommandLineParser clp2 = new CommandLineParser(fownl);
        Assert.assertFalse(clp2.parseOptions(System.err, args));

        //verify that positional arg can't be set to null
        args[2] = "TRUTHINESS=False";
        args[4] = "null";
        final CommandLineParser clp3 = new CommandLineParser(fownl);
        Assert.assertFalse(clp3.parseOptions(System.err, args));

    }


    @Test
    public void testOptionsFile() throws Exception {
        final File optionsFile = File.createTempFile("clp.", ".options");
        optionsFile.deleteOnExit();
        final PrintWriter writer = new PrintWriter(optionsFile);
        writer.println("T=18");
        writer.println("TRUTHINESS=True");
        writer.println("SHMIGGLE_TYPE=shmiggle0");
        writer.println("STRANGE_OPTION=shmiggle0");
        writer.close();
        final String[] args = {
                "OPTIONS_FILE=" + optionsFile.getPath(),
                // Multiple options files are allowed
                "OPTIONS_FILE=" + optionsFile.getPath(),
                "T=17",
                "FROBNICATION_FLAVOR=BAR",
                "TRUTHINESS=False",
                "SHMIGGLE_TYPE=shmiggle1",
                "positional1",
                "positional2",
        };
        final FrobnicateOptions fo = new FrobnicateOptions();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertTrue(clp.parseOptions(System.err, args));
        Assert.assertEquals(fo.positionalArguments.size(), 2);
        final File[] expectedPositionalArguments = {new File("positional1"), new File("positional2")};
        Assert.assertEquals(fo.positionalArguments.toArray(), expectedPositionalArguments);
        Assert.assertEquals(fo.FROBNICATION_THRESHOLD.intValue(), 17);
        Assert.assertEquals(fo.FROBNICATION_FLAVOR, FrobnicationFlavor.BAR);
        Assert.assertEquals(fo.SHMIGGLE_TYPE.size(), 3);
        final String[] expectedShmiggleTypes = {"shmiggle0", "shmiggle0", "shmiggle1"};
        Assert.assertEquals(fo.SHMIGGLE_TYPE.toArray(), expectedShmiggleTypes);
        Assert.assertFalse(fo.TRUTHINESS);
    }


    /**
     * In an options file, should not be allowed to override an option set on the command line
     *
     * @throws Exception
     */
    @Test
    public void testOptionsFileWithDisallowedOverride() throws Exception {
        final File optionsFile = File.createTempFile("clp.", ".options");
        optionsFile.deleteOnExit();
        final PrintWriter writer = new PrintWriter(optionsFile);
        writer.println("T=18");
        writer.close();
        final String[] args = {
                "T=17",
                "OPTIONS_FILE=" + optionsFile.getPath()
        };
        final FrobnicateOptions fo = new FrobnicateOptions();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertFalse(clp.parseOptions(System.err, args));
    }

    @DataProvider(name = "mutexScenarios")
    public Object[][] mutexScenarios() {
        return new Object[][]{
                {"pass", new String[]{"A=1", "B=2"}, true},
                {"no args", new String[0], false},
                {"1 of group required", new String[]{"A=1"}, false},
                {"mutex", new String[]{"A=1", "Y=3"}, false},
                {"mega mutex", new String[]{"A=1", "B=2", "Y=3", "Z=1", "M=2", "N=3"}, false}
        };
    }

    @Test(dataProvider = "mutexScenarios")
    public void testMutex(final String testName, final String[] args, final boolean expected) {
        final MutexOptions o = new MutexOptions();
        final CommandLineParser clp = new CommandLineParser(o);
        Assert.assertEquals(clp.parseOptions(System.err, args), expected);
    }

    class UninitializedCollectionOptions {
        @Option
        public List<String> LIST;
        @Option
        public ArrayList<String> ARRAY_LIST;
        @Option
        public HashSet<String> HASH_SET;
        @PositionalArguments
        public Collection<File> COLLECTION;

    }

    @Test
    public void testUninitializedCollections() {
        final UninitializedCollectionOptions o = new UninitializedCollectionOptions();
        final CommandLineParser clp = new CommandLineParser(o);
        final String[] args = {"LIST=L1", "LIST=L2", "ARRAY_LIST=S1", "HASH_SET=HS1", "P1", "P2"};
        Assert.assertTrue(clp.parseOptions(System.err, args));
        Assert.assertEquals(o.LIST.size(), 2);
        Assert.assertEquals(o.ARRAY_LIST.size(), 1);
        Assert.assertEquals(o.HASH_SET.size(), 1);
        Assert.assertEquals(o.COLLECTION.size(), 2);
    }

    class UninitializedCollectionThatCannotBeAutoInitializedOptions {
        @Option
        public Set<String> SET;
    }

    @Test(expectedExceptions = CommandLineParserDefinitionException.class)
    public void testCollectionThatCannotBeAutoInitialized() {
        final UninitializedCollectionThatCannotBeAutoInitializedOptions o = new UninitializedCollectionThatCannotBeAutoInitializedOptions();
        new CommandLineParser(o);
        Assert.fail("Exception should have been thrown");
    }

    class CollectionWithDefaultValuesOptions {
        @Option
        public List<String> LIST = CollectionUtil.makeList("foo", "bar");
    }

    @Test
    public void testClearDefaultValuesFromListOption() {
        final CollectionWithDefaultValuesOptions o = new CollectionWithDefaultValuesOptions();
        final CommandLineParser clp = new CommandLineParser(o);
        final String[] args = {"LIST=null"};
        Assert.assertTrue(clp.parseOptions(System.err, args));
        Assert.assertEquals(o.LIST.size(), 0);
    }

    @Test
    public void testClearDefaultValuesFromListOptionAndAddNew() {
        final CollectionWithDefaultValuesOptions o = new CollectionWithDefaultValuesOptions();
        final CommandLineParser clp = new CommandLineParser(o);
        final String[] args = {"LIST=null", "LIST=baz", "LIST=frob"};
        Assert.assertTrue(clp.parseOptions(System.err, args));
        Assert.assertEquals(o.LIST, CollectionUtil.makeList("baz", "frob"));
    }

    @Test
    public void testAddToDefaultValuesListOption() {
        final CollectionWithDefaultValuesOptions o = new CollectionWithDefaultValuesOptions();
        final CommandLineParser clp = new CommandLineParser(o);
        final String[] args = {"LIST=baz", "LIST=frob"};
        Assert.assertTrue(clp.parseOptions(System.err, args));
        Assert.assertEquals(o.LIST, CollectionUtil.makeList("foo", "bar", "baz", "frob"));
    }

    @CommandLineProgramProperties(
            usage = "Class with nested option",
            usageShort = "Class with nested option"
    )
    class OptionsWithNested {
        @Option
        public Integer AN_INT;
        @NestedOptions(doc = "Doc for FROB")
        public OptionsWithoutPositional FROB = new OptionsWithoutPositional();
        @NestedOptions
        public OptionsWithNestedAgain NESTED = new OptionsWithNestedAgain();
        @Option
        public String A_STRING;
    }

    class OptionsWithNestedAgain {
        @NestedOptions(doc = "Doc for inner FROB")
        public OptionsWithoutPositional FROB = new OptionsWithoutPositional();
    }


    @Test
    public void testStaticNestedOptions() {
        final OptionsWithNested o = new OptionsWithNested();
        final CommandLineParser clp = new CommandLineParser(o);
        clp.usage(System.out, false);
        clp.htmlUsage(System.out, "testStaticNestedOptions", false);
        final int outerInt = 123;
        final String outerString = "outerString";
        final FrobnicationFlavor outerFlavor = FrobnicationFlavor.BAR;
        final FrobnicationFlavor innerFlavor = FrobnicationFlavor.BAZ;
        final boolean outerTruthiness = true;
        final boolean innerTruthiness = false;
        final int innerThreshold = 10;
        final String[] outerShmiggle = {"shmiggle1", "shmiggle2"};
        final String[] innerShmiggle = {"inner1", "inner2", "inner3"};

        final List<String> args = new ArrayList<String>();
        args.add("AN_INT=" + outerInt);
        args.add("A_STRING=" + outerString);
        args.add("frob.truThIness=" + outerTruthiness);
        args.add("FrOb.FROBNICATION_FLAVOR=" + outerFlavor);
        for (final String shmiggle : outerShmiggle) {
            args.add("FROB.SHMIGGLE_TYPE=" + shmiggle);
        }
        args.add("NeStEd.Frob.FROBNICATION_THRESHOLD=" + innerThreshold);
        args.add("NeStEd.Frob.FROBNICATION_FLAVOR=" + innerFlavor);
        args.add("NeStEd.Frob.truthIness=" + innerTruthiness);
        for (final String shmiggle : innerShmiggle) {
            args.add("NESTED.FROB.SHMIGGLE_TYPE=");
            args.add(shmiggle);
        }
        Assert.assertTrue(clp.parseOptions(System.err, args.toArray(new String[args.size()])));
        System.out.println(clp.getCommandLine());
        Assert.assertEquals(o.AN_INT.intValue(), outerInt);
        Assert.assertEquals(o.A_STRING, outerString);
        Assert.assertEquals(o.FROB.FROBNICATION_FLAVOR, outerFlavor);
        Assert.assertEquals(o.FROB.FROBNICATION_THRESHOLD.intValue(), OptionsWithoutPositional.DEFAULT_FROBNICATION_THRESHOLD);
        Assert.assertEquals(o.FROB.SHMIGGLE_TYPE, Arrays.asList(outerShmiggle));
        Assert.assertEquals(o.FROB.TRUTHINESS.booleanValue(), outerTruthiness);
        Assert.assertEquals(o.NESTED.FROB.SHMIGGLE_TYPE, Arrays.asList(innerShmiggle));
        Assert.assertEquals(o.NESTED.FROB.FROBNICATION_THRESHOLD.intValue(), innerThreshold);
        Assert.assertEquals(o.NESTED.FROB.FROBNICATION_FLAVOR, innerFlavor);
        Assert.assertEquals(o.NESTED.FROB.TRUTHINESS.booleanValue(), innerTruthiness);
    }

    @Test
    void testStaticNestedNegative() {
        final OptionsWithNested o = new OptionsWithNested();
        final CommandLineParser clp = new CommandLineParser(o);
        final int outerInt = 123;
        final String outerString = "outerString";
        final FrobnicationFlavor outerFlavor = FrobnicationFlavor.BAR;
        final boolean outerTruthiness = true;
        final String[] outerShmiggle = {"shmiggle1", "shmiggle2"};

        final List<String> args = new ArrayList<String>();
        args.add("AN_INT=" + outerInt);
        args.add("A_STRING=" + outerString);
        Assert.assertFalse(clp.parseOptions(System.err, args.toArray(new String[args.size()])));
        System.out.println(clp.getCommandLine());
    }

    @CommandLineProgramProperties(
            usage = "Class with nested options.",
            usageShort = "Class with nested options",
            programGroup = Testing.class,
            omitFromCommandLine = true
    )
    class ClpOptionsWithNested extends CommandLineProgram {
        @Option
        public Integer AN_INT;
        @NestedOptions(doc = "This will be ignored")
        public OptionsWithoutPositional FROB = new OptionsWithoutPositional();

        @Option
        public String A_STRING;

        private final ClpOptionsWithNestedAgain NESTED = new ClpOptionsWithNestedAgain();

        @Override
        public Map<String, Object> getNestedOptions() {
            final Map<String, Object> ret = new LinkedHashMap<String, Object>();
            ret.put("CLP_NESTED", NESTED);
            ret.put("FRAB", FROB);
            return ret;
        }

        @Override
        public Map<String, Object> getNestedOptionsForHelp() {
            final Map<String, Object> ret = new LinkedHashMap<String, Object>();
            ret.put("CLP_NESTED", NESTED);
            return ret;
        }

        @Override
        protected int doWork() {
            return 0;
        }
    }

    @CommandLineProgramProperties(
            usage = "Class with nested options again.",
            usageShort = "Class with nested options again",
            programGroup = Testing.class,
            omitFromCommandLine = true
    )
    class ClpOptionsWithNestedAgain extends CommandLineProgram {
        private final OptionsWithoutPositional FROB = new OptionsWithoutPositional();

        @Override
        public Map<String, Object> getNestedOptions() {
            final Map<String, Object> ret = new LinkedHashMap<String, Object>();
            ret.put("FROB_NESTED", FROB);
            return ret;
        }

        @Override
        protected int doWork() {
            return 0;
        }
    }

    @Test
    public void testDynamicNestedOptions() {
        final ClpOptionsWithNested o = new ClpOptionsWithNested();
        final int outerInt = 123;
        final String outerString = "aString";
        final int outerThreshold = 456;
        final FrobnicationFlavor outerFlavor = FrobnicationFlavor.FOO;
        final List<String> outerShmiggleType = Arrays.asList("shmiggle1");
        final boolean outerTruthiness = true;
        final int innerThreshold = -1000;
        final FrobnicationFlavor innerFlavor = FrobnicationFlavor.BAZ;
        final List<String> innerShmiggleType = Arrays.asList("innershmiggle1", "skeezwitz");
        final boolean innerTruthiness = false;

        final List<String> args = new ArrayList<String>();
        args.add("AN_INT=" + outerInt);
        args.add("A_STRING=" + outerString);
        args.add("FRAB.FROBNICATION_THRESHOLD=" + outerThreshold);
        args.add("FRAB.FROBNICATION_FLAVOR=" + outerFlavor);
        args.add("FRAB.SHMIGGLE_TYPE=" + outerShmiggleType.get(0));
        args.add("FRAB.TRUTHINESS=" + outerTruthiness);
        args.add("CLP_NESTED.FROB_NESTED.FROBNICATION_THRESHOLD=" + innerThreshold);
        args.add("CLP_NESTED.FROB_NESTED.FROBNICATION_FLAVOR=" + innerFlavor);
        for (final String ist : innerShmiggleType) {
            args.add("CLP_NESTED.FROB_NESTED.SHMIGGLE_TYPE=" + ist);
        }
        args.add("CLP_NESTED.FROB_NESTED.TRUTHINESS=" + innerTruthiness);

        Assert.assertTrue(o.parseArgs(args.toArray(new String[args.size()])));
        System.out.println(o.getCommandLine());
        Assert.assertEquals(o.AN_INT.intValue(), outerInt);
        Assert.assertEquals(o.A_STRING, outerString);
        Assert.assertEquals(o.FROB.FROBNICATION_THRESHOLD.intValue(), outerThreshold);
        Assert.assertEquals(o.FROB.FROBNICATION_FLAVOR, outerFlavor);
        Assert.assertEquals(o.FROB.SHMIGGLE_TYPE, outerShmiggleType);
        Assert.assertEquals(o.FROB.TRUTHINESS.booleanValue(), outerTruthiness);
        Assert.assertEquals(o.NESTED.FROB.FROBNICATION_THRESHOLD.intValue(), innerThreshold);
        Assert.assertEquals(o.NESTED.FROB.FROBNICATION_FLAVOR, innerFlavor);
        Assert.assertEquals(o.NESTED.FROB.SHMIGGLE_TYPE, innerShmiggleType);
        Assert.assertEquals(o.NESTED.FROB.TRUTHINESS.booleanValue(), innerTruthiness);
        Assert.assertFalse(new ClpOptionsWithNested().parseArgs(new String[]{"-h"}));
        new CommandLineParser(o).htmlUsage(System.err, o.getClass().getSimpleName(), false);
    }

    class StaticPropagationParent {
        @Option
        public String STRING1 = "String1ParentDefault";

        @Option
        public String STRING2 = "String2ParentDefault";

        @Option(overridable = true)
        public String STRING3 = "String3ParentDefault";

        @Option
        public String STRING4 = "String4ParentDefault";

        @Option
        public String STRING5;

        @Option
        public List<String> COLLECTION;

        @NestedOptions
        public PropagationChild CHILD = new PropagationChild();
    }

    class PropagationChild {
        // Parent has default, child does not, should propagate
        @Option
        public String STRING1;

        // Parent and child have default, should not propagate
        @Option
        public String STRING2 = "String2ChildDefault";

        // Parent has explicitly set value, child has default, should propagate
        @Option
        public String STRING3 = "String3ChildDefault";

        // Parent has default, child has explicitly set value, should not propagate
        @Option
        public String STRING4;

        // Parent and child have explicitly set value, should not propagate
        @Option
        public String STRING5;

        // Parent has explicitly set value, but collection should not propagate
        @Option
        public List<String> COLLECTION;
    }

    @Test
    public void testStaticPropagation() {
        final StaticPropagationParent o = new StaticPropagationParent();
        final CommandLineParser clp = new CommandLineParser(o);
        clp.usage(System.out, false);
        clp.htmlUsage(System.out, "testStaticPropagation", false);

        final List<String> args = new ArrayList<String>();
        args.add("STRING3=String3Parent");
        args.add("CHILD.STRING4=String4Child");
        args.add("STRING5=String5Parent");
        args.add("CHILD.STRING5=String5Child");
        args.add("COLLECTION=CollectionParent");

        Assert.assertTrue(clp.parseOptions(System.err, args.toArray(new String[args.size()])));
        System.out.println(clp.getCommandLine());

        Assert.assertEquals(o.CHILD.STRING1, "String1ParentDefault");
        Assert.assertEquals(o.CHILD.STRING2, "String2ChildDefault");
        Assert.assertEquals(o.CHILD.STRING3, "String3Parent");
        Assert.assertEquals(o.CHILD.STRING4, "String4Child");
        Assert.assertEquals(o.CHILD.STRING5, "String5Child");
        Assert.assertEquals(o.CHILD.COLLECTION, new ArrayList<String>());
    }

    @CommandLineProgramProperties(
            usage = "",
            usageShort = "",
            programGroup = Testing.class,
            omitFromCommandLine = true
    )
    class DynamicPropagationParent extends CommandLineProgram {
        @Option
        public String STRING1 = "String1ParentDefault";

        @Option
        public String STRING2 = "String2ParentDefault";

        @Option
        public String STRING3 = "String3ParentDefault";

        @Option
        public String STRING4 = "String4ParentDefault";

        @Option
        public String STRING5;

        @Option
        public List<String> COLLECTION;

        public PropagationChild CHILD = new PropagationChild();

        @Override
        protected int doWork() {
            return 0;
        }

        @Override
        public Map<String, Object> getNestedOptions() {
            final Map<String, Object> ret = new HashMap<String, Object>();
            ret.put("CHILD", CHILD);
            return ret;
        }
    }

    @Test
    public void testDynamicPropagation() {
        final DynamicPropagationParent o = new DynamicPropagationParent();

        final List<String> args = new ArrayList<String>();
        args.add("STRING3=String3Parent");
        args.add("CHILD.STRING4=String4Child");
        args.add("STRING5=String5Parent");
        args.add("CHILD.STRING5=String5Child");
        args.add("COLLECTION=CollectionParent");

        Assert.assertTrue(o.parseArgs(args.toArray(new String[args.size()])));
        System.out.println(o.getCommandLine());
        Assert.assertFalse(new DynamicPropagationParent().parseArgs(new String[]{"-h"}));
        new CommandLineParser(o).htmlUsage(System.err, o.getClass().getSimpleName(), false);

        Assert.assertEquals(o.CHILD.STRING1, "String1ParentDefault");
        Assert.assertEquals(o.CHILD.STRING2, "String2ChildDefault");
        Assert.assertEquals(o.CHILD.STRING3, "String3Parent");
        Assert.assertEquals(o.CHILD.STRING4, "String4Child");
        Assert.assertEquals(o.CHILD.STRING5, "String5Child");
        Assert.assertEquals(o.CHILD.COLLECTION, new ArrayList<String>());
    }

    class NegativePropagationParent {
        @Option
        public int STRING1 = 1;

        @Option
        public List<String> COLLECTION;

        @NestedOptions
        public PropagationChild CHILD = new PropagationChild();
    }

    /** parent and child option of the same name are not type-compatible. */
    @Test(expectedExceptions = {IllegalArgumentException.class})
    public void testStaticPropagationNegative() {
        final NegativePropagationParent o = new NegativePropagationParent();
        final CommandLineParser clp = new CommandLineParser(o);
        clp.usage(System.out, false);

        clp.parseOptions(System.err, new String[0]);
    }

    class StaticParent {

        @Option
        public String STRING1 = "String1ParentDefault";

        @Option
        public String STRING2 = "String2ParentDefault";

        @Option(overridable = true)
        public String STRING3 = "String3ParentDefault";

        public void doSomething() {
            System.out.println(STRING3);
        }

    }

    class OverridePropagation extends StaticParent {
        @Option
        public String STRING3 = "String3Overriden";
    }

    @Test
    public void testOveriddenOptions() {
        final OverridePropagation overridden = new OverridePropagation();
        final CommandLineParser overrideClp = new CommandLineParser(overridden);

        overrideClp.parseOptions(System.err, new String[0]);

        final OverridePropagation props = (OverridePropagation) overrideClp.getCallerOptions();
        Assert.assertTrue(props.STRING3.equals("String3Overriden"));
        Assert.assertTrue(((StaticParent) props).STRING3.equals("String3Overriden"));

        overrideClp.parseOptions(System.err, new String[]{"STRING3=String3Supplied"});

        final OverridePropagation propsSet = (OverridePropagation) overrideClp.getCallerOptions();
        Assert.assertTrue(propsSet.STRING3.equals("String3Supplied"));
        Assert.assertTrue(((StaticParent) propsSet).STRING3.equals("String3Supplied"));
    }


    @DataProvider(name = "testHtmlEscapeData")
    public Object[][] testHtmlEscapeData() {
        final List<Object[]> retval = new ArrayList<>();

        retval.add(new Object[]{"<", "&lt;"});
        retval.add(new Object[]{"x<y", "x&lt;y"});
        retval.add(new Object[]{"x<y<z", "x&lt;y&lt;z"});
        retval.add(new Object[]{"\n", "<p>"});
        retval.add(new Object[]{"<html> x<y </html> y< <strong> x </strong>","<html> x&lt;y </html> y&lt; <strong> x </strong>"});

        return retval.toArray(new Object[0][]);
    }

//    @Test(dataProvider = "testHtmlEscapeData")
//    public void testHtmlEscape(final String text, final String expected) {
//        Assert.assertEquals(CommandLineParser.htmlEscape(text), expected, "problems");
//    }

    @Test(dataProvider = "testHtmlEscapeData")
    public void testHtmlUnescape(final String expected, final String html) {
        Assert.assertEquals(CommandLineParser.htmlUnescape(html), expected, "problems");
    }

    @DataProvider(name = "testHTMLConverter")
    public Object[][] testHTMLConverterData() {
        final List<Object[]> retval = new ArrayList<>();

        retval.add(new Object[]{"hello", "hello"});
        retval.add(new Object[]{"", ""});
        retval.add(new Object[]{"hi</th>bye", "hi\tbye"});
        retval.add(new Object[]{"hi<th>bye", "hibye"});
        retval.add(new Object[]{"hi<li>bye", "hi - bye"});
        retval.add(new Object[]{"hi<NOT_A_REAL_TAG>bye", "hibye"});
        retval.add(new Object[]{"</h4><pre>", "\n\n"});
        retval.add(new Object[]{"<a href=\"http://go.here.org\"> string</ a >", " string (http://go.here.org)"});
        retval.add(new Object[]{"<a href=\"http://go.here.org\" > string</ a>", " string (http://go.here.org)"});
        retval.add(new Object[]{"< a href=\"http://go.here.org\"> string<a />", " string (http://go.here.org)"});


        //for some reason, the next test seems to break intelliJ, but it works on the commandline
        retval.add(new Object[]{"hi</li>bye", "hi\nbye"});

        retval.add(new Object[]{"Using read outputs from high throughput sequencing (HTS) technologies, this tool provides " +
                "metrics regarding the quality of read alignments to a reference sequence, as well as the proportion of the reads " +
                "that passed machine signal-to-noise threshold quality filters (Illumina)." +
                "<h4>Usage example:</h4>" +
                "<pre>" +
                "    java -jar picard.jar CollectAlignmentSummaryMetrics \\<br />" +
                "          R=reference_sequence.fasta \\<br />" +
                "          I=input.bam \\<br />" +
                "          O=output.txt" +
                "</pre>" +
                "Please see <a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#AlignmentSummaryMetrics'>" +
                "the AlignmentSummaryMetrics documentation</a> for detailed explanations of each metric. <br /> <br />" +
                "Additional information about Illumina's quality filters can be found in the following documents on the Illumina website:" +
                "<ul><li><a href=\"http://support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/hiseq-x-percent-pf-technical-note-770-2014-043.pdf\">" +
                "hiseq-x-percent-pf-technical-note</a></li>" +
                "<li><a href=\"http://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/hiseqx/hiseq-x-system-guide-15050091-d.pdf\">" +
                "hiseq-x-system-guide</a></li></ul>" +
                "<hr />",

                "Using read outputs from high throughput sequencing (HTS) technologies, this tool provides " +
                        "metrics regarding the quality of read alignments to a reference sequence, as well as the proportion of the reads " +
                        "that passed machine signal-to-noise threshold quality filters (Illumina)." +
                        "\nUsage example:\n" +
                        "\n" +
                        "    java -jar picard.jar CollectAlignmentSummaryMetrics \\\n" +
                        "          R=reference_sequence.fasta \\\n" +
                        "          I=input.bam \\\n" +
                        "          O=output.txt" +
                        "\n" +
                        "Please see the AlignmentSummaryMetrics documentation (http://broadinstitute.github.io/picard/picard-metric-definitions.html#AlignmentSummaryMetrics) for detailed explanations of each metric. \n \n" +
                        "Additional information about Illumina's quality filters can be found in the following documents on the Illumina website:" +
                        "\n" +
                        " - hiseq-x-percent-pf-technical-note (http://support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/hiseq-x-percent-pf-technical-note-770-2014-043.pdf)\n" +
                        " - hiseq-x-system-guide (http://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/hiseqx/hiseq-x-system-guide-15050091-d.pdf)\n\n" +
                        "\n"});

        return retval.toArray(new Object[0][]);
    }

    @Test(dataProvider = "testHTMLConverter")
    public void testHTMLConverter(String input, String expected) {
        final String converted = CommandLineParser.convertFromHtml(input);
        Assert.assertEquals(converted, expected, "common part:\"" + expected.substring(0, lengthOfCommonSubstring(converted, expected)) + "\"\n\n");
    }

    @CommandLineProgramProperties(
            usage = TestParserFail.USAGE_SUMMARY + TestParserFail.USAGE_DETAILS,
            usageShort = TestParserFail.USAGE_SUMMARY,
            programGroup = Testing.class
    )
    protected class TestParserFail extends CommandLineProgram {

        static public final String USAGE_DETAILS = "blah &blah; blah ";
        static public final String USAGE_SUMMARY = "This tool offers.....";

        @Override
        protected int doWork() {return 0;}
    }

    @CommandLineProgramProperties(
            usage = TestParserSucceed.USAGE_SUMMARY + TestParserSucceed.USAGE_DETAILS,
            usageShort = TestParserSucceed.USAGE_SUMMARY,
            programGroup = Testing.class
    )

    protected class TestParserSucceed extends CommandLineProgram {

        static public final String USAGE_DETAILS = "This is the first row <p>And this is the second";
        static public final String USAGE_SUMMARY = " X &lt; Y ";

        @Override
        protected int doWork() {return 0;}
    }

    @Test(expectedExceptions = AssertionError.class)
    public void testNonAsciiAssertion() {
        TestParserFail testparserFail = new TestParserFail();
        testparserFail.parseArgs(new String[]{});

        PrintStream stream = new PrintStream(new NullOutputStream());
        testparserFail.getCommandLineParser().usage(stream, true);
    }

    @Test
    public void testNonAsciiConverted() {
        TestParserSucceed testparserSucceed = new TestParserSucceed();
        testparserSucceed.parseArgs(new String[]{});

        ByteArrayOutputStream byteArrayOutputStream = new ByteArrayOutputStream();
        PrintStream stream = new PrintStream(byteArrayOutputStream);
        testparserSucceed.getCommandLineParser().usage(stream, true);

        String expected = "USAGE: TestParserSucceed [options]\n" +
                "\n" +
                "Documentation: http://broadinstitute.github.io/picard/command-line-overview.html#TestParserSucceed\n" +
                "\n" +
                " X < Y This is the first row \n" +
                "And this is the second";
        Assert.assertEquals(byteArrayOutputStream.toString().substring(0, expected.length()), expected);
    }

    static private int lengthOfCommonSubstring(String lhs, String rhs) {
        int i = 0;
        while (i < Math.min(lhs.length(), rhs.length()) && lhs.charAt(i) == rhs.charAt(i)) i++;

        return i;
    }

    private class NullOutputStream extends OutputStream {
        @Override
        public void write(final int b) throws IOException {

        }
    }
}
