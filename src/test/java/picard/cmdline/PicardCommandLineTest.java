package picard.cmdline;

import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import java.util.*;

public class PicardCommandLineTest {

    // cache all the command line programs once so we don't have to rediscover them for each test
    private static Map<Class<CommandLineProgram>, CommandLineProgramProperties> allCLPS = new HashMap<>();

    @BeforeClass
    public static void getAllCommandLineProgramClasses() {
        // cache all the command line programs once so we don't have to rediscover them for each test
        PicardCommandLine.processAllCommandLinePrograms(
                PicardCommandLine.getPackageList(),
                (Class<CommandLineProgram> clazz, CommandLineProgramProperties clProperties) -> {
                    allCLPS.put(clazz, clProperties);
                });
        Assert.assertTrue(allCLPS.size() > 0);
    }

    @Test
    public void TestPicardPublic() { // this tests fails if any CLP in picard is missing its @CommandLineProgramProperties annotation
        PicardCommandLine picardCommandLine = new PicardCommandLine();
        picardCommandLine.instanceMain(new String[]{""});
    }

    @Test
    public void testProcessCommandLinePrograms() {
        allCLPS.forEach((clazz, properties) -> {
           Assert.assertTrue(CommandLineProgram.class.isAssignableFrom(clazz));
        });
    }

    // Find every CommandLineProgram in picard, instantiate it, and initialize a Barclay parser using the
    // resulting object. This catches any CommandLineProgram that depends on using duplicate (of overrideable) arguments,
    // which are not accepted by the Barclay parser, and have no other tests where this would surface.
    //
    // NOTE that it does NOT test that those tools actually work correctly, only that they aren't immediately
    // rejected with a CommandLineParserInternalException from the Barclay parser.
    @Test
    public void testLaunchAllCommandLineProgramsWithBarclayParser() {
        allCLPS.forEach((Class<CommandLineProgram> clazz, CommandLineProgramProperties clProperties) -> {
            // Check for missing annotations
            Assert.assertNotNull(clProperties);
            try {
                final Object commandLineProgram = clazz.newInstance();
                try {
                    new CommandLineArgumentParser(commandLineProgram);
                } catch (CommandLineException.CommandLineParserInternalException e) {
                    throw new RuntimeException("Barclay command line parser internal exception parsing class: " + clazz.getName(), e);
                }
            } catch (IllegalAccessException | InstantiationException e) {
                throw new RuntimeException("Failure instantiating command line program: " + clazz.getName(), e);
            }
        });
    }

    @Test
    public void testPrintUsage() {
        Assert.assertEquals(new PicardCommandLine().instanceMain(new String[]{"-h"}), 1);
    }

    @Test
    public void testCommandlineProgramPropertiesOneLineSummaryLength() {
        // Find and test each command line tool to ensure that one line
        // summaries don't exceed the maximum allowable length
        allCLPS.forEach(
            (Class<picard.cmdline.CommandLineProgram> clazz, CommandLineProgramProperties clProperties) -> {
                    Assert.assertNotNull(clProperties);
                    Assert.assertTrue(
                        clProperties.oneLineSummary().length() <= CommandLineProgram.MAX_ALLOWABLE_ONE_LINE_SUMMARY_LENGTH,
                        String.format("One line summary for tool '%s' exceeds allowable length of %d",
                                clazz.getCanonicalName(),
                                CommandLineProgram.MAX_ALLOWABLE_ONE_LINE_SUMMARY_LENGTH));
            });
    }

    @DataProvider(name="isLegacyPicardStyleTests")
    public final Object[][] getIsLegacyPicardStyle() {
        return new Object[][] {
                //arg list, is legacy style

                // legacy base cases
                {Arrays.asList("--INPUT", "path/to/some.bam"), false},
                {Arrays.asList("--INPUT", "path/to/some.bam", "--VALIDATION_STRINGENCY", "LENIENT"), false},

                // posix base cases
                {Arrays.asList("INPUT=path/to/some.bam"), true },
                {Arrays.asList("INPUT=path/to/some.bam", "VALIDATION_STRINGENCY=LENIENT"), true},

                // mixed syntax cases

                // APPEARS to isLegacyPicardStyle to contain a mix of styles, but is actually (in theory) a legitimate
                // posix style arg list, so select the posix parser, but issue a warning about possible mixed
                // args set
                {Arrays.asList("--INPUT", "path/to/some.bam", "--SOME_ARG", "date=01/01/2022"), false},

                // appears to isLegacyPicardStyle to contain a mix of styles, but is probably not valid, so select the
                // posix parser, issue a warning, and let the parser decide if its legitimate
                {Arrays.asList("--INPUT", "path/to/some.bam", "VALIDATION_STRINGENCY=LENIENT"), false},

                // appears to isLegacyPicardStyle to contain a mix of styles, but is probably not valid, so select the
                // posix parser, issue a warning, and let the parser decide if its legitimate
                {Arrays.asList("INPUT=path/to/some.bam", "--ARG=somevalue"), false},
        };
    }

    @Test(dataProvider="isLegacyPicardStyleTests")
    public void testIsLegacyPicardStyle(final List<String> args, final boolean isLegacy) {
        Assert.assertEquals(CommandLineSyntaxTranslater.isLegacyPicardStyle(args.toArray(new String[0])), isLegacy);
    }
}