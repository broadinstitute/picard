package picard.cmdline;

import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
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

}