package picard.cmdline;

import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.testng.Assert;
import org.testng.annotations.Test;
import java.util.*;

/**
 * Created by farjoun on 9/10/15.
 */
public class PicardCommandLineTest {

    @Test
    public void TestPicardPublic() { // this tests fails if any CLP in picard is missing its @CommandLineProgramProperties annotation
        PicardCommandLine picardCommandLine = new PicardCommandLine();
        picardCommandLine.instanceMain(new String[]{""});
    }

    @Test
    public void testProcessCommandLinePrograms() {
        final List<Class<CommandLineProgram>> allCLPs = new ArrayList<>();

        PicardCommandLine.processAllCommandLinePrograms(
                Collections.singletonList("picard"),
                (Class<CommandLineProgram> clazz, CommandLineProgramProperties clProperties) -> {
                    allCLPs.add(clazz);
                });
        Assert.assertTrue(allCLPs.size() > 0);
        allCLPs.forEach(clazz -> Assert.assertTrue(CommandLineProgram.class.isAssignableFrom(clazz)));
    }

    // Find every CommandLineProgram in picard, instantiate it, and initialize a Barclay parser using the
    // resulting object. This catches any CommandLineProgram that depends on using duplicate (of overrideable) arguments,
    // which are not accepted by the Barclay parser, and have no other tests where this would surface.
    //
    // NOTE that it does NOT test that those tools actually work correctly, only that they aren't immediately
    // rejected with a CommandLineParserInternalException from the Barclay parser.
    @Test
    public void testLaunchAllCommandLineProgramsWithBarclayParser() {
        PicardCommandLine.processAllCommandLinePrograms(
                Collections.singletonList("picard"),
                (Class<CommandLineProgram> clazz, CommandLineProgramProperties clProperties) -> {
                    // Check for missing annotations
                    if (null != clProperties) {
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
                    }
                }
        );
    }

    @Test
    public void testPrintUsage() {
        Assert.assertEquals(new PicardCommandLine().instanceMain(new String[]{"-h"}), 1);
    }
}