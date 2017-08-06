package picard.cmdline;

import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.testng.annotations.Test;

import java.lang.reflect.Modifier;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import static picard.cmdline.PicardCommandLine.getProgramProperty;

/**
 * Created by farjoun on 9/10/15.
 */
public class PicardCommandLineTest {

    @Test
    public void TestPicardPublic() { // this tests fails if any CLP in picard is missing its @CommandLineProgramProperties annotation
        PicardCommandLine picardCommandLine = new PicardCommandLine();
        picardCommandLine.instanceMain(new String[]{""});
    }

    // Find every CommandLineProgram in picard, instantiate it, and initialize a Barclay parser using the
    // resulting object. This catches any CommandLineProgram that depends on using duplicate (of overridable) arguments,
    // which are not accepted by the Barclay parser, and have no other tests where this would surface.
    //
    // NOTE that it does NOT test that those tools actually work correctly, only that they aren't immediately
    // rejected with a CommandLineParserInternalException from the Barclay parser.
    @Test
    public void testLaunchAllCommandLineProgramsWithBarclayParser() {
        final ClassFinder classFinder = new ClassFinder();
        for (final String pkg : Collections.singletonList("picard")) {
            classFinder.find(pkg, CommandLineProgram.class);
        }

        for (final Class clazz : classFinder.getClasses()) {
            // No interfaces, synthetic, primitive, local, or abstract classes.
            if (!clazz.isInterface() && !clazz.isSynthetic() && !clazz.isPrimitive() && !clazz.isLocalClass()
                    && !Modifier.isAbstract(clazz.getModifiers())) {
                final CommandLineProgramProperties clProperties = getProgramProperty(clazz);
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
        }
    }

}