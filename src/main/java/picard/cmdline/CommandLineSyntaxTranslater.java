package picard.cmdline;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Class for handling translation of Picard-style command line argument syntax to POSIX-style argument syntax;
 * used for running tests written with Picard style syntax against the Barclay command line parser.
 */
public class CommandLineSyntaxTranslater {

    public static String[] translatePicardStyleToPosixStyle(final String argv[]) {
        final List<String> convertedArgs = Arrays.stream(argv).flatMap(
            originalArgPair -> {
                final String[] splitArgPair = originalArgPair.split("=", -1);
                if (splitArgPair.length == 1) {   // assume positional arg
                    return Arrays.stream(new String[]{ originalArgPair });
                } else if (splitArgPair.length == 2) {
                    // it doesn't matter whether we use the short short name token ("-") or the long name token
                    // ("--"), so just treat everything as if it were a short name, since the CLP will accept either
                    return Arrays.stream(new String[]{"-" + splitArgPair[0], splitArgPair[1]});
                }
                else {
                    throw new RuntimeException(
                            "Argument syntax conversion failed. Too many \"=\" separated tokens to translate: " + originalArgPair);
                }
            }
        ).collect(Collectors.toList());
        return convertedArgs.toArray(new String[convertedArgs.size()]);

    }

}
