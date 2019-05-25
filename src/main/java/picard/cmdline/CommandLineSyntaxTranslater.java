package picard.cmdline;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Class for handling translation of Picard-style command line argument syntax to POSIX-style argument syntax;
 * used for running tests written with Picard style syntax against the Barclay command line parser.
 */
public class CommandLineSyntaxTranslater {

    // Prefixes used by the Barclay parser for short/long prefixes
    private static final String BARCLAY_SHORT_OPTION_PREFIX = "-";
    private static final String BARCLAY_LONG_OPTION_PREFIX = "--";
    // Separator used by the legacy parser for name/value arguments
    private static final String LEGACY_VALUE_SEPARATOR = "=";

    // Return true when the command line arguments appear to use Picard's legacy syntax.
    public static boolean isLegacyPicardStyle(final String argv[]) {
        return Arrays.stream(argv).anyMatch(
                putativeLegacyArg ->
                        !putativeLegacyArg.startsWith(BARCLAY_SHORT_OPTION_PREFIX) &&
                           !putativeLegacyArg.startsWith(BARCLAY_LONG_OPTION_PREFIX) &&
                                putativeLegacyArg.contains(LEGACY_VALUE_SEPARATOR)
        );
    }

    public static String[] convertPicardStyleToPosixStyle(final String argv[]) {
        final List<String> convertedArgs = Arrays.stream(argv).flatMap(
            originalArgPair -> {
                final String[] splitArgPair = originalArgPair.split(LEGACY_VALUE_SEPARATOR, 2);
                if (splitArgPair.length == 1) {   // assume positional arg
                    return Arrays.stream(new String[]{ originalArgPair });
                } else if (splitArgPair.length == 2) {
                    //deal with EXTRA_ARGUMENT in CollectMultipleMetrics
                    if (splitArgPair[0].equals("EXTRA_ARGUMENT")) {
                        splitArgPair[1] = splitArgPair[1]
                                .replace("::", "::" + BARCLAY_LONG_OPTION_PREFIX)
                                .replace(LEGACY_VALUE_SEPARATOR, " ");
                    }
                    // it doesn't matter whether we use the short short name token ("-") or the long name token
                    // ("--"), so just treat everything as if it were a short name, since the CLP will accept either
                    return Arrays.stream(new String[]{BARCLAY_SHORT_OPTION_PREFIX + splitArgPair[0], splitArgPair[1]});
                }
                else {
                    throw new RuntimeException("Cannot convert this argument: " + originalArgPair);
                }
            }
        ).collect(Collectors.toList());
        return convertedArgs.toArray(new String[convertedArgs.size()]);
    }
}
