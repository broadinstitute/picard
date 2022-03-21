package picard.cmdline;

import htsjdk.samtools.util.Log;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Class for handling translation of Picard-style command line argument syntax to POSIX-style argument syntax;
 * used for running tests written with Picard style syntax against the Barclay command line parser.
 */
public class CommandLineSyntaxTranslater {
    private final static Log log = Log.getInstance(CommandLineSyntaxTranslater.class);

    // Prefixes used by the Barclay parser for short/long prefixes
    private static final String BARCLAY_SHORT_OPTION_PREFIX = "-";
    private static final String BARCLAY_LONG_OPTION_PREFIX = "--";
    // Separator used by the legacy parser for name/value arguments
    private static final String LEGACY_VALUE_SEPARATOR = "=";

    // Return true when the command line arguments appear to use Picard's legacy syntax.
    public static boolean isLegacyPicardStyle(final String argv[]) {
        final boolean anyLegacy = Arrays.stream(argv).anyMatch(
                arg -> !arg.startsWith(BARCLAY_SHORT_OPTION_PREFIX) &&
                        !arg.startsWith(BARCLAY_LONG_OPTION_PREFIX) &&
                        arg.contains(LEGACY_VALUE_SEPARATOR)
        );
        if (anyLegacy && Arrays.stream(argv).anyMatch(
                arg -> arg.startsWith(BARCLAY_SHORT_OPTION_PREFIX) || arg.startsWith(BARCLAY_LONG_OPTION_PREFIX))) {
            // There appear to be both legacy and posix style args. Prefer/choose posix in this case since there are
            // legitimate cases where argument values might contain embedded "=" (i.e,
            // "--INPUT path/to/some.bam --SOME_ARG date=01/01/2022"), which makes them appear to be
            // legacy style args, even though they are not), whereas its very unlikely to encounter a legitimate
            // legacy option that starts with a posix prefix ("--" or "-")
            log.warn("!!!!!!Possible mixed (legacy and new style) arguments detected!!!!!!!\n"
                    + "Assuming new-style arguments are intended. See: " + CommandLineProgram.SYNTAX_TRANSITION_URL);
            return false;
        }
        return anyLegacy;
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
