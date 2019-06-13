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

import com.intel.gkl.compression.IntelDeflaterFactory;
import com.intel.gkl.compression.IntelInflaterFactory;
import htsjdk.samtools.Defaults;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMFileWriterImpl;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.metrics.StringHeader;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.BlockGunzipper;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.barclay.argparser.CommandLineParserOptions;
import org.broadinstitute.barclay.argparser.LegacyCommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.SpecialArgumentsCollection;
import picard.cmdline.argumentcollections.OptionalReferenceArgumentCollection;
import picard.cmdline.argumentcollections.ReferenceArgumentCollection;
import picard.cmdline.argumentcollections.RequiredReferenceArgumentCollection;
import picard.nio.PathProvider;
import picard.util.PropertyUtils;

import java.io.File;
import java.net.InetAddress;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Properties;
import java.util.stream.Collectors;

/**
 * Abstract class to facilitate writing command-line programs.
 *
 * To use:
 *
 * 1. Extend this class with a concrete class that is annotated with @COmmandLineProgramProperties, and has data members
 * annotated with @Argument, @PositionalArguments, and/or @ArgumentCollection annotations.
 *
 * 2. If there is any custom command-line validation, override customCommandLineValidation().  When this method is
 * called, the command line has been parsed and set into the data members of the concrete class.
 *
 * 3. Implement a method doWork().  This is called after successful command-line processing.  The value it returns is
 * the exit status of the program.  It is assumed that the concrete class emits any appropriate error message before
 * returning non-zero.  doWork() may throw unchecked exceptions, which are caught and reported appropriately.
 *
 */
public abstract class CommandLineProgram {
    // Picard CmdLine properties file resource path, placed at this path by the gradle build script
    private static String PICARD_CMDLINE_PROPERTIES_FILE = "picard/picardCmdLine.properties";
    private static String PROPERTY_USE_LEGACY_PARSER = "picard.useLegacyParser";
    private static String PROPERTY_CONVERT_LEGACY_COMMAND_LINE = "picard.convertCommandLine";
    private static Boolean useLegacyParser;

    /**
     * CommandLineProgramProperties oneLineSummary attribute must be shorted than this in order to maintain
     * reasonable help output formatting.
     */
    public static int MAX_ALLOWABLE_ONE_LINE_SUMMARY_LENGTH = 120;

    @Argument(doc="One or more directories with space available to be used by this program for temporary storage of working files",
            common=true, optional=true)
    public List<File> TMP_DIR = new ArrayList<>();

    @Argument(doc = "Control verbosity of logging.", common=true)
    public Log.LogLevel VERBOSITY = Log.LogLevel.INFO;

    @Argument(doc = "Whether to suppress job-summary info on System.err.", common=true)
    public Boolean QUIET = false;

    @Argument(doc = "Validation stringency for all SAM files read by this program.  Setting stringency to SILENT " +
            "can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) " +
            "do not otherwise need to be decoded.", common=true)
    public ValidationStringency VALIDATION_STRINGENCY = ValidationStringency.DEFAULT_STRINGENCY;

    @Argument(doc = "Compression level for all compressed files created (e.g. BAM and VCF).", common=true)
    public int COMPRESSION_LEVEL = Defaults.COMPRESSION_LEVEL;

    @Argument(doc = "When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. " +
            "Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.", optional=true, common=true)
    public Integer MAX_RECORDS_IN_RAM = SAMFileWriterImpl.getDefaultMaxRecordsInRam();

    @Argument(doc = "Whether to create a BAM index when writing a coordinate-sorted BAM file.", common=true)
    public Boolean CREATE_INDEX = Defaults.CREATE_INDEX;

    @Argument(doc="Whether to create an MD5 digest for any BAM or FASTQ files created.  ", common=true)
    public boolean CREATE_MD5_FILE = Defaults.CREATE_MD5;

    @ArgumentCollection
    public ReferenceArgumentCollection referenceSequence = makeReferenceArgumentCollection();

    // This is retained for compatibility with existing code that depends on accessing it, and is populated
    // after argument parsing using the value established by the user in the referenceSequence argument collection.
    protected File REFERENCE_SEQUENCE = Defaults.REFERENCE_FASTA;

    @Argument(doc="Google Genomics API client_secrets.json file path.", common = true)
    public String GA4GH_CLIENT_SECRETS="client_secrets.json";

    @ArgumentCollection(doc="Special Arguments that have meaning to the argument parsing system.  " +
                "It is unlikely these will ever need to be accessed by the command line program")
    public Object specialArgumentsCollection = useLegacyParser(getClass()) ?
            new Object() : // legacy parser does not require these
            new SpecialArgumentsCollection();

    @Argument(shortName = "use_jdk_deflater", doc = "Use the JDK Deflater instead of the Intel Deflater for writing compressed output", common = true)
    public Boolean USE_JDK_DEFLATER = false;

    @Argument(shortName = "use_jdk_inflater", doc = "Use the JDK Inflater instead of the Intel Inflater for reading compressed input", common = true)
    public Boolean USE_JDK_INFLATER = false;

    private static final String[] PACKAGES_WITH_WEB_DOCUMENTATION = {"picard"};

    static {
      // Register custom reader factory for reading data from Google Genomics
      // implementation of GA4GH API.
      // With this it will be possible to pass these urls as INPUT params.
      // E.g. java -jar dist/picard.jar ViewSam \
      //    INPUT=https://www.googleapis.com/genomics/v1beta2/readgroupsets/CK256frpGBD44IWHwLP22R4/ \
      //    GA4GH_CLIENT_SECRETS=../client_secrets.json
      if (System.getProperty("samjdk.custom_reader") == null) {
        System.setProperty("samjdk.custom_reader",
            "https://www.googleapis.com/genomics," +
            "com.google.cloud.genomics.gatk.htsjdk.GA4GHReaderFactory");
      }
    }

    /**
    * Initialized in parseArgs.  Subclasses may want to access this to do their
    * own validation, and then print usage using commandLineParser.
    */
    private CommandLineParser commandLineParser;

    private final List<Header> defaultHeaders = new ArrayList<>();

    /**
    * The reconstructed commandline used to run this program. Used for logging
    * and debugging.
    */
    private String commandLine;

    /**
    * Do the work after command line has been parsed. RuntimeException may be
    * thrown by this method, and are reported appropriately.
    * @return program exit status.
    */
    protected abstract int doWork();

    protected boolean requiresReference() {
        return false;
    }

    protected ReferenceArgumentCollection makeReferenceArgumentCollection() {
        return requiresReference() ?
                new RequiredReferenceArgumentCollection() :
                new OptionalReferenceArgumentCollection();
    }

    public void instanceMainWithExit(final String[] argv) {
        System.exit(instanceMain(argv));
    }

    public int instanceMain(final String[] argv) {
        String actualArgs[] = argv;

        if (System.getProperty(PROPERTY_CONVERT_LEGACY_COMMAND_LINE, "false").equals("true")) {
            actualArgs = CommandLineSyntaxTranslater.convertPicardStyleToPosixStyle(argv);
        } else if (CommandLineSyntaxTranslater.isLegacyPicardStyle(argv)) {
            final String[] messageLines = new String[] {
                "", "",
                "********** NOTE: Picard's command line syntax is changing.",
                "**********",
                "********** For more information, please see:",
                "********** https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)",
                "**********",
                "********** The command line looks like this in the new syntax:",
                "**********",
                "**********    %s %s",
                "**********",
                "", ""
            };
            final String message = String.join("\n", messageLines);
            final String syntax  = String.join(" ", CommandLineSyntaxTranslater.convertPicardStyleToPosixStyle(argv));
            final String info    = String.format(message, this.getClass().getSimpleName(), syntax);
            Log.getInstance(this.getClass()).info(info);
        }
        if (!parseArgs(actualArgs)) {
            return 1;
        }

        // Provide one temp directory if the caller didn't
        if (this.TMP_DIR == null) this.TMP_DIR = new ArrayList<>();
        if (this.TMP_DIR.isEmpty()) TMP_DIR.add(IOUtil.getDefaultTmpDir());

        // Build the default headers
        final Date startDate = new Date();
        this.defaultHeaders.add(new StringHeader(commandLine));
        this.defaultHeaders.add(new StringHeader("Started on: " + startDate));

        Log.setGlobalLogLevel(VERBOSITY);
        if (System.getProperty("ga4gh.client_secrets") == null) {
          System.setProperty("ga4gh.client_secrets", GA4GH_CLIENT_SECRETS);
        }
        SamReaderFactory.setDefaultValidationStringency(VALIDATION_STRINGENCY);
        BlockCompressedOutputStream.setDefaultCompressionLevel(COMPRESSION_LEVEL);

        if (VALIDATION_STRINGENCY != ValidationStringency.STRICT) VariantContextWriterBuilder.setDefaultOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);

        if (MAX_RECORDS_IN_RAM != null) {
            SAMFileWriterImpl.setDefaultMaxRecordsInRam(MAX_RECORDS_IN_RAM);
        }

        if (CREATE_INDEX) {
            SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(true);
        }

        SAMFileWriterFactory.setDefaultCreateMd5File(CREATE_MD5_FILE);

        for (final File f : TMP_DIR) {
            // Intentionally not checking the return values, because it may be that the program does not
            // need a tmp_dir. If this fails, the problem will be discovered downstream.
            if (!f.exists()) f.mkdirs();
            f.setReadable(true, false);
            f.setWritable(true, false);
            System.setProperty("java.io.tmpdir", f.getAbsolutePath()); // in loop so that last one takes effect
        }

        if (!USE_JDK_DEFLATER) {
            BlockCompressedOutputStream.setDefaultDeflaterFactory(new IntelDeflaterFactory());
        }

        if (!USE_JDK_INFLATER) {
            BlockGunzipper.setDefaultInflaterFactory(new IntelInflaterFactory());
        }

        if (!QUIET) {
            System.err.println("[" + new Date() + "] " + commandLine);

            // Output a one liner about who/where and what software/os we're running on
            try {
                final String pathProvidersMessage =
                        Arrays.stream(PathProvider.values())
                                .map(provider -> String.format("Provider %s is%s available;", provider.name(), provider.isAvailable ? "" : " not"))
                                .collect(Collectors.joining(" "));

                final boolean usingIntelDeflater = (BlockCompressedOutputStream.getDefaultDeflaterFactory() instanceof IntelDeflaterFactory &&
                        ((IntelDeflaterFactory)BlockCompressedOutputStream.getDefaultDeflaterFactory()).usingIntelDeflater());
                final boolean usingIntelInflater = (BlockGunzipper.getDefaultInflaterFactory() instanceof IntelInflaterFactory &&
                        ((IntelInflaterFactory)BlockGunzipper.getDefaultInflaterFactory()).usingIntelInflater());
                final String msg = String.format(
                    "[%s] Executing as %s@%s on %s %s %s; %s %s; Deflater: %s; Inflater: %s; %s Picard version: %s",
                    new Date(), System.getProperty("user.name"), InetAddress.getLocalHost().getHostName(),
                    System.getProperty("os.name"), System.getProperty("os.version"), System.getProperty("os.arch"),
                    System.getProperty("java.vm.name"), System.getProperty("java.runtime.version"),
                    usingIntelDeflater ? "Intel" : "Jdk", usingIntelInflater ? "Intel" : "Jdk",
                        pathProvidersMessage,
                    getCommandLineParser().getVersion());
                System.err.println(msg);
            }
            catch (Exception e) { /* Unpossible! */ }
        }

        int ret = -1;
        try {
            ret = doWork();
        } finally {
            try {
                // Emit the time even if program throws
                if (!QUIET) {
                    final Date endDate = new Date();
                    final double elapsedMinutes = (endDate.getTime() - startDate.getTime()) / (1000d * 60d);
                    final String elapsedString  = new DecimalFormat("#,##0.00").format(elapsedMinutes);
                    System.err.println("[" + endDate + "] " + getClass().getName() + " done. Elapsed time: " + elapsedString + " minutes.");
                    System.err.println("Runtime.totalMemory()=" + Runtime.getRuntime().totalMemory());
                    if (ret != 0 && hasWebDocumentation(this.getClass())) System.err.println(getFaqLink());
                }
            }
            catch (Throwable e) {
                // do nothing
            }
        }
        return ret;
    }

    /**
    * Put any custom command-line validation in an override of this method.
    * clp is initialized at this point and can be used to print usage and access argv.
     * Any options set by command-line parser can be validated.
    * @return null if command line is valid.  If command line is invalid, returns an array of error message
    * to be written to the appropriate place.
    */
    protected String[] customCommandLineValidation() {
        return null;
    }

    /**
    *
    * @return true if command line is valid
    */
    protected boolean parseArgs(final String[] argv) {

        commandLineParser = getCommandLineParser();

        boolean ret;
        try {
            ret = commandLineParser.parseArguments(System.err, argv);
        } catch (CommandLineException e) {
            // Barclay command line parser throws on parsing/argument errors
            System.err.println(commandLineParser.usage(false,false));
            System.err.println(e.getMessage());
            ret = false;
        }

        commandLine = commandLineParser.getCommandLine();
        if (!ret) {
            return false;
        }
        REFERENCE_SEQUENCE = referenceSequence.getReferenceFile();

        final String[] customErrorMessages = customCommandLineValidation();
        if (customErrorMessages != null) {
            System.err.print(commandLineParser.usage(false, false));
            for (final String msg : customErrorMessages) {
                System.err.println(msg);
            }
            return false;
        }
        return true;
    }

    /** Gets a MetricsFile with default headers already written into it. */
    protected <A extends MetricBase,B extends Comparable<?>> MetricsFile<A,B> getMetricsFile() {
        final MetricsFile<A,B> file = new MetricsFile<>();
        for (final Header h : this.defaultHeaders) {
            file.addHeader(h);
        }

        return file;
    }

    public String getStandardUsagePreamble() {
        return getCommandLineParser().getStandardUsagePreamble(getClass());
    }

    /**
     * @return Return the command line parser to be used.
     */
    public CommandLineParser getCommandLineParser() {
        if (commandLineParser == null) {
            commandLineParser = useLegacyParser(getClass()) ?
                        new LegacyCommandLineArgumentParser(this) :
                        new CommandLineArgumentParser(this,
                            Collections.EMPTY_LIST,
                            new HashSet<>(Collections.singleton(CommandLineParserOptions.APPEND_TO_COLLECTIONS)));
        }
        return commandLineParser;
    }

    /**
     * Return true if the Picard command line parser should be used in place of the Barclay command line parser,
     * otherwise, false.
     *
     * The Barclay parser is enabled by opt-in only, via the presence of a (true-valued) boolean property
     * "picard.useLegacyParser", either as a System property, or property in a file called "picard.properties".
     * System property value takes precedence.
     *
     * @return true if the legacy parser should be used
     */
    public static boolean useLegacyParser(final Class<?> clazz) {
        if (useLegacyParser == null) {
            String legacyPropertyValue = legacyPropertyValue = System.getProperty(PROPERTY_USE_LEGACY_PARSER);
            if (null == legacyPropertyValue){
                final Properties props = PropertyUtils.loadPropertiesFile(PICARD_CMDLINE_PROPERTIES_FILE, clazz);
                if (props != null) {
                    legacyPropertyValue = props.getProperty(PROPERTY_USE_LEGACY_PARSER);
                }
            }
            // remember the value so we only have to load the properties file once
            useLegacyParser = legacyPropertyValue == null ? true : Boolean.parseBoolean(legacyPropertyValue);
        }
        return useLegacyParser;
    }

    /**
     * @return Version stored in the manifest of the jarfile.
     */
    public String getVersion() {
        return getCommandLineParser().getVersion();
    }

    public String getCommandLine() {
        return commandLine;
    }

    public void setDefaultHeaders(final List<Header> headers) {
        this.defaultHeaders.clear();
        this.defaultHeaders.addAll(headers);
    }

    public List<Header> getDefaultHeaders() {
        return this.defaultHeaders;
    }

    /**
     * A typical command line program will call this to get the beginning of the usage message,
     * and then typically append a description of the program, like this:
     *
     * public String USAGE = getStandardUsagePreamble(getClass()) + "Frobnicate the freebozle."
     */
    public static String getStandardUsagePreamble(final Class<?> mainClass) {
        return "USAGE: " + mainClass.getSimpleName() +" [options]\n\n" +
                (hasWebDocumentation(mainClass) ?
                        "Documentation: http://broadinstitute.github.io/picard/command-line-overview.html" +
                                mainClass.getSimpleName() + "\n\n" :
                        "");
    }

    /**
     * Determine if a class has web documentation based on its package name
     *
     * @param clazz
     * @return true if the class has web documentation
     */
    public static boolean hasWebDocumentation(final Class<?> clazz){
        for (final String pkg: PACKAGES_WITH_WEB_DOCUMENTATION) {
            if (clazz.getPackage().getName().startsWith(pkg)) {
                return true;
            }
        }
        return false;
    }

    /**
     * @return the link to a FAQ
     */
    public static String getFaqLink() {
        return "To get help, see http://broadinstitute.github.io/picard/index.html#GettingHelp";
    }
}
