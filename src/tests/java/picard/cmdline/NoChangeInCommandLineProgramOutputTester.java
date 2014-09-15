package picard.cmdline;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProcessExecutor;
import htsjdk.samtools.util.StringUtil;
import org.testng.Assert;
import picard.PicardException;
import picard.analysis.CollectAlignmentSummaryMetrics;
import picard.analysis.CollectWgsMetrics;
import sun.tools.jar.resources.jar;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * This Tester is for when you need to modify or refactor an existing CommandLineProgram (CLP) without changing its output relative to what
 * it was before. The Tester compares the CLP's output for a given set of input files and/or command line args to the output of the clean
 * latest-picard-release version of this CLP. It works by running the clean CLP .jar on the inputs, then running the current CLP code
 * (eg. from your development branch) on the same inputs, and then comparing the outputs in an intelligent way (eg. ignoring trivial
 * differences in the file headers, or unordered-field value order, while asserting no changes in anything else). It supports comparison of
 * various output formats (BAM, SAM, text file), and can be extended to support additional formats.
 *
 * The Tester is different from other Test classes such as {@link picard.sam.CompareSAMsTest} and {@link picard.util.TextFileParsersTest}
 * in that it works with existing input files rather than synthetic data (making it easier to set up some realistic tests prior to starting
 * development). Also, instead of comparing the test output to a pre-generated output file, its will execute the released CLP and generate
 * the output files for you. This makes it simple to re-test with additional input files or with different command line args without having
 * to manually generate the expected outputs for each test case.
 */
public abstract class NoChangeInCommandLineProgramOutputTester {

    private static final Log log = Log.getInstance(CollectAlignmentSummaryMetrics.class);

    public abstract File getReleasedCommandLineProgramJarPath();
    public abstract CommandLineProgram getProgram();
    public abstract InputFileArg[] getInputFileArgs();
    public abstract Arg[] getOtherArgs();
    public abstract OutputFileArg[] getOutputFileArgs();

    private final File outputDir;
    private final boolean deleteOnExit;

    /*
    Maps lower-case format string (eg. "bam") to a comparator for comparing text files. Format string is arbitrary (eg. "genotype.vcf").
    Comparators can be added or replaced by subclasses to add special comparisons for outputs of particular CLPs.
    */
    protected final Map<String, OutputFormatComparator> formatToDefaultComparatorMap = new HashMap<String, OutputFormatComparator>();

    public NoChangeInCommandLineProgramOutputTester() {
        this(true);
    }

    public NoChangeInCommandLineProgramOutputTester(boolean deleteOutputsOnExit) {
        this.deleteOnExit = deleteOutputsOnExit;
        this.outputDir = IOUtil.createTempDir(this.getClass().getName() + ".", ".tmp");
        if(deleteOnExit){
            outputDir.deleteOnExit();
        }

        formatToDefaultComparatorMap.put("txt", new TextOutputFileComparator());
    }

    /**
     * Sets up the basic command line arguments for input and output and runs instanceMain.
     */
    public void runTest() {
        //validate release CLP jar
        final File jarPath = getReleasedCommandLineProgramJarPath();
        if(jarPath == null)
            throw new IllegalStateException("Must specify jar file for the clean released version of this CLP");
        IOUtil.assertFileIsReadable(jarPath);

        //accumulate args
        final ArrayList<String> args = new ArrayList<String>();

        //validate and process input files
        long newestInputModTime = jarPath.lastModified(); //for finding when any of the CLP inputs were last modified (including the jar)
        for(final InputFileArg input : getInputFileArgs()) {
            IOUtil.assertFileIsReadable(input.inputFilePath);
            newestInputModTime = Math.max(newestInputModTime, input.inputFilePath.lastModified());
            args.add(input.key + "=" + input.inputFilePath);
        }

        //process regular args
        for(Arg arg : getOtherArgs()) {
            args.add(arg.key+"="+arg.value);
        }

        //use hash to make sure that if any of the inputs change, new expected output files will be generated.
        final String inputsHash = Integer.toString((jarPath + " " + StringUtil.join("   ", args)).hashCode());

        final ArrayList<String> releasedCLPArgs = new ArrayList<String>(args); //create different set of args for released CLP to allow for different output file names

        //validate and process output files and check whether they need to be regenerated.
        //if all output files exit and are newer than both the released jar file and all input files, then no need to run the released CLP.
        final OutputFileArg[] outputs = getOutputFileArgs();
        if(outputs == null || outputs.length == 0)
            throw new IllegalStateException("Must define at least one output file");

        boolean needToRegenerateExpectedOutputs = false;
        final List<OutputFileComparison> comparisonsToDoAfterRun = new LinkedList<OutputFileComparison>();
        for(final OutputFileArg output : outputs) {
            //check for unknown output file format here before running anything
            final OutputFormatComparator comparator;
            if(output.comparator != null) {
                comparator = output.comparator;
            } else {
                comparator = formatToDefaultComparatorMap.get(output.fileFormat);
                if(comparator == null)
                    throw new IllegalStateException("No comparator found for output format: "+output.fileFormat+" in formatToComparatorMap");
            }


            final String cleanedUpKey = output.key.replaceAll("[`~!@#$%^&*()+<>?\\[\\]|{},./]", "_"); //remove special chars

            final File expectedOutput = new File(outputDir, cleanedUpKey + ".expected_output."+inputsHash+"."+output.fileFormat);
            releasedCLPArgs.add(output.key+"="+expectedOutput);
            final File actualOutput = new File(outputDir, cleanedUpKey + ".actual_output."+inputsHash+"."+output.fileFormat);
            args.add(output.key + "=" + actualOutput);

            comparisonsToDoAfterRun.add(new OutputFileComparison(expectedOutput, actualOutput, comparator));

            if(!needToRegenerateExpectedOutputs && (!expectedOutput.canRead() || expectedOutput.lastModified() < newestInputModTime)) {
                log.info("Output file " + expectedOutput.getName() + " doesn't exist or is out-of-date compared to some of the inputs or " +
                         "the released CLP jar. Will regenerate it by running released CLP "+getReleasedCommandLineProgramJarPath() +
                        " with args: " + StringUtil.join("  ", releasedCLPArgs.toArray()));
                needToRegenerateExpectedOutputs = true;
                //don't break here so all output args will be collected
            }
        }

        //run the clean/released version of the CLP to generate expected outputs if necessary
        if(needToRegenerateExpectedOutputs) {
            final String command = "java -jar " + jarPath + " " + StringUtil.join(" ", releasedCLPArgs.toArray());
            ProcessExecutor.ExitStatusAndOutput result = ProcessExecutor.executeAndReturnInterleavedOutput(command);
            log.info(result.stdout);
            Assert.assertEquals(result.exitStatus, 0, "Release CLP command exit status should be 0. Command: " + command);
        }

        //run the dev version of the CLP being tested
        log.info("Running current version of " + getProgram() + " with args: " + StringUtil.join("  ", args.toArray()));
        int exitStatus = getProgram().instanceMain(args.toArray(new String[args.size()]));
        Assert.assertEquals(exitStatus, 0);

        //compare expected vs actual output
        for(OutputFileComparison comparison : comparisonsToDoAfterRun) {
            comparison.doComparison();
        }
    }


    /**
     * Represents a regular Picard command line arg of the form key=value (eg. SORT_ORDER=coordinate)
     */
    public static class Arg {
        public final String key;
        public final String value;
        public Arg(final String key, final String value) {
            this.key = key;
            this.value = value;
        }
    }


    /** Arg whose value is an input file path (eg. INPUT=foo.bam) */
    public static class InputFileArg {
        public final String key;
        public final File inputFilePath;

        public InputFileArg(final String key, final String inputFilePath) {
            this(key, new File(inputFilePath));
        }

        public InputFileArg(final String key, final File inputFilePath) {
            this.key = key;
            this.inputFilePath = inputFilePath;
        }
    }

    /** Arg whose value is an output file path (eg. OUTPUT=bar.bam) */
    public static class OutputFileArg {
        private final String key;
        private final String fileFormat;
        private final OutputFormatComparator comparator;

        public OutputFileArg(final String key, final String fileFormat) {
            this(key, fileFormat, null);
        }

        public OutputFileArg(final String key, final String fileFormat, OutputFormatComparator c) {
            this.key = key;
            this.fileFormat = fileFormat;
            this.comparator = c;
        }
    }

    /** Struct-like class used to remember to compare two specific files using the appropriate comparator  */
    private static class OutputFileComparison {
        private final File expected;
        private final File actual;
        private final OutputFormatComparator comparator;

        public OutputFileComparison(final File expected, final File actual, final OutputFormatComparator comparator) {
            this.expected = expected;
            this.actual = actual;
            this.comparator = comparator;
        }

        public void doComparison() {
            log.info("Comparing " + expected.getName() + " vs. " + actual.getName());
            comparator.assertAreIdentical(expected, actual);
        }
    }

    /**
     * Parent for classes that know how to compare expected vs. actual files for a particular file format.
     */
    public static interface OutputFormatComparator {

        public abstract void assertAreIdentical(File expected, File actual);
    }

    /**
     * Implements a generic comparator for text file formats.
     */
    public static class TextOutputFileComparator implements OutputFormatComparator {

        private final String[] linesToIgnoreRegexps;

        /**
         * Default constructor.
         */
        public TextOutputFileComparator() {
            this(new String[0]);
        }

        /**
         * Constructor.
         * @param linesToIgnoreRegexps A set of regular expressions where lines that matches any of these regular expressions are ignored
         *                             and don't generate errors if they differ between the expected & actual output files.
         */
        public TextOutputFileComparator(final String[] linesToIgnoreRegexps) {
            this.linesToIgnoreRegexps = linesToIgnoreRegexps;  //TODO should probably use Pattern to pre-compile the regexps
        }

        @Override
        public void assertAreIdentical(final File expected, final File actual) {
            try {
                final BufferedReader expectedStream = new BufferedReader(new FileReader(expected));
                final BufferedReader actualStream = new BufferedReader(new FileReader(actual));
                String expectedLine = expectedStream.readLine(), actualLine = actualStream.readLine();
                while(expectedLine != null && actualLine != null) {
                    boolean ignoreThisLine = false;
                    for(final String regexp : linesToIgnoreRegexps) {
                        if(expectedLine.matches(regexp) && actualLine.matches(regexp)) {
                            ignoreThisLine = true;
                            break;
                        }
                    }
                    if(!ignoreThisLine)
                        compareLine(expectedLine, actualLine);

                    expectedLine = expectedStream.readLine();
                    actualLine = actualStream.readLine();
                }
                expectedStream.close();
                actualStream.close();

                Assert.assertFalse(expectedLine == null && actualLine != null, "Actual output file "+actual+" has more lines than expected output file "+expected);
                Assert.assertFalse(actualLine == null && expectedLine != null, "Expected output file "+expected+" has more lines than the actual output file "+actual);
            } catch(IOException e) {
                throw new PicardException("Unexpected I/O exception while reading files: " + expected + ", " + actual, e);
            }
        }

        protected void compareLine(final String expectedLine, final String actualLine) {
            Assert.assertEquals(expectedLine, actualLine);
        }
    }


    public static class BamOutputFileComparator implements OutputFormatComparator {
        @Override
        public void assertAreIdentical(final File expected, final File actual) {
            throw new IllegalStateException("Not yet implemented");
        }
    }


}