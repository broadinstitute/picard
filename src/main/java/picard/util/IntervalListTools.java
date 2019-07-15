/*
 * The MIT License
 *
 * Copyright (c) 2010-2018 The Broad Institute
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

package picard.util;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineParser.ClpEnum;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.IntervalsManipulationProgramGroup;
import picard.util.IntervalList.IntervalListScatter;
import picard.util.IntervalList.IntervalListScatterMode;
import picard.util.IntervalList.IntervalListScatterer;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.BinaryOperator;

/**
 * Performs various {@link IntervalList} manipulations.
 * <p>
 * <h3>Summary</h3>
 * This tool offers multiple interval list file manipulation capabilities, including: sorting,
 * merging, subtracting, padding, and other set-theoretic operations. The default action
 * is to merge and sort the intervals provided in the {@link #INPUT}s. Other options, e.g. interval subtraction, are
 * controlled by the arguments.
 * <br />
 * Both {@link IntervalList} and VCF files are accepted as input. {@link IntervalList} should be denoted with the extension
 * {@value htsjdk.samtools.util.IOUtil#INTERVAL_LIST_FILE_EXTENSION}, while a VCF must have one of {@value htsjdk.samtools.util.IOUtil#VCF_FILE_EXTENSION}, {@value htsjdk.samtools.util.IOUtil#COMPRESSED_VCF_FILE_EXTENSION},
 * {@value htsjdk.samtools.util.IOUtil#BCF_FILE_EXTENSION}. When VCF file is used as input, each variant is translated into an using its reference allele or the END
 * INFO annotation (if present) to determine the extent of the interval.
 * <p>
 * {@link IntervalListTools} can also "scatter" the resulting interval-list into many interval-files. This can be useful
 * for creating multiple interval lists for scattering an analysis over.
 * <p>
 * <h3>Details</h3>
 * The IntervalList file format is designed to help the users avoid mixing references when supplying intervals and
 * other genomic data to a single tool. A SAM style header must be present at the top of the file. After the header,
 * the file then contains records, one per line in text format with the following
 * values tab-separated:
 * <pre>
 * <ul>
 * <li>Sequence name (SN)</li>
 * <li>Start position (1-based)</li>
 * <li>End position (1-based, end inclusive)</li>
 * <li>Strand (either + or -)</li>
 * <li>Interval name (ideally unique names for intervals)</li>
 * </ul>
 * </pre>
 * The coordinate system is 1-based, closed-ended, so that the first base in a sequence is at position 1, and both the start
 * and the end positions are included in an interval.
 * <p>
 * For Example:
 * <pre>
 * \@HD	VN:1.0
 * \@SQ	SN:chr1	LN:501
 * \@SQ	SN:chr2	LN:401
 * chr1	1	100	+	starts at the first base of the contig and covers 100 bases
 * chr2	100	100	+	interval with exactly one base
 * </pre>
 * <p>
 * <h3>Usage examples</h3>
 * <h4>1. Combine the intervals from two interval lists:</h4>
 * <pre>
 * java -jar picard.jar IntervalListTools \\
 *       ACTION=CONCAT \\
 *       I=input.interval_list \\
 *       I=input_2.interval_list \\
 *       O=new.interval_list
 * </pre>
 * <p>
 * <h4>2. Combine the intervals from two interval lists, sorting the resulting in list and merging overlapping and abutting
 * intervals:</h4>
 * <pre>
 * java -jar picard.jar IntervalListTools \\
 *       ACTION=CONCAT \\
 *       SORT=true \\
 *       UNIQUE=true \\
 *       I=input.interval_list \\
 *       I=input_2.interval_list \\
 *       O=new.interval_list
 * </pre>
 * <p>
 * <h4>3. Subtract the intervals in SECOND_INPUT from those in INPUT:</h4>
 * <pre>
 * java -jar picard.jar IntervalListTools \\
 *       ACTION=SUBTRACT \\
 *       I=input.interval_list \\
 *       SI=input_2.interval_list \\
 *       O=new.interval_list
 * </pre>
 * <p>
 * <h4>4. Find bases that are in either input1.interval_list or input2.interval_list, and also in input3.interval_list:</h4>
 * <pre>
 * java -jar picard.jar IntervalListTools \\
 *       ACTION=INTERSECT \\
 *       I=input1.interval_list \\
 *       I=input2.interval_list \\
 *       SI=input3.interval_list \\
 *       O=new.interval_list
 * </pre>
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = IntervalListTools.USAGE_SUMMARY + IntervalListTools.USAGE_DETAILS,
        oneLineSummary = IntervalListTools.USAGE_SUMMARY,
        programGroup = IntervalsManipulationProgramGroup.class
)
@DocumentedFeature
public class IntervalListTools extends CommandLineProgram {
    static final String USAGE_SUMMARY = "A tool for performing various IntervalList manipulations";
    static final String USAGE_DETAILS =
            " <h3>Summary</h3>" +
                    "This tool offers multiple interval list file manipulation capabilities, including: sorting, " +
                    "merging, subtracting, padding, and other set-theoretic operations. The default action " +
                    "is to merge and sort the intervals provided in the INPUTs. Other options, e.g. interval subtraction, are " +
                    "controlled by the arguments." +
                    "<br />" +
                    "Both IntervalList and VCF files are accepted as input. IntervalList should be denoted with the extension " +
                    ".interval_list, while a VCF must have one of .vcf, .vcf.gz, .bcf " +
                    "When VCF file is used as input, each variant is translated into an using its reference allele or the END " +
                    "INFO annotation (if present) to determine the extent of the interval. " +
                    "\n" +
                    "IntervalListTools can also \"scatter\" the resulting interval-list into many interval-files. This can be useful " +
                    "for creating multiple interval lists for scattering an analysis over.\n" +
                    "\n" +
                    " <h3>Details</h3> " +
                    "The IntervalList file format is designed to help the users avoid mixing references when supplying intervals and " +
                    "other genomic data to a single tool. A SAM style header must be present at the top of the file. After the header, " +
                    "the file then contains records, one per line in text format with the following" +
                    "values tab-separated: \n" +
                    "\n" +
                    " - Sequence name (SN) \n" +
                    " - Start position (1-based)\n" +
                    " - End position (1-based, inclusive)\n" +
                    " - Strand (either + or -)\n" +
                    " - Interval name (ideally unique names for intervals)\n" +
                    "\n" +
                    "The coordinate system is 1-based, closed-ended so that the first base in a sequence has position 1, and both the start " +
                    "and the end positions are included in an interval.\n" +
                    "\n" +
                    "Example interval list file" +
                    "<pre>" +
                    "@HD\tVN:1.0\n" +
                    "@SQ\tSN:chr1\tLN:501\n" +
                    "@SQ\tSN:chr2\tLN:401\n" +
                    "chr1\t1\t100\t+\tstarts at the first base of the contig and covers 100 bases\n" +
                    "chr2\t100\t100\t+\tinterval with exactly one base\n" +
                    "</pre>" +
                    "\n" +
                    "\n" +
                    "<h3>Usage Examples</h3>" +
                    "<h4>1. Combine the intervals from two interval lists:</h4>" +
                    "<pre>" +
                    "java -jar picard.jar IntervalListTools \\\n" +
                    "      ACTION=CONCAT \\\n" +
                    "      I=input.interval_list \\\n" +
                    "      I=input_2.interval_list \\\n" +
                    "      O=new.interval_list" +
                    "</pre>" +
                    "" +
                    " <h4>2. Combine the intervals from two interval lists, sorting the resulting in list and merging overlapping and abutting " +
                    "intervals:</h4>" +
                    " <pre>" +
                    " java -jar picard.jar IntervalListTools \\\n" +
                    "       ACTION=CONCAT \\\n" +
                    "       SORT=true \\\n" +
                    "       UNIQUE=true \\\n" +
                    "       I=input.interval_list \\\n" +
                    "       I=input_2.interval_list \\\n" +
                    "       O=new.interval_list" +
                    " </pre>" +
                    "" +
                    " <h4>3. Subtract the intervals in SECOND_INPUT from those in INPUT</h4>" +
                    " <pre>" +
                    " java -jar picard.jar IntervalListTools \\\n" +
                    "       ACTION=SUBTRACT \\\n" +
                    "       I=input.interval_list \\\n" +
                    "       SI=input_2.interval_list \\\n" +
                    "       O=new.interval_list" +
                    " </pre>" +
                    "" +
                    " <h4>4. Find bases that are in either input1.interval_list or input2.interval_list, and also in input3.interval_list:</h4>" +
                    " <pre>" +
                    " java -jar picard.jar IntervalListTools \\\n" +
                    "       ACTION=INTERSECT \\\n" +
                    "       I=input1.interval_list \\\n" +
                    "       I=input2.interval_list \\\n" +
                    "       SI=input3.interval_list \\\n" +
                    "       O=new.interval_list" +
                    " </pre>" +
                    "" +
                    "";

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "One or more interval lists. If multiple interval lists are provided the output is the" +
                    "result of merging the inputs. Supported formats are interval_list and VCF.", minElements = 1)
    public List<File> INPUT;

    @Argument(doc = "The output interval list file to write (if SCATTER_COUNT == 1) or the directory into which " +
            "to write the scattered interval sub-directories (if SCATTER_COUNT > 1).", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional = true)
    public File OUTPUT;

    @Argument(doc = "The amount to pad each end of the intervals by before other operations are undertaken. Negative numbers are allowed " +
            "and indicate intervals should be shrunk. Resulting intervals < 0 bases long will be removed. Padding is applied to the " +
            "interval lists (both INPUT and SECOND_INPUT, if provided) <b> before </b> the ACTION is performed.", optional = true)
    public int PADDING = 0;

    @Argument(doc = "If true, merge overlapping and adjacent intervals to create a list of unique intervals. Implies SORT=true.")
    public boolean UNIQUE = false;

    @Argument(doc = "If true, sort the resulting interval list by coordinate.")
    public boolean SORT = true;

    @Argument(doc = "Action to take on inputs.")
    public Action ACTION = Action.CONCAT;

    @Argument(shortName = "SI", doc = "Second set of intervals for SUBTRACT and DIFFERENCE operations.", optional = true)
    public List<File> SECOND_INPUT;

    @Argument(doc = "One or more lines of comment to add to the header of the output file (as @CO lines in the SAM header).", optional = true)
    public List<String> COMMENT = null;

    @Argument(doc = "The number of files into which to scatter the resulting list by locus; in some situations, fewer intervals may be emitted.  ")
    public int SCATTER_COUNT = 1;

    @Argument(doc = "When scattering with this argument, each of the resultant files will (ideally) have this amount of 'content', which " +
            " means either base-counts or interval-counts depending on SUBDIVISION_MODE. When provided, overrides SCATTER_COUNT", optional = true)
    public Integer SCATTER_CONTENT = null;

    @Argument(doc = "Whether to include filtered variants in the vcf when generating an interval list from vcf.", optional = true)
    public boolean INCLUDE_FILTERED = false;

    @Argument(shortName = "BRK", doc = "If set to a positive value will create a new interval list with the original intervals" +
            " broken up at integer multiples of this value. Set to 0 to NOT break up intervals.", optional = true)
    public int BREAK_BANDS_AT_MULTIPLES_OF = 0;

    @Argument(shortName = "M", doc = "The mode used to scatter the interval list.")
    public IntervalListScatterMode SUBDIVISION_MODE = IntervalListScatterMode.INTERVAL_SUBDIVISION;

    @Argument(doc = "Produce the inverse list of intervals, that is, the regions in the genome that are <br>not</br> covered " +
            "by any of the input intervals. Will merge abutting intervals first. Output will be sorted.", optional = true)
    public boolean INVERT = false;

    @Argument(doc = "What value to output to COUNT_OUTPUT file or stdout (for scripting).  If COUNT_OUTPUT is provided, this parameter must not be NONE.")
    public Output OUTPUT_VALUE = Output.NONE;

    @Argument(doc = "File to which to print count of bases or intervals in final output interval list.  When not set, value indicated by OUTPUT_VALUE will be printed to stdout.  " +
            "If this parameter is set, OUTPUT_VALUE must not be NONE.", optional = true)
    public File COUNT_OUTPUT;

    enum Output {
        NONE {
            void output(final long totalBaseCount, final long intervalCount, final PrintStream stream){
                if (stream != System.out) {
                    throw new PicardException("Asked to output NONE with COUNT_OUTPUT set.");
                }
            }
        },
        BASES {
            void output(final long totalBaseCount, final long intervalCount, final PrintStream stream) {
                stream.println(totalBaseCount);
            }
        },
        INTERVALS {
            void output(final long totalBaseCount, final long intervalCount, final PrintStream stream) {
                stream.println(intervalCount);
            }
        };

        abstract void output(final long totalBaseCount, final long intervalCount, final PrintStream stream);
    }

    private static final Log LOG = Log.getInstance(IntervalListTools.class);

    public enum Action implements ClpEnum {

        CONCAT("The concatenation of all the intervals in all the INPUTs, no sorting or merging of overlapping/abutting " +
                "intervals implied. Will result in a possibly unsorted list unless requested otherwise.", false),

        UNION("Like CONCATENATE but with UNIQUE and SORT implied, the result being the set-wise union of all INPUTS, " +
                "with overlapping and abutting intervals merged into one.", false) {
            @Override
            IntervalList act(final IntervalList firstList, final IntervalList ignored) {
                return super.act(firstList.sorted().uniqued(), ignored);
            }
        },

        INTERSECT("The sorted and merged set of all loci that are contained in all of the INPUTs.", false) {
            @Override
            IntervalList reduceEach(final IntervalList list1, final IntervalList list2) {
                return IntervalList.intersection(CollectionUtil.makeList(list1, list2));
            }
        },
        SUBTRACT("Subtracts the intervals in SECOND_INPUT from those in INPUT. The resulting loci are those in INPUT that are not in SECOND_INPUT.", true) {
            @Override
            IntervalList act(final IntervalList list1, final IntervalList list2) {
                return IntervalList.subtract(list1, list2);
            }
        },
        SYMDIFF("Results in loci that are in INPUT or SECOND_INPUT but are not in both.", true) {
            @Override
            IntervalList act(final IntervalList list1, final IntervalList list2) {
                //todo: fix this when new htsjdk
                return IntervalList.difference(Collections.singleton(list1), Collections.singleton(list2));
            }
        },
        OVERLAPS("Outputs the entire intervals from INPUT that have bases which overlap any interval from SECOND_INPUT. " +
                "Note that this is different than INTERSECT in that each original interval is either emitted in its entirety, or not at all.", true) {
            @Override
            IntervalList act(final IntervalList list1, final IntervalList list2) {
                return IntervalList.overlaps(list1, list2);
            }
        };

        final String helpdoc;
        final boolean takesSecondInput;

        Action(final String helpdoc, boolean takesSecondInput) {
            this.helpdoc = helpdoc;
            this.takesSecondInput = takesSecondInput;
        }

        @Override
        public String getHelpDoc() {
            return helpdoc;
        }

        IntervalList act(final IntervalList firstList, final IntervalList secondList) {
            return firstList;
        }

        IntervalList reduceEach(final IntervalList list1, final IntervalList list2){
                return IntervalList.concatenate(CollectionUtil.makeList(list1, list2));
        }
    }

    @Override
    protected int doWork() {
        // Check inputs
        IOUtil.assertFilesAreReadable(INPUT);
        IOUtil.assertFilesAreReadable(SECOND_INPUT);
        if (COUNT_OUTPUT != null) {
            IOUtil.assertFileIsWritable(COUNT_OUTPUT);
        }

        // Read in the interval lists and apply any padding
        final IntervalList lists = openIntervalLists(INPUT, ACTION::reduceEach);

        // same for the second list
        final IntervalList secondLists = openIntervalLists(SECOND_INPUT, ACTION::reduceEach);

        if (UNIQUE && !SORT) {
            LOG.warn("UNIQUE=true requires sorting but SORT=false was specified.  Results will be sorted.");
        }

        final IntervalList result = ACTION.act(lists, secondLists);

        if (INVERT) {
            SORT = false; // no need to sort, since return will be sorted by definition.
            UNIQUE = true;
        }

        final IntervalList possiblySortedResult = SORT ? result.sorted() : result;
        final IntervalList possiblyInvertedResult = INVERT ? IntervalList.invert(possiblySortedResult) : possiblySortedResult;

        //only get unique if this has been asked unless inverting (since the invert will return a unique list)
        List<Interval> finalIntervals = UNIQUE ? possiblyInvertedResult.uniqued().getIntervals() : possiblyInvertedResult.getIntervals();

        if (BREAK_BANDS_AT_MULTIPLES_OF > 0) {
            finalIntervals = IntervalList.breakIntervalsAtBandMultiples(finalIntervals, BREAK_BANDS_AT_MULTIPLES_OF);
        }

        // Decide on a PG ID and make a program group
        final SAMFileHeader header = result.getHeader();
        final Set<String> pgs = new HashSet<>();
        for (final SAMProgramRecord pg : header.getProgramRecords()) pgs.add(pg.getId());
        for (int i = 1; i < Integer.MAX_VALUE; ++i) {
            if (!pgs.contains(String.valueOf(i))) {
                final SAMProgramRecord pg = new SAMProgramRecord(String.valueOf(i));
                pg.setCommandLine(getCommandLine());
                pg.setProgramName(getClass().getSimpleName());
                header.addProgramRecord(pg);
                break;
            }
        }

        // Add any comments
        if (COMMENT != null) {
            for (final String comment : COMMENT) {
                header.addComment(comment);
            }
        }

        final IntervalList output = new IntervalList(header);
        finalIntervals.forEach(output::add);

        final ScatterSummary resultIntervals;

        if (SCATTER_CONTENT != null) {
            final long listSize = SUBDIVISION_MODE.make().listWeight(output);
            SCATTER_COUNT = (int) listSize / SCATTER_CONTENT;
            LOG.info(String.format("Using SCATTER_CONTENT = %d and an interval of size %d, attempting to scatter into %s intervals.", SCATTER_CONTENT, listSize, SCATTER_COUNT));
        }

        if (OUTPUT != null && SCATTER_COUNT > 1) {

            IOUtil.assertDirectoryIsWritable(OUTPUT);

            final ScatterSummary scattered = writeScatterIntervals(output);
            LOG.info(String.format("Wrote %s scatter subdirectories to %s.", scattered.size, OUTPUT));
            if (scattered.size != SCATTER_COUNT) {
                LOG.warn(String.format(
                        "Requested scatter width of %s, but only emitted %s.  (This may be an expected consequence of running in %s mode.)",
                        SCATTER_COUNT,
                        scattered.size,
                        SUBDIVISION_MODE
                ));
            }
            resultIntervals = scattered;

        } else {
            if (OUTPUT != null) {
                output.write(OUTPUT);
            }
            resultIntervals = new ScatterSummary();
            resultIntervals.size = 1;
            resultIntervals.intervalCount = output.getIntervals().size();
            resultIntervals.baseCount = output.getBaseCount();
        }

        LOG.info("Produced " + resultIntervals.intervalCount + " intervals totalling " + resultIntervals.baseCount + " bases.");
        if (COUNT_OUTPUT != null) {
            try (final PrintStream countStream = new PrintStream(COUNT_OUTPUT)) {
                OUTPUT_VALUE.output(resultIntervals.baseCount, resultIntervals.intervalCount, countStream);
            }
            catch (final IOException e) {
                throw new PicardException("There was a problem writing count to " + COUNT_OUTPUT.getAbsolutePath());
            }
        } else {
            OUTPUT_VALUE.output(resultIntervals.baseCount, resultIntervals.intervalCount, System.out);
        }
        return 0;
    }

    private IntervalList openIntervalLists(final List<File> files, BinaryOperator<IntervalList> accumulator ) {
        return files.stream()
                .map(f->IntervalListInputType.getIntervalList(f, INCLUDE_FILTERED).padded(PADDING))
                .reduce(accumulator)
                .orElse(null);
    }

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errorMsgs = new ArrayList<>();
        if (SCATTER_COUNT < 1) {
            errorMsgs.add("SCATTER_COUNT must be greater than 0.");
        }
        if (BREAK_BANDS_AT_MULTIPLES_OF < 0) {
            errorMsgs.add("BREAK_BANDS_AT_MULTIPLES_OF must be greater than or equal to 0.");
        }
        if ((SECOND_INPUT == null || SECOND_INPUT.isEmpty()) && ACTION.takesSecondInput) {
            errorMsgs.add("SECOND_INPUT was not provided but action " + ACTION + " requires a second input.");
        }
        if ((SECOND_INPUT != null && !SECOND_INPUT.isEmpty()) && !ACTION.takesSecondInput) {
            errorMsgs.add("SECOND_INPUT was provided but action " + ACTION + " doesn't take a second input.");
        }
        if (COUNT_OUTPUT != null && OUTPUT_VALUE == Output.NONE) {
            errorMsgs.add("COUNT_OUTPUT was provided but OUTPUT_VALUE is set to NONE.");
        }

        return errorMsgs.isEmpty() ? null : errorMsgs.toArray(new String[0]);
    }


    static class ScatterSummary {
        public long size;
        public long baseCount;
        public long intervalCount;
    }
    /**
     * Method to scatter an interval list by locus.
     *
     * @param list The list of intervals to scatter
     * @return The scattered intervals, represented as a {@link List} of {@link IntervalList}
     */
    private  ScatterSummary writeScatterIntervals(final IntervalList list) {
        final IntervalListScatterer scatterer = SUBDIVISION_MODE.make();
        final IntervalListScatter scatter = new IntervalListScatter(scatterer, list, SCATTER_COUNT);

        final DecimalFormat fileNameFormatter = new DecimalFormat("0000");
        int fileIndex = 1;
        final ScatterSummary summary = new ScatterSummary();
        for (final IntervalList intervals : scatter) {
            summary.size++;
            summary.baseCount += intervals.getBaseCount();
            summary.intervalCount += intervals.getIntervals().size();
            intervals.write(createDirectoryAndGetScatterFile(OUTPUT, SCATTER_COUNT, fileNameFormatter.format(fileIndex++)));
        }

        return summary;
    }

    private static File getScatteredFileName(final File scatterDirectory, final long scatterTotal, final String formattedIndex) {
        return new File(scatterDirectory.getAbsolutePath() + "/temp_" + formattedIndex + "_of_" +
                scatterTotal + "/scattered" + IntervalList.INTERVAL_LIST_FILE_EXTENSION);
    }

    private static File createDirectoryAndGetScatterFile(final File outputDirectory, final long scatterCount, final String formattedIndex) {
        createDirectoryOrFail(outputDirectory);
        final File result = getScatteredFileName(outputDirectory, scatterCount, formattedIndex);
        createDirectoryOrFail(result.getParentFile());
        return result;
    }

    private static void createDirectoryOrFail(final File directory) {
        if (!directory.exists() && !directory.mkdir()) {
            throw new PicardException("Unable to create directory: " + directory.getAbsolutePath());
        }
    }

    enum IntervalListInputType {
        VCF(IOUtil.VCF_EXTENSIONS) {
            @Override
            protected IntervalList getIntervalListInternal(final File vcf, final boolean includeFiltered) {
                return VCFFileReader.fromVcf(vcf, includeFiltered);
            }
        },
        INTERVAL_LIST(IOUtil.INTERVAL_LIST_FILE_EXTENSION) {
            @Override
            protected IntervalList getIntervalListInternal(final File intervalList, final boolean includeFiltered) {
                return IntervalList.fromFile(intervalList);
            }
        };

        protected final Collection<String> applicableExtensions;

        IntervalListInputType(final String... s) {
            applicableExtensions = CollectionUtil.makeSet(s);
        }

        IntervalListInputType(final Collection<String> extensions) {
            applicableExtensions = extensions;
        }

        protected abstract IntervalList getIntervalListInternal(final File file, final boolean includeFiltered);

        static IntervalListInputType forFile(final File intervalListExtractable) {
            for (final IntervalListInputType intervalListInputType : IntervalListInputType.values()) {
                for (final String s : intervalListInputType.applicableExtensions) {
                    if (intervalListExtractable.getName().endsWith(s)) {
                        return intervalListInputType;
                    }
                }
            }
            throw new SAMException("Cannot figure out type of file " + intervalListExtractable.getAbsolutePath() + " from extension. Current implementation understands the following types: " + Arrays.toString(IntervalListInputType.values()));
        }

        public static IntervalList getIntervalList(final File file, final boolean includeFiltered) {
            return forFile(file).getIntervalListInternal(file, includeFiltered);
        }

        @Override
        public String toString() {
            return super.toString() + ": " + applicableExtensions.toString();
        }
    }
}
