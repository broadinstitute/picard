package picard.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import picard.PicardException;
import picard.cmdline.CommandLineParser;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.Usage;

import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Little class to aid working with interval lists.
 *
 * @author Tim Fennell
 */
public class IntervalListTools extends CommandLineProgram {
    @Usage
    public final String USAGE = getStandardUsagePreamble() + " General tool for manipulating interval lists, " +
            "including sorting, merging, padding, uniqueifying, and other set-theoretic operations. Default operation if given one or more inputs is to " +
            "merge and sort them.  Other options are controlled by arguments.";

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "One or more interval lists. If multiple interval lists are provided the output is the" +
                    "result of merging the inputs.", minElements = 1)
    public List<File> INPUT;

    @Option(doc = "The output interval list file to write (if SCATTER_COUNT is 1) or the directory into which " +
            "to write the scattered interval sub-directories (if SCATTER_COUNT > 1)", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional = true)
    public File OUTPUT;

    @Option(doc = "The amount to pad each end of the intervals by before other operations are undertaken. Negative numbers are allowed " +
            "and indicate intervals should be shrunk. Resulting intervals < 0 bases long will be removed. Padding is applied to the interval lists <b> before </b> the ACTION is performed.", optional = true)
    public int PADDING = 0;

    @Option(doc = "If true, merge overlapping and adjacent intervals to create a list of unique intervals. Implies SORT=true")
    public boolean UNIQUE = false;

    @Option(doc = "If true, sort the resulting interval list by coordinate.")
    public boolean SORT = true;

    @Option(doc = "Action to take on inputs.")
    public Action ACTION = Action.CONCAT;

    @Option(shortName = "SI", doc = "Second set of intervals for SUBTRACT and DIFFERENCE operations.", optional = true)
    public List<File> SECOND_INPUT;

    @Option(doc = "One or more lines of comment to add to the header of the output file.", optional = true)
    public List<String> COMMENT = null;

    @Option(doc = "The number of files into which to scatter the resulting list by locus; in some situations, fewer intervals may be emitted.")
    public int SCATTER_COUNT = 1;

    @Option(shortName = "M", doc = "Do not subdivide ")
    public IntervalListScatterer.Mode SUBDIVISION_MODE = IntervalListScatterer.Mode.INTERVAL_SUBDIVISION;

    @Option(doc = "Produce the inverse list", optional = true)
    public boolean INVERT = false;

    private static final Log LOG = Log.getInstance(IntervalListTools.class);

    public enum Action implements CommandLineParser.ClpEnum {

        CONCAT("The concatenation of all the INPUTs, no sorting or merging of overlapping/abutting intervals implied. Will result in an unsorted list unless requested otherwise.") {
            @Override
            IntervalList act(final List<IntervalList> list, final List<IntervalList> unused) {
                if (!unused.isEmpty())
                    throw new IllegalArgumentException(String.format("Second List found when action was %s. Ignoring second list.", this.name()));
                return IntervalList.concatenate(list);
            }
        },
        UNION("Like CONCATENATE but with UNIQUE and SORT implied, the result being the set-wise union of all INPUTS.") {
            @Override
            IntervalList act(final List<IntervalList> list, final List<IntervalList> unused) {
                if (!unused.isEmpty())
                    throw new IllegalArgumentException(String.format("Second List found when action was %s. Ignoring second list.", this.name()));
                return IntervalList.union(list);
            }
        },
        INTERSECT("The sorted, uniqued set of all loci that are contained in all of the INPUTs.") {
            @Override
            IntervalList act(final List<IntervalList> list, final List<IntervalList> unused) {
                if (!unused.isEmpty())
                    throw new IllegalArgumentException(String.format("Second List found when action was %s. Ignoring second list.", this.name()));
                return IntervalList.intersection(list);
            }
        },
        SUBTRACT("Subtracts SECOND_INPUT from INPUT. The resulting loci are there in INPUT that are not in SECOND_INPUT") {
            @Override
            IntervalList act(final List<IntervalList> list1, final List<IntervalList> list2) {
                return IntervalList.subtract(list1, list2);

            }
        },
        SYMDIFF("Find loci that are in INPUT or SECOND_INPUT but are not in both.") {
            @Override
            IntervalList act(final List<IntervalList> list1, final List<IntervalList> list2) {
                return IntervalList.difference(list1, list2);
            }
        };


        String helpdoc;

        Action(final String helpdoc) {
            this.helpdoc = helpdoc;
        }

        @Override
        public String getHelpDoc() {
            return helpdoc;
        }

        abstract IntervalList act(final List<IntervalList> list1, final List<IntervalList> list2);

    }

    // Stock main method
    public static void main(final String[] args) {
        new IntervalListTools().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        // Check inputs
        for (final File f : INPUT) IOUtil.assertFileIsReadable(f);
        for (final File f : SECOND_INPUT) IOUtil.assertFileIsReadable(f);

        if (OUTPUT != null) {
            if (SCATTER_COUNT == 1) {
                IOUtil.assertFileIsWritable(OUTPUT);
            } else {
                IOUtil.assertDirectoryIsWritable(OUTPUT);
            }
        }

        // Read in the interval lists and apply any padding
        final List<IntervalList> lists = new ArrayList<IntervalList>();
        for (final File f : INPUT) {
            final IntervalList list = IntervalList.fromFile(f);
            if (PADDING != 0) {
                final IntervalList out = new IntervalList(list.getHeader());
                for (final Interval i : list) {
                    final int start = i.getStart() - PADDING;
                    final int end = i.getEnd() + PADDING;
                    if (start <= end) {
                        final Interval i2 = new Interval(i.getSequence(), start, end, i.isNegativeStrand(), i.getName());
                        out.add(i2);
                    }
                }

                lists.add(out);
            } else {
                lists.add(list);
            }
        }

        // same for the second list
        final List<IntervalList> secondLists = new ArrayList<IntervalList>();
        for (final File f : SECOND_INPUT) {
            final IntervalList list = IntervalList.fromFile(f);
            if (PADDING != 0) {
                final IntervalList out = new IntervalList(list.getHeader());
                for (final Interval i : list) {
                    final int start = i.getStart() - PADDING;
                    final int end = i.getEnd() + PADDING;
                    if (start <= end) {
                        final Interval i2 = new Interval(i.getSequence(), start, end, i.isNegativeStrand(), i.getName());
                        out.add(i2);
                    }
                }

                secondLists.add(out);
            } else {
                secondLists.add(list);
            }
        }

        if (UNIQUE && !SORT) {
            LOG.warn("UNIQUE=true requires sorting but SORT=false was specified.  Results will be sorted!");
        }

        final IntervalList result = ACTION.act(lists, secondLists);

        if (INVERT) {
            SORT = false; // no need to sort, since return will be sorted by definition.
            UNIQUE = false; //no need to unique since invert will already return a unique list.
        }

        final IntervalList possiblySortedResult = SORT ? result.sorted() : result;
        final IntervalList possiblyInvertedResult = INVERT ? IntervalList.invert(possiblySortedResult) : possiblySortedResult;

        //only get unique if this has been asked unless inverting (since the invert will return a unique list)
        final List<Interval> finalIntervals = UNIQUE ? possiblyInvertedResult.uniqued().getIntervals() : possiblyInvertedResult.getIntervals();


        // Decide on a PG ID and make a program group
        final SAMFileHeader header = result.getHeader();
        final Set<String> pgs = new HashSet<String>();
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
        for (final Interval i : finalIntervals) {
            output.add(i);
        }

        final List<IntervalList> resultIntervals;
        if (OUTPUT != null) {
            if (SCATTER_COUNT == 1) {
                output.write(OUTPUT);
                resultIntervals = Arrays.asList(output);
            } else {
                final List<IntervalList> scattered = writeScatterIntervals(output);
                LOG.info(String.format("Wrote %s scatter subdirectories to %s.", scattered.size(), OUTPUT));
                if (scattered.size() != SCATTER_COUNT) {
                    LOG.warn(String.format(
                            "Requested scatter width of %s, but only emitted %s.  (This may be an expected consequence of running in %s mode.)",
                            SCATTER_COUNT,
                            scattered.size(),
                            IntervalListScatterer.Mode.BALANCING_WITHOUT_INTERVAL_SUBDIVISION
                    ));
                }
                resultIntervals = scattered;
            }
        } else {
            resultIntervals = Arrays.asList(output);
        }

        long totalUniqueBaseCount = 0;
        long intervalCount = 0;
        for (final IntervalList finalInterval : resultIntervals) {
            totalUniqueBaseCount = finalInterval.getUniqueBaseCount();
            intervalCount += finalInterval.size();    
        }

        LOG.info("Produced " + intervalCount + " intervals totalling " + totalUniqueBaseCount + " unique bases.");
        
        return 0;
    }

    @Override
    protected String[] customCommandLineValidation() {
        if (SCATTER_COUNT < 1) {
            return new String[]{"SCATTER_COUNT must be greater than 0."};
        }
        return null;
    }

    /**
     * Method to scatter an interval list by locus.
     *
     * @param list The list of intervals to scatter
     * @return The scattered intervals, represented as a {@link List} of {@link IntervalList}
     */
    private List<IntervalList> writeScatterIntervals(final IntervalList list) {
        final IntervalListScatterer scatterer = new IntervalListScatterer(SUBDIVISION_MODE);
        final List<IntervalList> scattered = scatterer.scatter(list, SCATTER_COUNT);

        final DecimalFormat fileNameFormatter = new DecimalFormat("0000");
        int fileIndex = 1;
        for (final IntervalList intervals : scattered) {
            intervals.write(createDirectoryAndGetScatterFile(OUTPUT, scattered.size(), fileNameFormatter.format(fileIndex++)));
        }

        return scattered;
    }

    public static File getScatteredFileName(final File scatterDirectory, final long scatterTotal, final String formattedIndex) {
        return new File(scatterDirectory.getAbsolutePath() + "/temp_" + formattedIndex + "_of_" +
                scatterTotal + "/scattered.intervals");

    }

    private static File createDirectoryAndGetScatterFile(final File outputDirectory, final long scatterCount, final String formattedIndex) {
        createDirectoryOrFail(outputDirectory);
        final File result = getScatteredFileName(outputDirectory, scatterCount, formattedIndex);
        createDirectoryOrFail(result.getParentFile());
        return result;
    }

    private static void createDirectoryOrFail(final File directory) {
        if (!directory.exists()) {
            if (!directory.mkdir()) {
                throw new PicardException("Unable to create directory: " + directory.getAbsolutePath());
            }
        }
    }
}
