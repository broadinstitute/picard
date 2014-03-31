package net.sf.picard.util;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineParser;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMProgramRecord;

import java.io.File;
import java.text.DecimalFormat;
import java.util.*;

/**
 * Little class to aid working with interval lists.
 *
 * @author Tim Fennell
 */
public class IntervalListTools extends CommandLineProgram {
    @Usage public final String USAGE = getStandardUsagePreamble() + " General tool for manipulating interval lists, " +
            "including sorting, merging, padding, uniqueifying, and other set-theoretic operations. Default operation if given one or more inputs is to " +
            "merge and sort them.  Other options are controlled by arguments.";

    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc="One or more interval lists. If multiple interval lists are provided the output is the" +
                "result of merging the inputs.", minElements = 1)
    public List<File> INPUT;

    @Option(doc="The output interval list file to write (if SCATTER_COUNT is 1) or the directory into which " +
            "to write the scattered interval sub-directories (if SCATTER_COUNT > 1)", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional = true)
    public File OUTPUT;

    @Option(doc="The amount to pad each end of the intervals by before other operations are undertaken. Negative numbers are allowed " +
            "and indicate intervals should be shrunk. Resulting intervals < 0 bases long will be removed. Padding is applied to the interval lists <b> before </b> the ACTION is performed.", optional=true)
    public int PADDING = 0;

    @Option(doc="If true, merge overlapping and adjacent intervals to create a list of unique intervals. Implies SORT=true")
    public boolean UNIQUE = false;

    @Option(doc="If true, sort the resulting interval list by coordinate.")
    public boolean SORT = true;

    @Option(doc = "Action to take on inputs.")
    public Action ACTION = Action.CONCAT;

    @Option(shortName = "SI", doc = "Second set of intervals for SUBTRACT and DIFFERENCE operations.", optional = true)
    public List<File> SECOND_INPUT;

    @Option(doc="One or more lines of comment to add to the header of the output file.", optional = true)
    public List<String> COMMENT = null;

    @Option(doc="The number of files into which to scatter the resulting list by locus.")
    public int SCATTER_COUNT = 1;

    @Option(doc = "Produce the inverse list", optional = true)
    public boolean INVERT = false;

    private static final Log log = Log.getInstance(IntervalListTools.class);

    public enum Action implements CommandLineParser.ClpEnum{

        CONCAT("The concatenation of all the INPUTs, no sorting or merging of overlapping/abutting intervals implied. Will result in an unsorted list unless requested otherwise.") {
            @Override
            IntervalList act(final List<IntervalList> list, final List<IntervalList> unused) {
                if(!unused.isEmpty()) throw new IllegalArgumentException(String.format("Second List found when action was %s. Ignoring second list.",this.name()));
                return IntervalList.concatenate(list);
            }
        },
        UNION ("Like CONCATENATE but with UNIQUE and SORT implied, the result being the set-wise union of all INPUTS.") {
            @Override
            IntervalList act(final List<IntervalList> list, final List<IntervalList> unused) {
                if(!unused.isEmpty()) throw new IllegalArgumentException(String.format("Second List found when action was %s. Ignoring second list.",this.name()));
                return IntervalList.union(list);
            }
        },
        INTERSECT ("The sorted, uniqued set of all loci that are contained in all of the INPUTs.") {
            @Override
            IntervalList act(final List<IntervalList> list, final List<IntervalList> unused) {
                if(!unused.isEmpty()) throw new IllegalArgumentException(String.format("Second List found when action was %s. Ignoring second list.",this.name()));
                return IntervalList.intersection(list);
            }
        },
       SUBTRACT ("Subtracts SECOND_INPUT from INPUT. The resulting loci are there in INPUT that are not in SECOND_INPUT") {
            @Override
            IntervalList act(final List<IntervalList> list1, final List<IntervalList> list2) {
                return IntervalList.subtract(list1, list2);

                }
        },
        SYMDIFF ("Find loci that are in INPUT or SECOND_INPUT but are not in both." ) {
            @Override
            IntervalList act(final List<IntervalList> list1, final List<IntervalList> list2) {
                return IntervalList.difference(list1, list2);
            }
        };


        String helpdoc;
        Action(final String helpdoc){
            this.helpdoc=helpdoc;
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
        for (final File f : INPUT) IoUtil.assertFileIsReadable(f);
        for (final File f : SECOND_INPUT) IoUtil.assertFileIsReadable(f);

        if (OUTPUT != null) {
            if (SCATTER_COUNT == 1) {
                IoUtil.assertFileIsWritable(OUTPUT);
            }
            else {
                IoUtil.assertDirectoryIsWritable(OUTPUT);
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
                    final int end   = i.getEnd()   + PADDING;
                    if (start <= end) {
                        final Interval i2 = new Interval(i.getSequence(), start, end, i.isNegativeStrand(), i.getName());
                        out.add(i2);
                    }
                }

                lists.add(out);
            }
            else {
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
                    final int end   = i.getEnd()   + PADDING;
                    if (start <= end) {
                        final Interval i2 = new Interval(i.getSequence(), start, end, i.isNegativeStrand(), i.getName());
                        out.add(i2);
                    }
                }

                secondLists.add(out);
            }
            else {
                secondLists.add(list);
            }
        }

        if (UNIQUE && !SORT ) {
            log.warn("UNIQUE=true requires sorting but SORT=false was specified.  Results will be sorted!");
        }

        final IntervalList result = ACTION.act(lists, secondLists);

        if(INVERT){
            SORT=false; // no need to sort, since return will be sorted by definition.
            UNIQUE=false; //no need to unique since invert will already return a unique list.
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
        long total = 0;
        for (final Interval i : finalIntervals) {
            output.add(i);
            total += i.length();
        }

        int intervals = 0;

        if (OUTPUT != null) {
            if (SCATTER_COUNT == 1) {
                intervals = output.size();
                output.write(OUTPUT);
            }
            else {
                intervals = scatterIntervals(output);
            }
        }

        log.info("Output " + intervals + " intervals totalling " + output.getUniqueBaseCount() + " unique bases.");

        return 0;
    }

    @Override
    protected String[] customCommandLineValidation() {
        if (SCATTER_COUNT < 1) {
            return new String[] {"SCATTER_COUNT must be greater than 0."};
        }
        return null;
    }

    /**
     * Method to scatter an interval list by locus.
     * @param list  The list of intervals to scatter
     * @return the number of intervals across the scattered lists (which may differ from the input
     *         as some may have been split)
     */
    private int scatterIntervals(final IntervalList list) {
        // algorithm (to match the GATK):
        // split = ()
        // set size = 0
        // pop the head H off locs.
        // If size + size(H) < splitSize:
        //      add H to split, continue
        // If size + size(H) == splitSize:
        //      done with split, put in splits, restart
        // if size + size(H) > splitSize:
        //      cut H into two pieces, first of which has splitSize - size bp
        //      push both pieces onto locs, continue
        // The last split is special -- when you have only one split left, it gets all of the remaining locs
        // to deal with rounding issues
        final IntervalList uniquedList=list.uniqued();
        final long idealSplitLength = Math.max((long)Math.floor(uniquedList.getBaseCount() / (1.0*SCATTER_COUNT)), 1);
        int splitLength = 0;
        IntervalList split = new IntervalList(uniquedList.getHeader());
        int index = 1;   // The index of the next split file to write
        final Iterator<Interval> it = uniquedList.iterator();
        int totalIntervals = 0;
        final DecimalFormat format = new DecimalFormat("0000");

        while (it.hasNext() && index < SCATTER_COUNT) {
            final Interval interval = it.next();
            int projectedSize = splitLength + interval.length();
            if (projectedSize < idealSplitLength) {
                split.add(interval);
                totalIntervals++;
                splitLength += interval.length();
            }
            else if (projectedSize == idealSplitLength) {
                split.add(interval);
                totalIntervals++;
                split.write(createDirectoryAndGetScatterFile(format.format(index++)));
                split = new IntervalList(uniquedList.getHeader());
                splitLength = 0;
            }
            else {
                int consumed = 0;
                while (projectedSize > idealSplitLength && index < SCATTER_COUNT) {
                    final int amountToConsume =(int)(idealSplitLength - splitLength);
                    final Interval partial = new Interval(interval.getSequence(), interval.getStart()+consumed,
                        interval.getStart()+consumed+amountToConsume-1, interval.isNegativeStrand(), interval.getName());
                    split.add(partial);
                    totalIntervals++;
                    split.write(createDirectoryAndGetScatterFile(format.format(index++)));
                    split = new IntervalList(uniquedList.getHeader());

                    consumed += amountToConsume;
                    splitLength = 0;
                    projectedSize = interval.length() - consumed;
                }

                // Add the remainder, if any, to the next split
                if (projectedSize > 0) {
                    final Interval remainder = new Interval(interval.getSequence(), interval.getStart()+consumed,
                            interval.getEnd(), interval.isNegativeStrand(), interval.getName());
                    split.add(remainder);
                    totalIntervals++;
                    splitLength = remainder.length();
                }
            }
        }
        // Write everything left to the last split
        while (it.hasNext()) {
            split.add(it.next());
            totalIntervals++;
        }
        split.write(createDirectoryAndGetScatterFile(format.format(index)));
        return totalIntervals;


    }

    public static File getScatteredFileName(final File scatterDirectory, final int scatterTotal, final String formattedIndex) {
        return new File(scatterDirectory.getAbsolutePath() + "/temp_" + formattedIndex + "_of_" +
                scatterTotal + "/scattered.intervals");

    }

    private File createDirectoryAndGetScatterFile(final String formattedIndex) {
        createDirectoryOrFail(OUTPUT);
        final File result = getScatteredFileName(OUTPUT, SCATTER_COUNT, formattedIndex);
        createDirectoryOrFail(result.getParentFile());
        return result;
    }

    private void createDirectoryOrFail(final File directory) {
        if (!directory.exists()) {
            if (!directory.mkdir()) {
                throw new PicardException("Unable to create directory: " + directory.getAbsolutePath());
            }
        }

    }

}
