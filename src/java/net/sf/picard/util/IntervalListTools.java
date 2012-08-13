package net.sf.picard.util;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.util.SequenceUtil;

import java.io.File;
import java.util.*;

import net.sf.picard.cmdline.CommandLineProgram;

/**
 * Little class to aid working with interval lists.
 *
 * @author Tim Fennell
 */
public class IntervalListTools extends CommandLineProgram {
    @Usage public final String USAGE = getStandardUsagePreamble() + " General tool for manipulating interval lists, " +
            "including sorting, merging, padding, uniqueifying. Default operation if given one or more inputs is to " +
            "merge and sort them.  Other options are controlled by arguments.";

    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc="One or more interval lists. If multiple interval lists are provided the output is the" +
                "result of merging the inputs.")
    public List<File> INPUT;

    @Option(doc="The output interval list file to write (if SCATTER_COUNT is 1) or the directory into which " +
            "to write the scattered interval sub-directories", shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional=true)
    public File OUTPUT;

    @Option(doc="The amount to pad each end of the intervals by before other operations are undertaken. Negative numbers are allowed " +
            "and indicate intervals should be shrunk. Resulting intervals < 0 bases long will be removed.", optional=true)
    public int PADDING = 0;

    @Option(doc="If true, merge overlapping and adjacent intervals to create a list of unique intervals. Implies SORT=true")
    public boolean UNIQUE = false;

    @Option(doc="If true, sort the resulting interval list by coordinate.")
    public boolean SORT = true;

    @Option(doc="One or more lines of comment to add to the header of the output file.", optional=true)
    public List<String> COMMENT = null;

    @Option(doc="The number of files into which to scatter the resulting list by locus,")
    public int SCATTER_COUNT = 1;

    private final Log log = Log.getInstance(IntervalListTools.class);

    // Stock main method
    public static void main(final String[] args) {
        new IntervalListTools().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        // Check inputs
        for (final File f : INPUT) IoUtil.assertFileIsReadable(f);
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

        if (UNIQUE && !SORT ) {
            log.warn("UNIQUE=true requires sorting but SORT=false was specified.  Sorting anyway!");
            SORT = true;
        }

        // Ensure that all the sequence dictionaries agree and merge the lists
        IntervalList merged= null;
        for (final IntervalList in : lists) {
            if (merged == null) {
                merged = in;
            }
            else {
                SequenceUtil.assertSequenceDictionariesEqual(merged.getHeader().getSequenceDictionary(),
                                                             in.getHeader().getSequenceDictionary());

                for (final Interval i : in) {
                    merged.add(i);
                }
            }
        }

        if (SORT) merged.sort();
        final List<Interval> finalIntervals = UNIQUE ? merged.getUniqueIntervals() : merged.getIntervals();

        // Decide on a PG ID and make a program group
        final SAMFileHeader header = merged.getHeader();
        final Set<String> pgs = new HashSet<String>();
        for (final SAMProgramRecord pg : header.getProgramRecords()) pgs.add(pg.getId());
        for (int i=1; i<Integer.MAX_VALUE; ++i) {
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
        final long idealSplitLength = Math.max((long)Math.floor(list.getUniqueBaseCount() / (1.0*SCATTER_COUNT)), 1);
        int splitLength = 0;
        IntervalList split = new IntervalList(list.getHeader());
        int index = 1;   // The index of the next split file to write
        final Iterator<Interval> it = list.iterator();
        int totalIntervals = 0;

        while (it.hasNext() && index < SCATTER_COUNT) {
            final Interval interval = it.next();
            final int projectedSize = splitLength + interval.length();
            if (projectedSize < idealSplitLength) {
                split.add(interval);
                totalIntervals++;
                splitLength += interval.length();
            }
            else if (projectedSize == idealSplitLength) {
                split.add(interval);
                totalIntervals++;
                split.write(createDirectoryAndGetScatterFile(index++));
                split = new IntervalList(list.getHeader());
                splitLength = 0;
            }
            else {
                final int diff = (int)(idealSplitLength - splitLength);
                // Make one interval to get to the right amount of territory,
                // add it to the current split and write it out
                final Interval firstHalf = new Interval(interval.getSequence(), interval.getStart(),
                        interval.getStart()+diff-1, interval.isNegativeStrand(), interval.getName());
                split.add(firstHalf);
                totalIntervals++;
                split.write(createDirectoryAndGetScatterFile(index++));
                // Add the remainder to the next split
                split = new IntervalList(list.getHeader());
                final Interval secondHalf = new Interval(interval.getSequence(), interval.getStart()+diff,
                        interval.getEnd(), interval.isNegativeStrand(), interval.getName());
                split.add(secondHalf);
                totalIntervals++;
                splitLength = secondHalf.length();
            }
        }
        // Write everything left to the last split
        while (it.hasNext()) {
            split.add(it.next());
            totalIntervals++;
        }
        split.write(createDirectoryAndGetScatterFile(index));
        return totalIntervals;


    }

    public static File getScatteredFileName(final File scatterDirectory, final int scatterTotal, final int index) {
        return new File(scatterDirectory.getAbsolutePath() + "/temp_" + index + "_of_" +
                scatterTotal + "/scattered.intervals");

    }

    private File createDirectoryAndGetScatterFile(final int index) {
        createDirectoryOrFail(OUTPUT);
        final File result = getScatteredFileName(OUTPUT, SCATTER_COUNT, index);
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
