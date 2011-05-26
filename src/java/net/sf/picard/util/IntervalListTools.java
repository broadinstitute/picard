package net.sf.picard.util;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.util.SequenceUtil;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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

    @Option(doc="The output interval list file to write", shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional=true)
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

    private final Log log = Log.getInstance(IntervalListTools.class);

    // Stock main method
    public static void main(final String[] args) {
        new IntervalListTools().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        // Check inputs
        for (final File f : INPUT) IoUtil.assertFileIsReadable(f);
        if (OUTPUT != null) IoUtil.assertFileIsWritable(OUTPUT);

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
        for (final Interval i : finalIntervals) {
            output.add(i);
        }

        if (OUTPUT != null) output.write(OUTPUT);

        long total = 0;
        for (final Interval i : output) {
            total += i.length();
        }

        log.info("Output " + output.size() + " intervals totalling " + total + " bases.");

        return 0;
    }
}
