package picard.util;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.programgroups.Intervals;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.StringUtil;

import java.io.File;
import java.lang.Boolean;import java.lang.Override;import java.lang.String;import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;


/**
 * A CLP for breaking up a reference into intervals of Ns and ACGTs bases.
 * Used for creating a broken-up interval list for calling WGS.
 *
 * @author Yossi Farjoun
 */

@CommandLineProgramProperties(
        usage = ScatterIntervalsByNs.USAGE_SUMMARY + ScatterIntervalsByNs.USAGE_DETAILS,
        usageShort = ScatterIntervalsByNs.USAGE_SUMMARY,
        programGroup = Intervals.class
)
public class ScatterIntervalsByNs extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Writes an interval list based on splitting the reference by Ns.  ";
    static final String USAGE_DETAILS = "This tool identifies positions in the reference where the basecalls are Ns and writes out an " +
            "interval list using the resulting coordinates (excluding the N bases). This can be used to create an interval list for " +
            "whole genome sequence (WGS) for e.g. scatter-gather purposes, as an alternative to using fixed-length intervals. The number " +
            "of contiguous Ns that can be tolerated before creating a break is adjustable from the command line.<br />" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar ScatterIntervalsByNs \\<br />" +
            "      R=reference_sequence.fasta \\<br />" +
            "      OT=BOTH \\<br />" +
            "      O=output.interval_list" +
            "</pre>" +
            "<hr />";
    @Option(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence to use.")
    public File REFERENCE;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output file for interval list.")
    public File OUTPUT;

    @Option(shortName = "OT", doc = "Type of intervals to output.", optional = true)
    public OutputType OUTPUT_TYPE = OutputType.BOTH;

    @Option(shortName = "N", doc = "Maximal number of contiguous N bases to tolerate, thereby continuing the current ACGT interval.", optional = true)
    public int MAX_TO_MERGE = 1;

    //not using an enum since Interval.name is a String, and am using that to define the type of the Interval
    static final String
            ACGTmer = "ACGTmer",
            Nmer    = "Nmer";

    //an enum to determine which types of intervals get outputted
    private enum OutputType {
        N(Nmer),
        ACGT(ACGTmer),
        BOTH(Nmer, ACGTmer);

        private final Set acceptedTypes;

        public Boolean accepts(final String string) {return acceptedTypes.contains(string);}

        OutputType(final String... strings) {
            acceptedTypes = new HashSet<String>();
            Collections.addAll(acceptedTypes, strings);
        }
    }

    private static final Log log = Log.getInstance(ScatterIntervalsByNs.class);
    final ProgressLogger locusProgress = new ProgressLogger(log, (int) 1e7, "examined", "loci");
    final ProgressLogger intervalProgress = new ProgressLogger(log, (int) 10, "found", "intervals");

    public static void main(final String[] args) {
        new ScatterIntervalsByNs().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(REFERENCE);
        IOUtil.assertFileIsWritable(OUTPUT);

        final ReferenceSequenceFile refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE, true);

        // get the intervals
        final IntervalList intervals = segregateReference(refFile, MAX_TO_MERGE);

        log.info(String.format("Found %d intervals in %d loci during %s seconds", intervalProgress.getCount(), locusProgress.getCount(), locusProgress.getElapsedSeconds()));

        /**********************************
         * Now output regions for calling *
         **********************************/

        final IntervalList outputIntervals = new IntervalList(intervals.getHeader().clone());
        log.info(String.format("Collecting requested type of intervals (%s)", OUTPUT_TYPE));

        for (final Interval i : intervals.getIntervals()) {
            if (OUTPUT_TYPE.accepts(i.getName())) {
                outputIntervals.add(i);
            }
        }

        log.info("Writing Intervals.");
        outputIntervals.write(OUTPUT);

        log.info(String.format("Execution ending. Total time %d seconds", locusProgress.getElapsedSeconds()));

        return 0;
    }

    /**
     * ****************************************************************
     * Generate an interval list that alternates between Ns and ACGTs *
     * ****************************************************************
     */
    public static IntervalList segregateReference(final ReferenceSequenceFile refFile, final int maxNmerToMerge) {
        final List<Interval> preliminaryIntervals = new LinkedList<Interval>();
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(refFile.getSequenceDictionary());
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final IntervalList finalIntervals = new IntervalList(header);

        //iterate over all the sequences in the dictionary
        for (final SAMSequenceRecord rec : refFile.getSequenceDictionary().getSequences()) {
            final ReferenceSequence ref = refFile.getSequence(rec.getSequenceName());
            final byte[] bytes = ref.getBases();
            StringUtil.toUpperCase(bytes);

            boolean nBlockIsOpen = (bytes[0] == 'N');
            int start = 0;

            for (int i = 0; i < bytes.length; ++i) {
                final boolean currentBaseIsN = (bytes[i] == 'N');

                //create intervals when switching, i.e "nBlockIsOpen" disagrees with "currentBaseIsN"
                if (nBlockIsOpen != currentBaseIsN) {
                    preliminaryIntervals.add(new Interval(rec.getSequenceName(), start + 1, i, false, nBlockIsOpen ? Nmer : ACGTmer));
                    start = i;
                    nBlockIsOpen = !nBlockIsOpen;
                }
            }
            // Catch the last block of chromosome
            preliminaryIntervals.add(new Interval(rec.getSequenceName(), start + 1, bytes.length, false, nBlockIsOpen ? Nmer : ACGTmer));
        }

        // now that we have the whole list, we need to remove the short Nmers.
        // process the list, replacing trios with short Nmers in the middle with longer intervals:
        while (!preliminaryIntervals.isEmpty()) {

            //if top trio match the bill, replace them with a merged interval,
            // and push it back the top of the list (we expect alternating Nmers and ACGTmers, but
            // not assuming it in the logic)

            //(I want this to be fast and the strings are all copies of the static prototypes Nmer and ACGTmer )
            //noinspection StringEquality
            if (preliminaryIntervals.size() >= 3 && // three or more intervals
                    preliminaryIntervals.get(0).getName() == ACGTmer &&   //an N-mer
                    preliminaryIntervals.get(1).getName() == Nmer &&      //between two
                    preliminaryIntervals.get(2).getName() == ACGTmer &&   //ACGT-mers
                    preliminaryIntervals.get(0).abuts(preliminaryIntervals.get(1)) && // all abutting
                    preliminaryIntervals.get(1).abuts(preliminaryIntervals.get(2)) && // each other (there are many contigs...)
                    preliminaryIntervals.get(1).length() <= maxNmerToMerge) //and the N-mer is of length N or less
            {
                // create the new ACGTmer interval
                final Interval temp = new Interval(
                        preliminaryIntervals.get(0).getSequence(),
                        preliminaryIntervals.get(0).getStart(),
                        preliminaryIntervals.get(2).getEnd(), false, ACGTmer);

                //remove the first 3 elements of the list
                for (int i = 0; i < 3; ++i) {
                    preliminaryIntervals.remove(0);
                }
                //and replace them with the newly created one
                preliminaryIntervals.add(0, temp);
            } else { //if cannot merge top three intervals, transfer the top intervals to finalIntervals
                finalIntervals.add(preliminaryIntervals.remove(0));
            }
        }
        return finalIntervals;
    }
}
