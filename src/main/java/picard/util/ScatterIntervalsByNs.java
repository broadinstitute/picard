package picard.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.*;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.argumentcollections.ReferenceArgumentCollection;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.File;
import java.util.*;


/**
 * A Tool for breaking up a reference into intervals of alternating regions of N and ACGT bases.
 *
 * <br/>
 * <br/>
 * <h3>Summary</h3>
 * Used for creating a broken-up interval list that can be used for scattering a variant-calling pipeline in a way that
 * will not cause problems at the edges of the intervals. By using large enough N blocks (so that the tools will not be
 * able to anchor on both sides) we can be assured that the results of scattering and gathering the variants with the
 * resulting interval list will be the same as calling with one large region.
 * <br/>
 * <h3>Input</h3>
 * <il>
 *     <li>A reference file to use for creating the intervals</li>
 *     <li>Which type of intervals to emit in the output (Ns only, ACGT only or both).</li>
 *     <li>An integer indicating the largest number of Ns in a contiguous block that will be "tolerated" and not
 *     converted into an N block.</li>
 * </il>
 * <br/>
 * <h3>Output</h3>
 * <br/>
 * An interval list (with a SAM header) where the names of the intervals are labeled (either N-block or ACGT-block) to
 * indicate what type of block they define.
 *
 *
 * <h3>Usage example</h3>
 * <h4>Create an interval list of intervals that do not contain any N blocks for use with haplotype caller on short reads</h4>
 * <pre>
 * java -jar picard.jar ScatterIntervalsByNs \
 *       R=reference_sequence.fasta \
 *       OT=BOTH \
 *       O=output.interval_list
 * </pre>
 *
 * @author Yossi Farjoun
 **/


@CommandLineProgramProperties(
        summary = ScatterIntervalsByNs.USAGE_SUMMARY + ScatterIntervalsByNs.USAGE_DETAILS,
        oneLineSummary = ScatterIntervalsByNs.USAGE_SUMMARY,
        programGroup = ReferenceProgramGroup.class
)
@DocumentedFeature
public class ScatterIntervalsByNs extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Writes an interval list created by splitting a reference at Ns.";
    static final String USAGE_DETAILS = "A Program for breaking up a reference into intervals of alternating regions of N and ACGT bases." +
            "<br/>" +
            "<br/>" +
            "<br/>" +
            "Used for creating a broken-up interval list that can be used for scattering a variant-calling pipeline in a way that will not cause problems at the edges of the intervals. " +
            "By using large enough N blocks (so that the tools will not be able to anchor on both sides) we can be assured that the results of scattering and gathering the variants with " +
            "the resulting interval list will be the same as calling with one large region.\n" +
            "<br/>" +
            "<h3>Input</h3>" +
            "- A reference file to use for creating the intervals (needs to have index and dictionary next to it.)\n" +
            "- Which type of intervals to emit in the output (Ns only, ACGT only or both.)\n" +
            "- An integer indicating the largest number of Ns in a contiguous block that will be \"tolerated\" and not converted into an N block.\n" +
            "\n" +
            "<h3>Output</h3>" +
            "- An interval list (with a SAM header) where the names of the intervals are labeled (either N-block or ACGT-block) to indicate what type of block they define.\n" +
            "\n" +
            "<h3>Usage example</h3>" +
            "<h4>Create an interval list of intervals that do not contain any N blocks for use with haplotype caller on short reads</h4>" +
            "<pre>" +
            "java -jar picard.jar ScatterIntervalsByNs \\\n" +
            "      REFERENCE=reference_sequence.fasta \\\n" +
            "      OUTPUT_TYPE=ACGT \\\n" +
            "      OUTPUT=output.interval_list\n" +
            "</pre>\n" +
            "\n";
    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output file for interval list.")
    public File OUTPUT;

    @Argument(shortName = "OT", doc = "Type of intervals to output.", optional = true)
    public OutputType OUTPUT_TYPE = OutputType.BOTH;

    @Argument(shortName = "N", doc = "Maximal number of contiguous N bases to tolerate, thereby continuing the current ACGT interval.", optional = true)
    public int MAX_TO_MERGE = 1;

    //not using an enum since Interval.name is a String, and am using that to define the type of the Interval
    private static final String
            ACGTmer = "ACGTmer",
            Nmer    = "Nmer";

    //an enum to determine which types of intervals get outputted
    private enum OutputType {
        N(Nmer),
        ACGT(ACGTmer),
        BOTH(Nmer, ACGTmer);

        private final Set<String> acceptedTypes;

        public Boolean accepts(final String string) {return acceptedTypes.contains(string);}

        OutputType(final String... strings) {
            acceptedTypes = new HashSet<>();
            Collections.addAll(acceptedTypes, strings);
        }
    }

    private static final Log log = Log.getInstance(ScatterIntervalsByNs.class);
    private static final ProgressLogger locusProgress = new ProgressLogger(log, (int) 1e7, "examined", "loci");
    private static final ProgressLogger intervalProgress = new ProgressLogger(log, (int) 10, "found", "intervals");

    // return a custom argument collection since this tool uses a (required) argument name
    // of "REFERENCE", not "REFERENCE_SEQUENCE"
    @Override
    protected ReferenceArgumentCollection makeReferenceArgumentCollection() {
        return new ScatterIntervalsByNReferenceArgumentCollection();
    }

    public static class ScatterIntervalsByNReferenceArgumentCollection implements ReferenceArgumentCollection {
        @Argument(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence to use. " +
                "Note: this tool requires that the reference fasta " +
                "has both an associated index and a dictionary.")
        public File REFERENCE;

        @Override
        public File getReferenceFile() {
            return REFERENCE;
        };
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
        IOUtil.assertFileIsWritable(OUTPUT);

        final ReferenceSequenceFile refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE, true);
        if (!refFile.isIndexed()) {
            throw new IllegalStateException("Reference file must be indexed, but no index file was found");
        }
        if (refFile.getSequenceDictionary() == null) {
            throw new IllegalStateException("Reference file must include a dictionary, but no dictionary file was found");
        }

        // get the intervals
        final IntervalList intervals = segregateReference(refFile, MAX_TO_MERGE);

        log.info(String.format("Found %d intervals in %d loci during %s seconds", intervalProgress.getCount(), locusProgress.getCount(), locusProgress.getElapsedSeconds()));

        /**********************************
         * Now output regions for calling *
         **********************************/

        final IntervalList outputIntervals = new IntervalList(intervals.getHeader().clone());
        log.info(String.format("Collecting requested type of intervals (%s)", OUTPUT_TYPE));

        intervals.getIntervals().stream().filter(i -> OUTPUT_TYPE.accepts(i.getName())).forEach(outputIntervals::add);

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
    static IntervalList segregateReference(final ReferenceSequenceFile refFile, final int maxNmerToMerge) {
        final List<Interval> preliminaryIntervals = new LinkedList<>();
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(refFile.getSequenceDictionary());
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final IntervalList finalIntervals = new IntervalList(header);

        //iterate over all the sequences in the dictionary
        for (final SAMSequenceRecord rec : refFile.getSequenceDictionary().getSequences()) {
            final ReferenceSequence ref = refFile.getSequence(rec.getSequenceName());
            final byte[] bytes = ref.getBases();
            StringUtil.toUpperCase(bytes);

            boolean nBlockIsOpen = SequenceUtil.isNoCall(bytes[0]);
            int start = 0;

            for (int i = 0; i < bytes.length; ++i) {
                locusProgress.record(rec.getSequenceName(), i);
                final boolean currentBaseIsN = SequenceUtil.isNoCall(bytes[i]);

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
                        preliminaryIntervals.get(0).getContig(),
                        preliminaryIntervals.get(0).getStart(),
                        preliminaryIntervals.get(2).getEnd(), false, ACGTmer);

                //remove the first 3 elements of the list
                for (int i = 0; i < 3; ++i) {
                    preliminaryIntervals.remove(0);
                }
                //and replace them with the newly created one
                preliminaryIntervals.add(0, temp);
            } else { //if cannot merge top three intervals, transfer the top intervals to finalIntervals
                final Interval remove = preliminaryIntervals.remove(0);
                finalIntervals.add(remove);
                intervalProgress.record(remove.getContig(),remove.getStart());
            }
        }
        return finalIntervals;
    }
}
