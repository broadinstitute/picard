package picard.sam.markduplicates;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;
import org.apache.commons.io.output.NullOutputStream;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;

@CommandLineProgramProperties(
        summary = CheckDuplicateMarking.USAGE_DETAILS,
        oneLineSummary = CheckDuplicateMarking.USAGE_SUMMARY,
        programGroup = DiagnosticsAndQCProgramGroup.class)

@DocumentedFeature(enable = false)
public class CheckDuplicateMarking extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Checks the consistency of duplicate markings.";
    static final String USAGE_DETAILS = "This tool checks that all reads with the same queryname have their duplicate marking flags set the same way. " +
            "NOTE: This tool does NOT check that the duplicate marking is correct. The ONLY thing that it checks is that the 0x400 bit-flags of records " +
            "with the same queryname are equal.";

    @Argument(doc = "Input BAM or SAM file to check.", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "Output file into which bad querynames will be placed (if not null).", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional = true)
    public File OUTPUT = null;

    @Argument(doc = "Which reads of the same name should be checked to have same duplicate marking.")
    public Mode MODE = Mode.ALL;

    private String currentReadName = "";
    private boolean currentReadDuplicateMarked = false;

    private static final Log log = Log.getInstance(CheckDuplicateMarking.class);
    private final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Checked.");

    private static final int NUM_WARNINGS = 100;
    private int numBadRecords = 0;

    public enum Mode implements CommandLineParser.ClpEnum {
        ALL("Check all reads."),
        PRIMARY_ONLY("Check primary alignments."),
        PRIMARY_MAPPED_ONLY("Check mapped alignments."),
        PRIMARY_PROPER_PAIR_ONLY("Check mapped alignments.");

        private final String message;

        Mode(final String message) {
            this.message = message;
        }

        public String getHelpDoc() {
            return message;
        }
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        if (OUTPUT != null) {
            IOUtil.assertFileIsWritable(OUTPUT);
        }

        try (SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT)) {

            checkDuplicateMarkingsInIterator(getSortedRecordsFromReader(reader));
            if (numBadRecords > 0) {
                log.error("Found " + numBadRecords + " records that do not agree on their duplicate flag.");
            } else {
                log.info("All records' duplicate markings agree.");
            }

        } catch (IOException e) {
            throw new PicardException("Error while reading input file " + INPUT, e);
        }

        return numBadRecords > 0 ? 1 : 0;
    }

    private Iterator<SAMRecord> getSortedRecordsFromReader(final SamReader reader) {
        if (reader.getFileHeader().getSortOrder() == SAMFileHeader.SortOrder.queryname) {
            return reader.iterator();
        }

        log.info("Input file isn't queryname sorted. Sorting into temp space.");

        final ProgressLogger sortProgress = new ProgressLogger(log, (int) 1e6, "Read into sorter");
        final SortingCollection<SAMRecord> alignmentSorter = SortingCollection.newInstance(SAMRecord.class,
                new BAMRecordCodec(reader.getFileHeader()),
                new SAMRecordQueryNameComparator(),
                MAX_RECORDS_IN_RAM);

        for (final SAMRecord rec : reader) {
            alignmentSorter.add(rec);
            sortProgress.record(rec);
        }
        final CloseableIterator<SAMRecord> iterator = alignmentSorter.iterator();

        alignmentSorter.cleanup();

        return iterator;
    }

    private boolean checkAndTallyRecordDuplicateMarking(final SAMRecord rec) {
        if (!rec.getReadName().equals(currentReadName)) {
            // this case the queryname changed, and thus there's no comparison to make
            currentReadName = rec.getReadName();
            currentReadDuplicateMarked = rec.getDuplicateReadFlag();
        } else if (rec.getDuplicateReadFlag() != currentReadDuplicateMarked) {
            // Here the current queryname is the same, but the duplicate flag doesn't match the first record with that queryname
            numBadRecords++;

            if (numBadRecords <= NUM_WARNINGS) {
                log.warn(() -> "Reads with queryname " + currentReadName + " have different duplicate flags (at " +
                        rec.getContig() + ":" + rec.getStart() + ")");
            }

            if (numBadRecords == NUM_WARNINGS) {
                log.warn("Further warnings will be suppressed.");
            }

            return false;
        }
        return true;
    }

    private void checkDuplicateMarkingsInIterator(final Iterator<SAMRecord> iterator) throws IOException {
        try (PrintWriter writer = OUTPUT == null ?
                new PrintWriter(NullOutputStream.NULL_OUTPUT_STREAM) :
                new PrintWriter(new FileWriter(OUTPUT))) {

            while (iterator.hasNext()) {
                final SAMRecord rec = iterator.next();

                if (MODE != Mode.ALL && rec.isSecondaryOrSupplementary()) {
                    continue;
                }

                if (MODE == Mode.PRIMARY_MAPPED_ONLY && rec.getReadUnmappedFlag()) {
                    continue;
                }

                if (MODE == Mode.PRIMARY_PROPER_PAIR_ONLY && !rec.getProperPairFlag()) {
                    continue;
                }

                if (!checkAndTallyRecordDuplicateMarking(rec)) {
                    writer.println(rec.getReadName());
                }
                progress.record(rec);
            }
        }
    }
}
