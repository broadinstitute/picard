package picard.sam.testers;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import picard.cmdline.CommandLineProgram;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Abstract class for doing basic on the fly SAM file testing.
 */
public abstract class SamFileTester {

    public static final String TEST_DATA_BASE_DIR = "testdata/picard/sam/";
    private final SAMRecordSetBuilder samRecordSetBuilder;
    protected final Map<String, Boolean> duplicateFlags = new HashMap<String, Boolean>();
    private File outputDir;
    private File output;
    private int readNameCounter = 0;
    private boolean noMateCigars = false;
    private boolean deleteOnExit = true;
    private final ArrayList<String> args = new ArrayList<String>();

    public SamFileTester(final int readLength, final boolean deleteOnExit, final int defaultChromosomeLength) {
        this.deleteOnExit = deleteOnExit;
        this.samRecordSetBuilder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate, true, defaultChromosomeLength);
        samRecordSetBuilder.setReadLength(readLength);
        setOutputDir();
    }

    public SamFileTester(final int readLength, final boolean deleteOnExit) {
        this.deleteOnExit = deleteOnExit;
        this.samRecordSetBuilder = new SAMRecordSetBuilder();
        samRecordSetBuilder.setReadLength(readLength);
        setOutputDir();
    }

    public File getOutput() {
        return output;
    }

    public void addArg(final String arg) {
        args.add(arg);
    }

    public File getOutputDir() {
        return outputDir;
    }

    private void setOutputDir() {
        this.outputDir = IOUtil.createTempDir(this.getClass().getSimpleName() + ".", ".tmp");
        if(deleteOnExit){
            outputDir.deleteOnExit();
        }
    }

    public void setNoMateCigars(final boolean value) {
        this.noMateCigars = value;
    }

    protected String samRecordToDuplicatesFlagsKey(final SAMRecord record) {
        String readName = record.getReadName()
                + "-"
                + record.getReadPairedFlag()
                + "-";
        if (record.getReadPairedFlag()) {
            readName += record.getFirstOfPairFlag()
                    + "-"
                    + record.getSecondOfPairFlag();
        } else {
            readName += "false-false";
        }
        return readName;
    }
    // Below are a bunch of utility methods for adding records to the SAMRecordSetBuilder
    public void addUnmappedFragment(final int referenceSequenceIndex,
                                    final int defaultQualityScore) {
        addFragment(referenceSequenceIndex, -1, true, false, null, null, defaultQualityScore);
    }

    public void addUnmappedFragment(final int referenceSequenceIndex,
                                    final String qualityString){
        addFragment(referenceSequenceIndex, -1, true, false, null, qualityString, -1);
    }

    public void addUnmappedPair(final int referenceSequenceIndex,
                                final int defaultQualityScore) {
        addMatePair(referenceSequenceIndex, -1, -1, true, true, false, false, null, null, false, false, false, defaultQualityScore);
    }

    public void addMappedFragment(final int referenceSequenceIndex, final int alignmentStart, final boolean isDuplicate,
                                  final int defaultQualityScore) {
        addFragment(referenceSequenceIndex, alignmentStart, false, isDuplicate, null, null, defaultQualityScore);
    }

    public void addMappedFragment(final int referenceSequenceIndex, final int alignmentStart, final boolean isDuplicate, final String cigar,
                                  final int defaultQualityScore) {
        addFragment(referenceSequenceIndex, alignmentStart, false, isDuplicate, cigar, null, defaultQualityScore);
    }

    public void addMappedFragment(final int referenceSequenceIndex, final int alignmentStart, final boolean isDuplicate, final String cigar,
                                  final String qualityString,
                                  final int defaultQualityScore) {
        addFragment(referenceSequenceIndex, alignmentStart, false, isDuplicate, cigar, qualityString, defaultQualityScore);
    }

    public void addMappedPair(final int referenceSequenceIndex,
                              final int alignmentStart1,
                              final int alignmentStart2,
                              final boolean isDuplicate1,
                              final boolean isDuplicate2,
                              final int defaultQualityScore) {
        addMappedPair(referenceSequenceIndex, alignmentStart1, alignmentStart2, isDuplicate1, isDuplicate2, null, null,
                false, defaultQualityScore);
    }

    public void addMappedPair(final int referenceSequenceIndex,
                              final int alignmentStart1,
                              final int alignmentStart2,
                              final boolean isDuplicate1,
                              final boolean isDuplicate2,
                              final String cigar1,
                              final String cigar2,
                              final boolean firstOnly,
                              final int defaultQualityScore) {
        addMappedPair(referenceSequenceIndex, alignmentStart1, alignmentStart2, isDuplicate1, isDuplicate2, cigar1,
                cigar2, false, true, firstOnly, defaultQualityScore);
    }

    public void addMappedPair(final int referenceSequenceIndex,
                              final int alignmentStart1,
                              final int alignmentStart2,
                              final boolean isDuplicate1,
                              final boolean isDuplicate2,
                              final String cigar1,
                              final String cigar2,
                              final boolean strand1,
                              final boolean strand2,
                              final boolean firstOnly,
                              final int defaultQualityScore) {
        addMatePair(referenceSequenceIndex, alignmentStart1, alignmentStart2, false, false, isDuplicate1, isDuplicate2, cigar1, cigar2,
                strand1, strand2, firstOnly, defaultQualityScore);
    }

    private void addFragment(final int referenceSequenceIndex, final int alignmentStart, final boolean recordUnmapped, final boolean isDuplicate, final String cigar,
                             final String qualityString, final int defaultQualityScore) {
        final SAMRecord record = samRecordSetBuilder.addFrag("READ" + readNameCounter++, referenceSequenceIndex, alignmentStart, false,
                recordUnmapped, cigar, qualityString, defaultQualityScore);

        this.duplicateFlags.put(samRecordToDuplicatesFlagsKey(record), isDuplicate);
    }

    public void addMatePair(final int referenceSequenceIndex,
                            final int alignmentStart1,
                            final int alignmentStart2,
                            final boolean record1Unmapped,
                            final boolean record2Unmapped,
                            final boolean isDuplicate1,
                            final boolean isDuplicate2,
                            final String cigar1,
                            final String cigar2,
                            final boolean strand1,
                            final boolean strand2,
                            final boolean firstOnly,
                            final int defaultQuality) {
        final List<SAMRecord> samRecordList = samRecordSetBuilder.addPair("READ" + readNameCounter++, referenceSequenceIndex, alignmentStart1, alignmentStart2,
                record1Unmapped, record2Unmapped, cigar1, cigar2, strand1, strand2, defaultQuality);

        final SAMRecord record1 = samRecordList.get(0);
        final SAMRecord record2 = samRecordList.get(1);

        if (this.noMateCigars) {
            record1.setAttribute("MC", null);
            record2.setAttribute("MC", null);
        }

        if (firstOnly) {
            samRecordSetBuilder.getRecords().remove(record2);
        }

        this.duplicateFlags.put(samRecordToDuplicatesFlagsKey(record1), isDuplicate1);
        this.duplicateFlags.put(samRecordToDuplicatesFlagsKey(record2), isDuplicate2);
    }

    protected abstract void test();

    protected abstract CommandLineProgram getProgram();

    /**
     * Sets up the basic command line arguments for input and output and runs instanceMain.
     */
    public void runTest() {
        if (getProgram() != null) {
            final File input = createInputFile();

            output = new File(outputDir, "output.sam");
            args.add("INPUT=" + input.getAbsoluteFile());
            args.add("OUTPUT=" + output.getAbsoluteFile());
            Assert.assertEquals(getProgram().instanceMain(args.toArray(new String[args.size()])), 0);
        }
        test();
    }

    private File createInputFile() {
        // Create the input file
        final File input = new File(outputDir, "input.sam");
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(samRecordSetBuilder.getHeader(), true, input);
        for (final SAMRecord record : samRecordSetBuilder.getRecords()) {
            writer.addAlignment(record);
        }
        writer.close();
        return input;
    }

    public SAMFileReader getInput(){
        return samRecordSetBuilder.getSamReader();
    }
}