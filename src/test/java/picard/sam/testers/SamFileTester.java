package picard.sam.testers;

import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaReferenceWriter;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Abstract class for doing basic on the fly SAM file testing.
 */
public abstract class SamFileTester extends CommandLineProgramTest {

    protected final SAMRecordSetBuilder samRecordSetBuilder;
    protected final Map<String, Boolean> duplicateFlags = new HashMap<>();
    private File outputDir;
    private File output;
    private int readNameCounter = 0;
    protected boolean noMateCigars = false;
    private boolean deleteOnExit = true;
    private final ArrayList<String> args = new ArrayList<>();

    public SamFileTester(final int readLength, final boolean deleteOnExit, final int defaultChromosomeLength, final ScoringStrategy duplicateScoringStrategy, final SAMFileHeader.SortOrder sortOrder, boolean recordsNeedSorting) {
        this.deleteOnExit = deleteOnExit;
        this.samRecordSetBuilder = new SAMRecordSetBuilder(recordsNeedSorting, sortOrder, true, defaultChromosomeLength, duplicateScoringStrategy);
        samRecordSetBuilder.setReadLength(readLength);
        setOutputDir();
    }

    public SamFileTester(final int readLength, final boolean deleteOnExit, final int defaultChromosomeLength, final ScoringStrategy duplicateScoringStrategy, final SAMFileHeader.SortOrder sortOrder) {
        this(readLength, deleteOnExit, defaultChromosomeLength, duplicateScoringStrategy, SAMFileHeader.SortOrder.coordinate, true);
    }

    public SamFileTester(final int readLength, final boolean deleteOnExit, final int defaultChromosomeLength, final ScoringStrategy duplicateScoringStrategy) {
        this(readLength, deleteOnExit, defaultChromosomeLength, duplicateScoringStrategy, SAMFileHeader.SortOrder.coordinate);
    }

    public SamFileTester(final int readLength, final boolean deleteOnExit, final int defaultChromosomeLength) {
        this(readLength, deleteOnExit, defaultChromosomeLength, SAMRecordSetBuilder.DEFAULT_DUPLICATE_SCORING_STRATEGY);
    }

    public SamFileTester(final int readLength, final boolean deleteOnExit) {
        this(readLength, deleteOnExit, SAMRecordSetBuilder.DEFAULT_CHROMOSOME_LENGTH);
    }

    public void setHeader(final SAMFileHeader header) {
        this.samRecordSetBuilder.setHeader(header);
    }

    public void addRecord(final SAMRecord record) {
        final String key = samRecordToDuplicatesFlagsKey(record);
        Assert.assertFalse(this.duplicateFlags.containsKey(key));
        this.duplicateFlags.put(key, record.getDuplicateReadFlag());
        this.samRecordSetBuilder.addRecord(record);
    }

    public int getNumberOfRecords() {
        return this.samRecordSetBuilder.size();
    }

    public CloseableIterator<SAMRecord> getRecordIterator() {
        return this.samRecordSetBuilder.iterator();
    }

    public File getOutput() {
        return output;
    }

    public void setOutput(final File output) {
        this.output = output;
    }


    public void addArg(final String arg) {
        args.add(arg);
    }

    public ArrayList<String> getArgs() {
        return args;
    }

    public File getOutputDir() {
        return outputDir;
    }

    private void setOutputDir() {
        this.outputDir = IOUtil.createTempDir(this.getClass().getSimpleName() + ".", ".tmp");
        if (deleteOnExit) {
            outputDir.deleteOnExit();
        }
    }

    public void setNoMateCigars(final boolean value) {
        this.noMateCigars = value;
    }

    public boolean getDeleteOnExit() {
        return deleteOnExit;
    }

    protected String samRecordToDuplicatesFlagsKey(final SAMRecord record) {
        final StringBuilder nameBuilder = new StringBuilder();
        nameBuilder.append(record.getReadName());
        nameBuilder.append("-");

        if (record.getReadUnmappedFlag()) {
            nameBuilder.append("Unmapped");
        } else {
            nameBuilder.append(record.getContig())
                    .append("-")
                    .append(record.getAlignmentStart());
        }
        nameBuilder.append("-")
                .append(record.getReadPairedFlag())
                .append("-").append(record.getNotPrimaryAlignmentFlag())
                .append("-");

        if (record.getReadPairedFlag()) {
            nameBuilder.append(record.getFirstOfPairFlag())
                    .append("-")
                    .append(record.getSecondOfPairFlag());
        } else {
            nameBuilder.append("false-false");
        }
        return nameBuilder.toString();
    }

    // Below are a bunch of utility methods for adding records to the SAMRecordSetBuilder
    public void addUnmappedFragment(final int referenceSequenceIndex,
                                    final int defaultQualityScore) {
        addFragment(referenceSequenceIndex, -1, true, false, null, null, defaultQualityScore, false);
    }

    public void addUnmappedFragment(final int referenceSequenceIndex,
                                    final String qualityString) {
        addFragment(referenceSequenceIndex, -1, true, false, null, qualityString, -1, false);
    }

    public void addUnmappedPair(final int referenceSequenceIndex,
                                final int defaultQualityScore) {
        addMatePair(referenceSequenceIndex, -1, -1, true, true, false, false, null, null, false, false, false, false, false, defaultQualityScore);
    }

    public void addMappedFragment(final int referenceSequenceIndex, final int alignmentStart, final boolean isDuplicate,
                                  final int defaultQualityScore) {
        addFragment(referenceSequenceIndex, alignmentStart, false, isDuplicate, null, null, defaultQualityScore, false);
    }

    public void addMappedFragment(final int referenceSequenceIndex, final int alignmentStart, final boolean isDuplicate,
                                  final int defaultQualityScore, final boolean isSecondary) {
        addFragment(referenceSequenceIndex, alignmentStart, false, isDuplicate, null, null, defaultQualityScore, isSecondary);
    }

    public void addMappedFragment(final int referenceSequenceIndex, final int alignmentStart, final boolean isDuplicate, final String cigar,
                                  final int defaultQualityScore) {
        addFragment(referenceSequenceIndex, alignmentStart, false, isDuplicate, cigar, null, defaultQualityScore, false);
    }

    public void addMappedFragment(final int referenceSequenceIndex, final int alignmentStart, final boolean isDuplicate, final String cigar,
                                  final String qualityString,
                                  final int defaultQualityScore) {
        addFragment(referenceSequenceIndex, alignmentStart, false, isDuplicate, cigar, qualityString, defaultQualityScore, false);
    }

    public void addMappedFragment(final String readName, final int referenceSequenceIndex, final int alignmentStart, final boolean isDuplicate, final String cigar,
                                  final String qualityString, final boolean isSecondary, final boolean isSupplementary,
                                  final int defaultQualityScore) {
        addFragment(readName, referenceSequenceIndex, alignmentStart, false, isDuplicate, cigar, qualityString, defaultQualityScore, isSecondary, isSupplementary);
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
                strand1, strand2, firstOnly, false, false, defaultQualityScore);
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
                            final boolean record1NonPrimary,
                            final boolean record2NonPrimary,
                            final int defaultQualityScore) {
        addMatePair("READ" + readNameCounter++, referenceSequenceIndex, alignmentStart1, alignmentStart2, record1Unmapped, record2Unmapped,
                isDuplicate1, isDuplicate2, cigar1, cigar2, strand1, strand2, firstOnly, record1NonPrimary, record2NonPrimary,
                defaultQualityScore);
    }

    private void addFragment(final int referenceSequenceIndex, final int alignmentStart, final boolean recordUnmapped, final boolean isDuplicate, final String cigar,
                             final String qualityString, final int defaultQualityScore, final boolean isSecondary) {
        addFragment("READ" + readNameCounter++, referenceSequenceIndex, alignmentStart, recordUnmapped, isDuplicate, cigar,
                qualityString, defaultQualityScore, isSecondary, false);
    }

    private void addFragment(final String readName, final int referenceSequenceIndex, final int alignmentStart, final boolean recordUnmapped, final boolean isDuplicate, final String cigar,
                             final String qualityString, final int defaultQualityScore, final boolean isSecondary, final boolean isSupplementary) {

        final SAMRecord record = samRecordSetBuilder.addFrag(readName, referenceSequenceIndex, alignmentStart, false,
                recordUnmapped, cigar, qualityString, defaultQualityScore, isSecondary, isSupplementary);

        final String key = samRecordToDuplicatesFlagsKey(record);
        Assert.assertFalse(this.duplicateFlags.containsKey(key));
        this.duplicateFlags.put(key, isDuplicate);
    }

    public void addMatePair(final String readName,
                            final int referenceSequenceIndex1,
                            final int referenceSequenceIndex2,
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
                            final boolean record1NonPrimary,
                            final boolean record2NonPrimary,
                            final int defaultQuality,
                            final String umi) {
        final List<SAMRecord> samRecordList = samRecordSetBuilder.addPair(readName, referenceSequenceIndex1, referenceSequenceIndex2, alignmentStart1, alignmentStart2,
                record1Unmapped, record2Unmapped, cigar1, cigar2, strand1, strand2, record1NonPrimary, record2NonPrimary, defaultQuality);

        final SAMRecord record1 = samRecordList.get(0);
        final SAMRecord record2 = samRecordList.get(1);

        if (this.noMateCigars) {
            record1.setAttribute("MC", null);
            record2.setAttribute("MC", null);
        }

        if (umi != null) {
            record1.setAttribute("RX", umi);
            record2.setAttribute("RX", umi);
        }

        if (firstOnly) {
            samRecordSetBuilder.getRecords().remove(record2);
        }

        final String key1 = samRecordToDuplicatesFlagsKey(record1);
        Assert.assertFalse(this.duplicateFlags.containsKey(key1));
        this.duplicateFlags.put(key1, isDuplicate1);

        final String key2 = samRecordToDuplicatesFlagsKey(record2);
        Assert.assertFalse(this.duplicateFlags.containsKey(key2));
        this.duplicateFlags.put(key2, isDuplicate2);
    }

    public void addMatePair(final String readName,
                            final int referenceSequenceIndex1,
                            final int referenceSequenceIndex2,
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
                            final boolean record1NonPrimary,
                            final boolean record2NonPrimary,
                            final int defaultQuality) {
        addMatePair(readName, referenceSequenceIndex1, referenceSequenceIndex2, alignmentStart1, alignmentStart2,
                record1Unmapped, record2Unmapped, isDuplicate1, isDuplicate2, cigar1, cigar2, strand1, strand2,
                firstOnly, record1NonPrimary, record2NonPrimary, defaultQuality, null);
    }
    public void addMatePair(final String readName,
                            final int referenceSequenceIndex,
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
                            final boolean record1NonPrimary,
                            final boolean record2NonPrimary,
                            final int defaultQuality) {
        addMatePair(readName, referenceSequenceIndex, referenceSequenceIndex, alignmentStart1, alignmentStart2, record1Unmapped, record2Unmapped,
                isDuplicate1, isDuplicate2, cigar1, cigar2, strand1, strand2, firstOnly, record1NonPrimary, record2NonPrimary, defaultQuality);
    }

    protected abstract void test() throws IOException;

    /**
     * Sets up the basic command line arguments for input and output and runs instanceMain.
     */
    public void runTest() {
        runTest(".sam");
    }

    /**
     * Sets up the basic command line arguments for input and output and runs instanceMain.
     */
    public void runTest(final String inputExtension) {
        try {
            final File input = createInputFile(inputExtension);

            output = new File(outputDir, "output.sam");
            args.add("INPUT=" + input.getAbsoluteFile());
            args.add("OUTPUT=" + output.getAbsoluteFile());
            if(inputExtension.equals(".cram")){
                args.add("REFERENCE_SEQUENCE=input.fasta");
            }
            Assert.assertEquals(runPicardCommandLine(args), 0);
            test();
        } catch (IOException ex) {
            // Rethrows as unchecked exception to avoid having to change existing tests
            throw new RuntimeException(ex);
        }
    }

    private File createInputFile(final String extension) throws IOException {
        // Create the input file
        final File input = new File(outputDir, "input" + extension);

        final SAMFileWriterFactory samFileWriterFactory = new SAMFileWriterFactory();

        final SAMFileWriter writer;
        if(extension.equals(".cram")) {
            final Path fasta = IOUtil.getPath("Testinput.fasta");
            IOUtil.deleteOnExit(fasta);
            IOUtil.deleteOnExit(ReferenceSequenceFileFactory.getFastaIndexFileName(fasta));
            IOUtil.deleteOnExit(ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(fasta));

            samRecordSetBuilder.writeRandomReference(fasta);
            writer = samFileWriterFactory.makeWriter(samRecordSetBuilder.getHeader(), true, input, fasta.toFile());
        } else {
            writer = samFileWriterFactory.makeWriter(samRecordSetBuilder.getHeader(), true, input, null);
        }

        samRecordSetBuilder.getRecords().forEach(writer::addAlignment);
        writer.close();
        return input;
    }

    public SamReader getInput() {
        return samRecordSetBuilder.getSamReader();
    }

    public SAMRecordSetBuilder getSamRecordSetBuilder() { return samRecordSetBuilder; }

}
