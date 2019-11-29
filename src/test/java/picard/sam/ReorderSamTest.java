package picard.sam;

import htsjdk.samtools.*;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.sam.testers.ValidateSamTester;
import picard.util.SequenceDictionaryUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;


public class ReorderSamTest extends CommandLineProgramTest {

    @Override
    public String getCommandLineProgramName() {
        return ReorderSam.class.getSimpleName();
    }

    @DataProvider(name = "testDictionaryData")
    public Iterator<Object[]> testDictionaryData() throws IOException {
        List<Object[]> tests = new ArrayList<>();

        final SAMRecordSetBuilder setBuilder = new SAMRecordSetBuilder();

        setBuilder.addPair("pair1", setBuilder.getHeader().getSequenceIndex("chr1"), 1, 100);
        setBuilder.addPair("pair2", setBuilder.getHeader().getSequenceIndex("chr2"), 1, 100);
        setBuilder.addPair("pair3", setBuilder.getHeader().getSequenceIndex("chr3"), 1, 100);
        setBuilder.addPair("pair4", setBuilder.getHeader().getSequenceIndex("chr4"), 1, 100);
        setBuilder.addPair("pair5", setBuilder.getHeader().getSequenceIndex("chr5"), 1, 100);

        setBuilder.addUnmappedFragment("unmapped_frag");
        setBuilder.addUnmappedPair("unmapped_pair");

        setBuilder.addPair("split_pair1",
                setBuilder.getHeader().getSequenceIndex("chr1"),
                setBuilder.getHeader().getSequenceIndex("chr6"),
                100, 200,
                false, false, "36M",
                "36M", true, false, false,
                false, 30);

        setBuilder.addPair("split_pair2",
                setBuilder.getHeader().getSequenceIndex("chr1"),
                setBuilder.getHeader().getSequenceIndex("chr6"),
                100, 200,
                false, true, "36M",
                "*", true, false, false,
                false, 30);

        { // same dictionary, but different order
            final List<SAMSequenceRecord> sequences = new ArrayList<>(new SAMRecordSetBuilder().getHeader().getSequenceDictionary().getSequences());
            Collections.shuffle(sequences, new Random(42));

            tests.add(new Object[]{setBuilder, sequences, true, 0});
            tests.add(new Object[]{setBuilder, sequences, false, 0});
        }
        { // same dictionary, but different lengths
            final List<SAMSequenceRecord> sequences = new ArrayList<>(new SAMRecordSetBuilder().getHeader().getSequenceDictionary().getSequences());
            sequences.forEach(s -> s.setSequenceLength(s.getSequenceLength() + 1));

            tests.add(new Object[]{setBuilder, sequences, false, 1});
            tests.add(new Object[]{setBuilder, sequences, true, 1});
        }

        { // dropped contigs with split reads

            final List<SAMSequenceRecord> sequences = new ArrayList<>(new SAMRecordSetBuilder().getHeader().getSequenceDictionary().getSequences());
            sequences.remove(setBuilder.getHeader().getSequenceDictionary().getSequence("chr6"));
            sequences.remove(setBuilder.getHeader().getSequenceDictionary().getSequence("chr7"));
            sequences.remove(setBuilder.getHeader().getSequenceDictionary().getSequence("chr8"));
            sequences.remove(setBuilder.getHeader().getSequenceDictionary().getSequence("chr9"));

            Collections.shuffle(sequences, new Random(42));

            tests.add(new Object[]{setBuilder, sequences, true, 0});
            tests.add(new Object[]{setBuilder, sequences, false, 1});
        }

        { // dropped contigs with no reads

            final List<SAMSequenceRecord> sequences = new ArrayList<>(new SAMRecordSetBuilder().getHeader().getSequenceDictionary().getSequences());
            sequences.remove(setBuilder.getHeader().getSequenceDictionary().getSequence("chr7"));
            sequences.remove(setBuilder.getHeader().getSequenceDictionary().getSequence("chr8"));
            sequences.remove(setBuilder.getHeader().getSequenceDictionary().getSequence("chr9"));

            Collections.shuffle(sequences, new Random(42));

            tests.add(new Object[]{setBuilder, sequences, true, 0});
            tests.add(new Object[]{setBuilder, sequences, false, 1});
        }

        { // added contigs

            final List<SAMSequenceRecord> sequences = new ArrayList<>(new SAMRecordSetBuilder().getHeader().getSequenceDictionary().getSequences());
            sequences.add(new SAMSequenceRecord("test1", 100));
            sequences.add(new SAMSequenceRecord("test2", 100));
            sequences.add(new SAMSequenceRecord("test3", 100));
            sequences.add(new SAMSequenceRecord("test4", 100));
            Collections.shuffle(sequences, new Random(42));

            tests.add(new Object[]{setBuilder, sequences, true, 0});
            tests.add(new Object[]{setBuilder, sequences, false, 0});
        }

        { //dropped contigs with reads

            final List<SAMSequenceRecord> sequences = new ArrayList<>(new SAMRecordSetBuilder().getHeader().getSequenceDictionary().getSequences());
            sequences.remove(setBuilder.getHeader().getSequenceDictionary().getSequence("chr1"));
            sequences.remove(setBuilder.getHeader().getSequenceDictionary().getSequence("chr2"));

            tests.add(new Object[]{setBuilder, sequences, true, 0});
            tests.add(new Object[]{setBuilder, sequences, false, 1});

        }
        return tests.iterator();
    }

    @Test(dataProvider = "testDictionaryData")
    public void TestsWithIndex(final SAMRecordSetBuilder builder, final List<SAMSequenceRecord> sequences, final boolean allowIncomplete, final int expected) throws IOException {

        final File dictionary = File.createTempFile("reorder", ".dict");
        dictionary.deleteOnExit();
        writeDictionary(dictionary, sequences);

        final File bam = File.createTempFile("reorderIN" , ".bam");
        bam.deleteOnExit();
        final File bamIndex = new File(bam + ".bai");
        bamIndex.deleteOnExit();

        final File bamOut = File.createTempFile("reorderOUT", ".bam");
        bamOut.deleteOnExit();
        final File bamOutIndex = new File(bamOut + ".bai");
        bamOutIndex.deleteOnExit();

        tester(bam, bamOut, builder, dictionary, allowIncomplete, expected);
    }

    @Test(dataProvider = "testDictionaryData")
    public void TestsWithoutIndex(final SAMRecordSetBuilder builder, final List<SAMSequenceRecord> sequences, final boolean allowIncomplete, final int expected) throws IOException {
        final File sam = File.createTempFile("reorder", ".sam");
        sam.deleteOnExit();

        final File dictionary = File.createTempFile("reorder", ".dict");
        dictionary.deleteOnExit();
        writeDictionary(dictionary, sequences);

        final File bamOut = File.createTempFile("reorder", ".bam");
        bamOut.deleteOnExit();
        final File bamOutIndex = new File(bamOut + ".bai");
        bamOutIndex.deleteOnExit();

        tester(sam, bamOut, builder, dictionary, allowIncomplete, expected);
    }

    private void tester(File input, File output, final SAMRecordSetBuilder builder, final File dictionary, final boolean allowIncomplete, final int expected) {

        try (SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true).makeWriter(builder.getHeader(), false, input, null)) {
            for (final SAMRecord record : builder) {
                writer.addAlignment(record);
            }
        }

        final ValidateSamTester inputValidator = new ValidateSamTester();
        final List<SAMValidationError.Type> ignoreList = Collections.singletonList(SAMValidationError.Type.MISSING_TAG_NM);
        inputValidator.setIgnoreError(ignoreList);
        inputValidator.assertSamValid(input);

        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + output.getAbsolutePath(),
                "SEQUENCE_DICTIONARY=" + dictionary.getAbsolutePath(),
                "ALLOW_INCOMPLETE_DICT_CONCORDANCE=" + allowIncomplete
        };

        Assert.assertEquals(runPicardCommandLine(args), expected);

        if (expected == 0) {
            Assert.assertTrue(SAMSequenceDictionaryExtractor.extractDictionary(output.toPath()).isSameDictionary(
                    SAMSequenceDictionaryExtractor.extractDictionary(dictionary.toPath())));

            Assert.assertFalse(SAMSequenceDictionaryExtractor.extractDictionary(output.toPath()).isSameDictionary(
                    SAMSequenceDictionaryExtractor.extractDictionary(input.toPath())));

            final ValidateSamTester outputValidator = new ValidateSamTester();
            outputValidator.setIgnoreError(ignoreList);
            outputValidator.assertSamValid(output);
        }
    }

    private static void writeDictionary(final File dictionary, Collection<SAMSequenceRecord> sequences) throws IOException {
        try (BufferedWriter bufWriter = new BufferedWriter(new FileWriter(dictionary))) {
            SequenceDictionaryUtils.encodeDictionary(bufWriter, sequences.iterator());
        }
    }
}