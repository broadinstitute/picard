package picard.sam;

import htsjdk.samtools.SamStreams;
import htsjdk.samtools.cram.CRAMException;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgram;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

public class CramCompatibilityTest {

    public static final String CRAM_FILE = "testdata/picard/sam/test_cram_file_coordinate_sorted.cram";
    public static final String CRAM_FILE_2 = "testdata/picard/sam/test_cram_file_header_only.cram";
    public static final String CRAM_FILE_ONE_PAIR_MC = "testdata/picard/sam/MarkDuplicates/one_pair_mc.cram";

    public static final String CRAM_FILE_QUERY_SORTED_UNMAPPED = "testdata/picard/sam/unmapped_queryname_sorted.cram";
    public static final String CRAM_FILE_QUERY_SORTED = "testdata/picard/sam/test_cram_file_query_sorted.cram";

    public static final String REFERENCE_FILE = "testdata/picard/sam/test_cram_file.ref.fa";
    public static final String FASTQ_FILE = "testdata/picard/sam/fastq2bam/fastq-sanger/5k-v1-Rhodobacter_LW1.sam.fastq.gz";

    public static final String CRAM_UNMAPPED = "testdata/picard/sam/SamFormatConverterTest/unmapped.cram";
    public static final String CRAM_UNMAPPED_WITH_OQ_TAG = "testdata/picard/sam/unmapped_with_oq_tag.cram";

    public static final String CRAM_UNMAPPED_PART_1 = "testdata/picard/sam/unmapped_part_1.cram";
    public static final String CRAM_UNMAPPED_PART_2 = "testdata/picard/sam/unmapped_part_2.cram";

    public static final String CRAM_SPLIT_UNMAPPED = "testdata/picard/sam/split_test_unmapped.cram";

    public static final String MBA_ALIGNED_CRAM = "testdata/picard/sam/MergeBamAlignment/cliptest.aligned.cram";
    public static final String MBA_UNMAPPED_CRAM = "testdata/picard/sam/MergeBamAlignment/cliptest.unmapped.cram";
    public static final String MBA_REFERENCE = "testdata/picard/sam/MergeBamAlignment/cliptest.fasta";

    private static final File outputDir = IOUtil.createTempDir("testdata/picard/sam/CramCompatibilityTest", ".tmp");

    @AfterTest
    public void tearDown() {
        IOUtil.recursiveDelete(outputDir.toPath());
    }

    @DataProvider(name = "programArgsForCRAMWithReference")
    public Object[][] getArgsForCRAMWithReference() {
        return new Object[][] {
                {"picard.sam.AddOrReplaceReadGroups",
                        "RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20",
                        CRAM_FILE,
                        REFERENCE_FILE
                },
                {"picard.sam.CleanSam", null, CRAM_FILE, REFERENCE_FILE},
                {"picard.sam.DownsampleSam", null, CRAM_FILE, REFERENCE_FILE},
                {"picard.sam.FixMateInformation", null, CRAM_FILE, REFERENCE_FILE},
                {"picard.sam.markduplicates.MarkDuplicates",
                        "M=" + createTempFile("MarkDuplicates", ".dir").getAbsolutePath(),
                        CRAM_FILE,
                        REFERENCE_FILE
                },
                {"picard.sam.MergeSamFiles", null, CRAM_FILE, REFERENCE_FILE},
                {"picard.sam.PositionBasedDownsampleSam", "FRACTION=0.5", CRAM_FILE, REFERENCE_FILE},
                {"picard.sam.SortSam", "SORT_ORDER=queryname", CRAM_FILE, REFERENCE_FILE},
                {"picard.sam.ReplaceSamHeader", "HEADER=" + CRAM_FILE_2, CRAM_FILE, REFERENCE_FILE},
                {"picard.sam.RevertOriginalBaseQualitiesAndAddMateCigar", "CREATE_INDEX=false", CRAM_FILE_QUERY_SORTED, REFERENCE_FILE},
                {"picard.sam.GatherBamFiles",
                        "I=" + new File(CRAM_UNMAPPED).getAbsolutePath(),
                        CRAM_FILE_QUERY_SORTED,
                        REFERENCE_FILE
                },
                {"picard.sam.markduplicates.MarkDuplicatesWithMateCigar",
                        "M=" + createTempFile("MarkDuplicatesWithMateCigar", ".txt").getAbsolutePath(),
                        CRAM_FILE,
                        REFERENCE_FILE
                },
                {"picard.sam.markduplicates.SimpleMarkDuplicatesWithMateCigar",
                        "M=" + createTempFile("SimpleMarkDuplicatesWithMateCigar", ".txt").getAbsolutePath(),
                        CRAM_FILE_ONE_PAIR_MC,
                        REFERENCE_FILE
                },
                {"picard.sam.ReorderSam",
                        "SEQUENCE_DICTIONARY=" + REFERENCE_FILE,
                        CRAM_FILE,
                        REFERENCE_FILE
                },
                {"picard.sam.SetNmMdAndUqTags", null, CRAM_FILE, REFERENCE_FILE},
                {"picard.sam.MergeBamAlignment",
                        "UNMAPPED=" + new File(MBA_UNMAPPED_CRAM).getAbsolutePath() +
                        " ALIGNED=" + new File(MBA_ALIGNED_CRAM).getAbsolutePath(),
                        null,
                        MBA_REFERENCE
                },
                {"picard.illumina.MarkIlluminaAdapters",
                        "METRICS=" + createTempFile("picard.illumina.MarkIlluminaAdapters", ".txt").getAbsolutePath(),
                        CRAM_FILE_QUERY_SORTED,
                        REFERENCE_FILE
                },
                {"picard.sam.SplitSamByLibrary", null, CRAM_FILE, REFERENCE_FILE}
        };
    }

    @Test(dataProvider = "programArgsForCRAMWithReference")
    public void testShouldWriteCRAMWhenCRAMWithReference(String program,
                                                         String parameters,
                                                         String cramFile,
                                                         String reference) throws IOException, IllegalAccessException, InstantiationException, ClassNotFoundException {
        if (!program.equals("picard.sam.SplitSamByLibrary")) {
            final File outputFile = createTempCram(program);
            launchProgram(program, cramFile, outputFile.getAbsolutePath(), parameters, reference);
            assertCRAM(outputFile);
        } else {
            final File tmpDir = IOUtil.createTempDir(outputDir.getAbsolutePath(), program);
            launchProgram(program, cramFile, tmpDir.getAbsolutePath(), parameters, reference);
            assertCRAMs(tmpDir);
        }
    }

    @DataProvider(name  = "programArgsForCRAMWithoutReferenceToFail")
    public Object[][] getArgsForCRAMWithoutReferenceToFail() {
        return new Object[][] {
                {"picard.sam.AddOrReplaceReadGroups",
                        "RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20",
                        CRAM_FILE
                },
                {"picard.sam.CleanSam", null, CRAM_FILE},
                {"picard.sam.DownsampleSam", null, CRAM_FILE},
                {"picard.sam.FixMateInformation", null, CRAM_FILE},
                {"picard.sam.markduplicates.MarkDuplicates",
                        "M=" + createTempFile("MarkDuplicates", ".dir").getAbsolutePath(),
                        CRAM_FILE
                },
                {"picard.sam.MergeSamFiles", null, CRAM_FILE},
                {"picard.sam.PositionBasedDownsampleSam", "FRACTION=0.5", CRAM_FILE},
                {"picard.sam.SortSam", "SORT_ORDER=queryname", CRAM_FILE},
                {"picard.sam.ReplaceSamHeader", "HEADER=" + CRAM_FILE_2, CRAM_FILE},
                {"picard.sam.RevertOriginalBaseQualitiesAndAddMateCigar", null, CRAM_FILE_QUERY_SORTED},
                {"picard.sam.GatherBamFiles",
                        "I=" + new File(CRAM_UNMAPPED).getAbsolutePath(),
                        CRAM_FILE_QUERY_SORTED
                },
                {"picard.sam.markduplicates.MarkDuplicatesWithMateCigar",
                        "M=" + createTempFile("MarkDuplicatesWithMateCigar", ".txt").getAbsolutePath(),
                        CRAM_FILE},
                {"picard.sam.markduplicates.SimpleMarkDuplicatesWithMateCigar",
                        "M=" + createTempFile("SimpleMarkDuplicatesWithMateCigar", ".txt").getAbsolutePath(),
                        CRAM_FILE_ONE_PAIR_MC},
                {"picard.illumina.MarkIlluminaAdapters",
                        "METRICS=" + createTempFile("picard.illumina.MarkIlluminaAdapters", ".txt").getAbsolutePath(),
                        CRAM_FILE_QUERY_SORTED,
                },
                {"picard.sam.SplitSamByLibrary", null, CRAM_FILE}
        };
    }

    @Test(dataProvider = "programArgsForCRAMWithoutReferenceToFail", expectedExceptions = {CRAMException.class, IllegalArgumentException.class})
    public void testShouldFailWhenCRAMWithoutReference(String program,
                                                       String parameters,
                                                       String cramFile) throws IOException, IllegalAccessException, InstantiationException, ClassNotFoundException {
        if (!program.equals("picard.sam.SplitSamByLibrary")) {
            final File outputFile = createTempCram(program);
            launchProgram(program, cramFile, outputFile.getAbsolutePath(), parameters, null);
            assertCRAM(outputFile);
        } else {
            final File tmpDir = IOUtil.createTempDir(outputDir.getAbsolutePath(), program);
            launchProgram(program, cramFile, tmpDir.getAbsolutePath(), parameters, null);
            assertCRAMs(tmpDir);
        }
    }

    // test with CRAMs that don't need reference (unmapped CRAMs for input or output)
    @DataProvider(name = "programArgsWithUnmappedCRAM")
    public Object[][] getArgsWithUnmappedCRAM() {
        return new Object[][] {
                {"picard.sam.AddOrReplaceReadGroups", "RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20", CRAM_UNMAPPED},
                {"picard.sam.CleanSam", null, CRAM_UNMAPPED},
                {"picard.sam.DownsampleSam", null, CRAM_UNMAPPED},
                {"picard.sam.FixMateInformation", null, CRAM_UNMAPPED},
                {"picard.sam.markduplicates.MarkDuplicates",
                        "M=" + createTempFile("MarkDuplicates", ".dir").getAbsolutePath(),
                        CRAM_UNMAPPED
                },
                {"picard.sam.MergeSamFiles", null, CRAM_UNMAPPED},
                {"picard.sam.PositionBasedDownsampleSam", "FRACTION=0.5", CRAM_UNMAPPED},
                {"picard.sam.SortSam", "SORT_ORDER=queryname", CRAM_UNMAPPED},
                {"picard.sam.ReplaceSamHeader", "HEADER=" + MBA_UNMAPPED_CRAM, CRAM_UNMAPPED},
                {"picard.sam.RevertOriginalBaseQualitiesAndAddMateCigar", "CREATE_INDEX=false", CRAM_UNMAPPED_WITH_OQ_TAG},
                {"picard.sam.GatherBamFiles",
                        "I=" + new File(CRAM_UNMAPPED_PART_2).getAbsolutePath(),
                        CRAM_UNMAPPED_PART_1
                },
                {"picard.sam.FastqToSam", "F1=" + FASTQ_FILE + " SAMPLE_NAME=s1", null},
                {"picard.illumina.IlluminaBasecallsToSam",
                        "BASECALLS_DIR=" + new File("testdata/picard/illumina/25T8B25T/Data/Intensities/BaseCalls") +
                        " LANE=1 READ_STRUCTURE=25S8S25T RUN_BARCODE=HiMom SAMPLE_ALIAS=HiDad LIBRARY_NAME=HelloWorld SEQUENCING_CENTER=BI" ,
                        null
                },
                {"picard.illumina.MarkIlluminaAdapters",
                        "METRICS=" + createTempFile("picard.illumina.MarkIlluminaAdapters", ".txt").getAbsolutePath(),
                        CRAM_FILE_QUERY_SORTED_UNMAPPED
                },
                {"picard.sam.SplitSamByLibrary", null, CRAM_SPLIT_UNMAPPED}
        };
    }

    @Test(dataProvider = "programArgsWithUnmappedCRAM")
    public void testShouldWriteCRAMWhenUnmappedCRAMWithoutReference(String program,
                                                                    String parameters,
                                                                    String cramFile) throws IOException, IllegalAccessException, InstantiationException, ClassNotFoundException {
        if (!program.equals("picard.sam.SplitSamByLibrary")) {
            final File outputFile = createTempCram(program);
            launchProgram(program, cramFile, outputFile.getAbsolutePath(), parameters, null);
            assertCRAM(outputFile);
        } else {
            final File tmpDir = IOUtil.createTempDir(outputDir.getAbsolutePath(), program);
            launchProgram(program, cramFile, tmpDir.getAbsolutePath(), parameters, null);
            assertCRAMs(tmpDir);
        }
    }

    private File createTempCram(String name) throws IOException {
        return createTempFile(name, ".cram");
    }

    private static File createTempFile(String name, String extension) {
        File file = null;
        try {
            file = File.createTempFile(name, extension, outputDir);
            file.deleteOnExit();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return file;
    }

    private void launchProgram(String programClassname,
                               String input,
                               String output,
                               String exParams,
                               String reference) throws ClassNotFoundException, IllegalAccessException, InstantiationException {
        final Collection<String> args = new ArrayList<>();

        if (input != null) {
            args.add("INPUT=" + new File(input).getAbsolutePath());
        }
        args.add("OUTPUT=" + output);

        if (exParams != null) {
            args.addAll(Arrays.asList(exParams.split(" ")));
        }

        if (reference != null) {
            args.add("REFERENCE_SEQUENCE=" + new File(reference).getAbsolutePath());
        }

        final CommandLineProgram program = (CommandLineProgram) Class.forName(programClassname).newInstance();
        program.instanceMain(args.toArray(new String[0]));
    }

    static void assertCRAM(final File outputFile) {
        Assert.assertTrue(outputFile.exists(), "Expected output file " + outputFile.getAbsolutePath() + " doesn't exist.");
        try (InputStream in = new FileInputStream(outputFile)) {
            Assert.assertTrue(SamStreams.isCRAMFile(new BufferedInputStream(in)), "File " + outputFile.getAbsolutePath() + " is not a CRAM.");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void assertCRAMs(File dir) {
        Arrays.stream(dir.listFiles()).filter(file -> file.getName().endsWith("cram")).forEach(CramCompatibilityTest::assertCRAM);
    }
}
