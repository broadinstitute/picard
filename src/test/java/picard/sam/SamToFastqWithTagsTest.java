package picard.sam;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class SamToFastqWithTagsTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/bam2fastq/paired");

    public String getCommandLineProgramName() {
        return SamToFastqWithTags.class.getSimpleName();
    }

    private static final String[] STG_ARRAY = {"UR,CR", "CR,UR"};
    @DataProvider(name = "argumentCombinationsdata")
    public Object[][] argumentCombinationsdata() {
        return new Object[][]{
                // not supplying any arguments
                {null, null, null, 1},
                // only supply STGs and no QTG or SEP
                {STG_ARRAY, null, null, 0},
                // not supplying the correct number of QTGs or SEPs
                {STG_ARRAY, new String[]{"UY"}, null, 1},
                {STG_ARRAY, null, new String[]{"AAAAA"}, 1},
                // supplying correct number of QTGs or SEPs
                {new String[]{"UR,CR"}, new String[]{"CY,UY"}, null, 0},
                {STG_ARRAY, null, new String[]{"AAA", "AAAAA"}, 0},
                {STG_ARRAY, new String[]{"CY", "CY,UY"}, new String[]{"AAA", "AAAAA"}, 0},

        };
    }

    @Test(dataProvider = "argumentCombinationsdata")
    public void argumentCombinations(final String[] stgArguments,
                                     final String[] qtgArguments,
                                     final String[] sepArguments,
                                     final int returnCode) throws IOException {

        final File temp_input = File.createTempFile("input", ".fastq");
        temp_input.deleteOnExit();
        final File temp_input_2 = File.createTempFile("input_2", ".fastq");
        temp_input_2.deleteOnExit();
        final File temp_cr = File.createTempFile("UR_CR", ".fastq");
        temp_cr.deleteOnExit();
        final File temp_cr_ur = File.createTempFile("CR_UR", ".fastq");
        temp_cr_ur.deleteOnExit();

        final ArrayList<String> args = new ArrayList<>();
        args.add("INPUT=" + TEST_DATA_DIR + "/ok/sorted-pair.sam");
        args.add("FASTQ=" + temp_input.getAbsolutePath());
        args.add("SECOND_END_FASTQ=" + temp_input_2.getAbsolutePath());

        if (stgArguments != null) {
            for (String value : stgArguments) {
                args.add("STG=" + value);
            }
        }

        if (qtgArguments != null) {
            for (String value : qtgArguments) {
                args.add("QTG=" + value);
            }
        }

        if (sepArguments != null) {
            for (String value : sepArguments) {
                args.add("SEP=" + value);
            }
        }

        // make sure program invocation returns the expected return code
        Assert.assertEquals(runPicardCommandLine(args), returnCode);
    }

    @DataProvider(name = "queryNonExistantTag")
    public Object[][] queryNonExistantTag() {
        return new Object[][] {
                {"ok/sorted-single.sam"}
        };
    }

    @Test (dataProvider = "queryNonExistantTag", expectedExceptions =  PicardException.class)
    public void testQueryNonExistantTag(final String samFilename) throws IOException {
        final File temp_input = File.createTempFile("input", ".fastq");
        temp_input.deleteOnExit();

        runPicardCommandLine(new String[]{
                "INPUT=" + TEST_DATA_DIR + "/" + samFilename,
                "FASTQ=" + temp_input.getAbsolutePath(),
                "STG=BR"
        });
    }

    @DataProvider(name = "okGroupedFiles")
    public Object[][] okGroupedFiles() {
        return new Object[][] {
                {"ok/grouped-last-pair-mates-flipped.sam", new String[]{"rg1","rg2"}},
        };
    }

    // make sure the interleaved fastq and sequence tag fastq contain the same mate pairs
    @Test(dataProvider = "okGroupedFiles")
    public void testOkGroupedFiles(final String samFilename, final String [] groupFiles) throws IOException {
        final File samFile = new File(TEST_DATA_DIR,samFilename);
        final Map<String, Set<String>> outputSets = new HashMap<>(groupFiles.length);

        final String tmpDir = IOUtil.getDefaultTmpDir().getAbsolutePath() + "/";
        final String [] args = {
                "INPUT=" + samFile.getAbsolutePath(),
                "OUTPUT_PER_RG=true",
                "OUTPUT_DIR=" + tmpDir,
                "INTERLEAVE=true",
                "STG=CR,UR",
                "QTG=CY,UY"
        };
        runPicardCommandLine(args);

        Set<String> outputHeaderSet1;
        Set<String> outputHeaderSet2;
        for (final String groupPUName : groupFiles) {
            String keyName1 = groupPUName + "_1";
            String keyName2 = groupPUName + "_CR_UR";
            String fname1 = tmpDir + "/" + keyName1 + ".fastq";
            String fname2 = tmpDir + "/" + keyName2 + ".fastq";
            File f1 = new File(fname1);
            File f2 = new File(fname2);
            f1.deleteOnExit();
            f2.deleteOnExit();
            IOUtil.assertFileIsReadable(f1);
            IOUtil.assertFileIsReadable(f2);

            // Check that paired fastq files are same size and store them for later comparison
            outputHeaderSet1 = SamToFastqTest.createFastqReadHeaderSet(f1);
            outputHeaderSet2 = SamToFastqTest.createFastqReadHeaderSet(f2);
            outputSets.put(keyName1 , outputHeaderSet1);
            outputSets.put(keyName2, outputHeaderSet2);
            Assert.assertEquals(outputHeaderSet1.size(), outputHeaderSet2.size());
        }

        // Create map of read groups and mate pairs from SAM records
        final Map<String, Map<String,SamToFastqTest.MatePair>> map = SamToFastqTest.createPUPairsMap(samFile);

        for(final Map.Entry<String, Map<String, SamToFastqTest.MatePair>> groupEntry : map.entrySet()) {
            // Ensure that for each group, each mate of each pair in the SAM file is in the correct fastq pair file
            for (final Map.Entry<String,SamToFastqTest.MatePair> entry : groupEntry.getValue().entrySet() ) {
                final SamToFastqTest.MatePair mpair = entry.getValue();
                outputHeaderSet1 = outputSets.get(groupEntry.getKey() + "_1");
                outputHeaderSet2 = outputSets.get(groupEntry.getKey() + "_CR_UR");

                Assert.assertNotNull(mpair.mate1); // ensure we have two mates
                Assert.assertNotNull(mpair.mate2);
                Assert.assertEquals(mpair.mate1.getReadName(),mpair.mate2.getReadName());
                final String readName = mpair.mate1.getReadName() ;
                Assert.assertTrue(outputHeaderSet1.contains(readName + "/1")); // ensure mate is in correct file
                Assert.assertTrue(outputHeaderSet2.contains(readName + "/2"));
            }
        }
    }

    @DataProvider(name = "testTagGroupFastqData")
    public Object[][] testTagGroupFastqData() {
        return new Object[][] {
                {new String[]{"CR"}, new String[]{"CY"}, null, "AAAAA", "11111", "CCCCC", "22222", null, null, null, null, "One Sequence Tag Group"},
                {new String[]{"CR,UR"}, new String[]{"CY,UY"}, null, "AAAAATTT", "11111222", "CCCCCGGG", "22222111", null, null, null, null, "One Sequence Tag Group with Two Tags"},
                {new String[]{"CR,UR", "UR,CR"}, new String[]{"CY,UY", "UY,CY"}, null, "AAAAATTT", "11111222", "CCCCCGGG", "22222111", "TTTAAAAA", "22211111", "GGGCCCCC", "11122222", "Two Sequence Tag Groups with Two Tags"},
                {new String[]{"CR,UR", "UR,CR"}, new String[]{"CY,UY", "UY,CY"}, new String[]{"ATGC", "CGTA"} , "AAAAAATGCTTT", "11111~~~~222", "CCCCCATGCGGG", "22222~~~~111", "TTTCGTAAAAAA", "222~~~~11111", "GGGCGTACCCCC", "111~~~~22222", "Two Sequence Tag Groups with Two Tags with Separators"},

        };
    }

    @Test(dataProvider = "testTagGroupFastqData")
    public void testTagGroupFastq(final String[] sequenceTags, final String[] qualityTags, final String[] sepArgs, final String bases_1_1, final String quals_1_1,
                             final String bases_1_2, final String quals_1_2, final String bases_2_1, final String quals_2_1, final String bases_2_2,
                             final String quals_2_2, final String testName) throws IOException {
        final File samFile = new File(TEST_DATA_DIR, "ok/sorted-single.sam") ;
        final File fastq = File.createTempFile("testtag", "fastq");
        final File f1 = new File(fastq.getParent(), sequenceTags[0].replace(",", "_") + ".fastq");
        File f2 = null;
        if (quals_2_1 != null) {
            f2 = new File(fastq.getParent(), sequenceTags[1].replace(",", "_") + ".fastq");
            f2.deleteOnExit();
        }
        f1.deleteOnExit();


        final ArrayList<String> args = new ArrayList<>();
        args.add("INPUT=" + samFile.getAbsolutePath());
        args.add("FASTQ=" + fastq.getAbsolutePath());

        for (String value : sequenceTags) {
            args.add("STG=" + value);
        }
        if (qualityTags != null) {
            for (String value : qualityTags) {
                args.add("QTG=" + value);
            }
        }
        if (sepArgs != null) {
            for (String value : sepArgs) {
                args.add("SEP=" + value);
            }
        }

        Assert.assertEquals(runPicardCommandLine(args), 0);

        Iterator<FastqRecord> it = new FastqReader(f1).iterator();
        FastqRecord first = it.next();
        Assert.assertEquals(first.getReadString(), bases_1_1, testName);
        Assert.assertEquals(first.getBaseQualityString(), quals_1_1, testName);
        FastqRecord second = it.next();
        Assert.assertEquals(second.getReadString(), bases_1_2, testName);
        Assert.assertEquals(second.getBaseQualityString(), quals_1_2, testName);
        if (bases_2_1 != null) {
            it = new FastqReader(f2).iterator();
            first = it.next();
            Assert.assertEquals(first.getReadString(), bases_2_1, testName);
            Assert.assertEquals(first.getBaseQualityString(), quals_2_1, testName);
            second = it.next();
            Assert.assertEquals(second.getReadString(), bases_2_2, testName);
            Assert.assertEquals(second.getBaseQualityString(), quals_2_2, testName);
        }
    }
}
