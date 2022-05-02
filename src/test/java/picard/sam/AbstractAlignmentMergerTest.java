package picard.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.cmdline.argumentcollections.RequiredReferenceArgumentCollection;
import picard.nio.PicardHtsPath;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Tests related to code in AbstractAlignmentMerger
 */
public class AbstractAlignmentMergerTest extends CommandLineProgramTest {

    // the actual data provider is overlapReadDataWithSwaps
    public Object[][] overlapReadData() {
        // The spaces here are deliberate to illustrate the region the two default reads match
        final String default120LongR1Bases = "ATCACACCAGTGTCTGCGTTCACAGCAGGCATCATCAGTAGCCTCCAGAGGCCTCAGGTCCAGTCTCTAAAAATATCTCAGGAGGCTGCAGTGGCTGACCAGATTCTCCTGTCAGTTTGC";
        final String default120LongR2Bases = "CGTTGGCAATGCCGGGCACAATCACACCAGTGTCTGCGTTCACAGCAGGCATCATCAGTAGCCTCCAGAGGCCTCAGGTCCAGTCTCTAAAAATATCTCAGGAGGCTGCAGTGGCTGACC";

        final String default110LongR1Bases = "ATCACACCAGTGTCTGCGTTCACAGCAGGCATCATCAGTAGCCTCCAGAGGCCTCAGGTCCAGTCTCTAAAAATATCTCAGGAGGCTGCAGTGGCTGACCAGATTCTCCT";
        final String default110LongR2Bases = "GCCGGGCACAATCACACCAGTGTCTGCGTTCACAGCAGGCATCATCAGTAGCCTCCAGAGGCCTCAGGTCCAGTCTCTAAAAATATCTCAGGAGGCTGCAGTGGCTGACC";

        final String sharedBases = "ATCACACCAGTGTCTGCGTTCACAGCAGGCATCATCAGTAGCCTCCAGAGGCCTCAGGTCCAGTCTCTAAAAATATCTCAGGAGGCTGCAGTGGCTGACC";

        final String default120LongR1ClippedBases = "AGATTCTCCTGTCAGTTTGC";
        final String default120LongR2ClippedBases = "TGTGCCCGGCATTGCCAACG";

        final String default110LongR1ClippedBases = "AGATTCTCCT";
        final String default110LongR2ClippedBases = "TGTGCCCGGC";

        final String default27LongR1Bases = "AGATTCTCCTTGTGCCCGGCAGATTCT";
        final String default27LongR2Bases = "TGTGCCCGGCAGATTCTCCTCTTGTGC";

        final String default27LongR1Qualities = "ABCDEFGHIJKLMNOPQRSTUVWXYZ.";
        final String default27LongR2Qualities = "abcdefghijklmnopqrstuvwxyz,";


        final String default120LongR1BaseQualities = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.FFF.FFF.FFF";
        final String default120LongR2BaseQualities = "FFFFFF.FFFFF.FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
        final String default110LongR1BaseQualities = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.FFF.FFF";
        final String default110LongR2BaseQualities = "FFFFFF.FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";

        final String sharedQualities = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";

        final String r1ClippedQualities10 = default120LongR1BaseQualities.substring(default120LongR1BaseQualities.length() - 10);
        final String r2ClippedQualities10 = StringUtil.reverseString(default120LongR2BaseQualities.substring(0, 10));
        final String r1ClippedQualities20 = default120LongR1BaseQualities.substring(default120LongR1BaseQualities.length() - 20);
        final String r2ClippedQualities20 = StringUtil.reverseString(default120LongR2BaseQualities.substring(0, 20));


        return new Object[][]{
                {110, false, 100, 200, "110M", "110M", false, true, 100, 200, "110M", "110M",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null}, // Non overhanging reads
                {110, false, 100, 200, "110M", "110M", false, true, 100, 200, "110M", "110M",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null},

                {110, false, 100, 300, "110M", "110M", false, true, 100, 300, "110M", "110M",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null}, // Non overlapping reads
                {110, false, 100, 300, "110M", "110M", false, true, 100, 300, "110M", "110M",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null},

                {110, false, 100, 210, "110M", "110M", true, false, 100, 210, "110M", "110M",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null}, // Non overlapping reads, outies, abutting
                {110, false, 100, 210, "110M", "110M", true, false, 100, 210, "110M", "110M",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null},

                {110, false, 100, 210, "110M", "110M", false, true, 100, 210, "110M", "110M",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null}, // Non overlapping reads, innies, abutting
                {110, false, 100, 210, "110M", "110M", false, true, 100, 210, "110M", "110M",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null},

                {110, false, 100, 209, "110M", "110M", true, false, 209, 209, "109S1M", "1M109S",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null}, // overlap by one base only
                {110, true, 100, 209, "110M", "110M", true, false, 209, 209, "109H1M", "1M109H",
                        default110LongR1Bases, default110LongR2Bases, "T", "G", SequenceUtil.reverseComplement(default110LongR1Bases.substring(0, default110LongR1Bases.length() - 1)), default110LongR2Bases.substring(1),
                        default110LongR1BaseQualities, default110LongR2BaseQualities, "F", "F", StringUtil.reverseString(default110LongR1BaseQualities.substring(0, default110LongR1BaseQualities.length() - 1)), default110LongR2BaseQualities.substring(1)},

                {110, true, 100, 200, "110M", "110M", false, true, 100, 200, "110M", "110M",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null},

                {110, false, 100, 200, "110M", "110M", false, false, 100, 200, "110M", "110M",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null}, // F1F2

                {110, true, 100, 200, "110M", "110M", false, false, 100, 200, "110M", "110M",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null},

                {110, false, 100, 200, "110M", "110M", true, true, 100, 200, "110M", "110M",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null}, // R1R2

                {110, true, 100, 200, "110M", "110M", true, true, 100, 200, "110M", "110M",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null},

                {110, false, 100, 300, "110M", "110M", true, false, 100, 300, "110M", "110M",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null}, // Non overlapping "outies"

                {110, true, 100, 300, "110M", "110M", true, false, 100, 300, "110M", "110M",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null},

                {110, false, 100, 90, "110M", "110M", false, true, 100, 100, "100M10S", "10S100M",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null}, // Basic overlapped read

                {110, true, 100, 90, "110M", "110M", false, true, 100, 100, "100M10H", "10H100M",
                        default110LongR1Bases, default110LongR2Bases, sharedBases, sharedBases, default110LongR1ClippedBases, default110LongR2ClippedBases,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, sharedQualities, sharedQualities, r1ClippedQualities10, r2ClippedQualities10},

                {120, false, 100, 95, "110M10S5H", "5H15S105M", false, true, 100, 100, "100M20S5H", "5H20S100M",
                        default120LongR1Bases, default120LongR2Bases, default120LongR1Bases, default120LongR2Bases, null, null,
                        default120LongR1BaseQualities, default120LongR2BaseQualities, default120LongR1BaseQualities, default120LongR2BaseQualities, null, null}, // Already hard and soft clipped

                {120, true, 100, 95, "110M10S5H", "5H15S105M", false, true, 100, 100, "100M25H", "25H100M",
                        default120LongR1Bases, default120LongR2Bases, sharedBases, sharedBases, default120LongR1ClippedBases, default120LongR2ClippedBases,
                        default120LongR1BaseQualities, default120LongR2BaseQualities, sharedQualities, sharedQualities, r1ClippedQualities20, r2ClippedQualities20},

                {120, false, 100, 95, "110M10S", "15S105M", false, true, 100, 100, "100M20S", "20S100M",
                        default120LongR1Bases, default120LongR2Bases, default120LongR1Bases, default120LongR2Bases, null, null,
                        default120LongR1BaseQualities, default120LongR2BaseQualities, default120LongR1BaseQualities, default120LongR2BaseQualities, null, null}, // Already soft clipped

                {120, true, 100, 95, "110M10S", "15S105M", false, true, 100, 100, "100M20H", "20H100M",
                        default120LongR1Bases, default120LongR2Bases, sharedBases, sharedBases, default120LongR1ClippedBases, default120LongR2ClippedBases,
                        default120LongR1BaseQualities, default120LongR2BaseQualities, sharedQualities, sharedQualities, r1ClippedQualities20, r2ClippedQualities20},

                {120, true, 100, 95, "95M25S", "15S105M", false, true, 100, 100, "95M5S20H", "20H100M",
                        default120LongR1Bases, default120LongR2Bases, sharedBases, sharedBases, default120LongR1ClippedBases, default120LongR2ClippedBases,
                        default120LongR1BaseQualities, default120LongR2BaseQualities, sharedQualities, sharedQualities, r1ClippedQualities20, r2ClippedQualities20},

                // 5' end soft clip tests
                /*
                                       SSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM->
                                 <-MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSSSSS

                                 soft clipping becomes:
                                       SSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSSSSSSSSS->
                                 <-SSSSSSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSSSSS

                                 hard clipping becomes:
                                       SSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSSSSSHHHH->
                                 <-HHHHSSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSSSSS
                */

                {110, false, 105, 90, "5S105M", "103M7S", false, true, 105, 105, "5S88M17S", "15S88M7S",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null}, // Already soft clipped at 5' end

                {110, true, 105, 90, "5S105M", "103M7S", false, true, 105, 105, "5S88M7S10H", "10H5S88M7S",
                        default110LongR1Bases, default110LongR2Bases, sharedBases, sharedBases, default110LongR1ClippedBases, default110LongR2ClippedBases,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, sharedQualities, sharedQualities, r1ClippedQualities10, r2ClippedQualities10},

                // 3' end soft clip tests
                /*
                                       SSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM>
                                 <SSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSS

                                 soft clipping becomes:
                                       SSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSSSSSS>
                                 <SSSSSSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSS

                                 hard clipping becomes:
                                       SSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSHHHHH>
                                 <HHHHHSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSS
                */

                {110, false, 105, 100, "10S100M", "10S95M5S", false, true, 105, 105, "10S90M10S", "15S90M5S",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null}, // Already soft clipped at 5' end

                {110, true, 105, 100, "10S100M", "10S95M5S", false, true, 105, 105, "10S90M5S5H", "5H10S90M5S",
                        default110LongR1Bases, default110LongR2Bases,
                        default110LongR1Bases.substring(0, 105), default110LongR2Bases.substring(5),
                        default110LongR1Bases.substring(105), SequenceUtil.reverseComplement(default110LongR2Bases.substring(0, 5)),
                        default110LongR1BaseQualities, default110LongR2BaseQualities,
                        default110LongR1BaseQualities.substring(0, 105), default110LongR2BaseQualities.substring(5),
                        default110LongR1BaseQualities.substring(105), StringUtil.reverseString(default110LongR2BaseQualities.substring(0, 5))}, // Already soft clipped at 5' end

                /*
                                SSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM->
                                     <-MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSSSSS

                                soft clipping becomes:
                                SSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSS->
                                     <-SSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSSSSS

                                hard-clipping results in the same as with soft-clipping in this case.

                 */

                {110, false, 105, 100, "10S100M", "103M7S", false, true, 105, 105, "10S98M2S", "5S98M7S",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null},

                {110, false, 105, 100, "10S100M", "103M7S", false, true, 105, 105, "10S98M2S", "5S98M7S",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null},

                /*
                                SSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSS->
                                     <-MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSSSSSSSSSSSSSSSSSS

                                soft clipping becomes:
                                SSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSSSSSSSSSS->
                                     <-SSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSSSSSSSSSSSSSSSSSS

                                hard-clipping results in the same as with soft-clipping in this case.

                 */
                {110, false, 105, 100, "10S97M3S", "99M11S", false, true, 105, 105, "10S94M6S", "5S94M11S",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null},

                {110, false, 105, 100, "10S97M3S", "99M11S", false, true, 105, 105, "10S94M6S", "5S94M11S",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null},

                /*
                                          SSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSS->
                                     <-SSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM

                                soft clipping becomes:
                                          SSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSS->
                                     <-SSSSSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM


                                hard clipping becomes:
                                          SSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSHHH->
                                     <-HHHSSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM

                 */

                {110, false, 105, 96, "12S80M18S", "13S97M", false, true, 105, 105, "12S80M18S", "22S88M",
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null},

                {110, true, 105, 96, "12S80M18S", "13S97M", false, true, 105, 105, "12S80M8S10H", "10H12S88M",
                        default110LongR1Bases, default110LongR2Bases, sharedBases, sharedBases, default110LongR1ClippedBases, default110LongR2ClippedBases,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, sharedQualities, sharedQualities, r1ClippedQualities10, r2ClippedQualities10},

                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //                  MMMMMMMMMMMMMMMMMMsssssssss>
                //     <ssssssssMMM-MMMMMMMMMMMMMMMM
                //          should become
                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //                  MMMMMMMMMMMMMMMMSSsssssssss>
                //      <ssssssssSSSMMMMMMMMMMMMMMMM


                {27, false, 18, 14, "18M9S", "8S3M1D16M", false, true, 18, 18, "16M11S", "11S16M",
                        default27LongR1Bases, default27LongR2Bases, default27LongR1Bases, default27LongR2Bases, null, null,
                        default27LongR1Qualities, default27LongR2Qualities, default27LongR1Qualities, default27LongR2Qualities, null, null},

                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //                  MMMMMMMMMMMMMMMMMMsssssssss>
                //   <ssssssssMMM---MMMMMMMMMMMMMMMM
                //          should become
                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //                  MMMMMMMMMMMMMMMMSSsssssssss>
                //      <ssssssssSSSMMMMMMMMMMMMMMMM


                {27, false, 18, 12, "18M9S", "8S3M3D16M", false, true, 18, 18, "16M11S", "11S16M",
                        default27LongR1Bases, default27LongR2Bases, default27LongR1Bases, default27LongR2Bases, null, null,
                        default27LongR1Qualities, default27LongR2Qualities, default27LongR1Qualities, default27LongR2Qualities, null, null},

                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //     <ssssssssMMMMMMMMMMMMMMMMMMM
                //                  MMMMMMMMMMMMMMM-MMMsssssssss>
                //          should become
                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //     <ssssssssSSSSMMMMMMMMMMMMMMM
                //                  MMMMMMMMMMMMMMMSSSsssssssss>


                {27, false, 14, 18, "8S19M", "15M1D3M9S", true, false, 18, 18, "12S15M", "15M12S",
                        default27LongR1Bases, default27LongR2Bases, default27LongR1Bases, default27LongR2Bases, null, null,
                        default27LongR1Qualities, default27LongR2Qualities, default27LongR1Qualities, default27LongR2Qualities, null, null},


                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //     <ssssssssMMMMMMMMMMMMMMMMMMM
                //                  MMMMMMMMMMMMMMM---MMMsssssssss>
                //          should become
                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //     <ssssssssSSSSMMMMMMMMMMMMMMM
                //                  MMMMMMMMMMMMMMMSSSsssssssss>


                {27, false, 14, 18, "8S19M", "15M3D3M9S", true, false, 18, 18, "12S15M", "15M12S",
                        default27LongR1Bases, default27LongR2Bases, default27LongR1Bases, default27LongR2Bases, null, null,
                        default27LongR1Qualities, default27LongR2Qualities, default27LongR1Qualities, default27LongR2Qualities, null, null},

                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //     <ssssssssMMMMMMMMMMMMMMMMMMM
                //                  MMMMMMMMMMMMMM---MMMMsssssssss>
                //          should become
                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //     <ssssssssSSSSMMMMMMMMMMMMMMM
                //                  MMMMMMMMMMMMMMSSSSsssssssss>
                {27, false, 14, 18, "8S19M", "14M3D4M9S", true, false, 18, 18, "12S15M", "14M13S",
                        default27LongR1Bases, default27LongR2Bases, default27LongR1Bases, default27LongR2Bases, null, null,
                        default27LongR1Qualities, default27LongR2Qualities, default27LongR1Qualities, default27LongR2Qualities, null, null},

                // 123456789.123456789.123456789.12-3456789.123456789.123456789.123456789.123456789
                //     <ssssssssMMMMMMMMMMMMMMMMMMM
                //      MMMMMMMMMMMMMMMMMMMMMMMMMMM-> ## deletion right at the end
                //          should remain
                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //     <ssssssssMMMMMMMMMMMMMMMMMMM
                //      MMMMMMMMMMMMMMMMMMMMMMMMMMM->  ## deletion right at the end
                {27, false, 14, 6, "8S19M", "27M1D", true, false, 14, 6, "8S19M", "27M1D",
                        default27LongR1Bases, default27LongR2Bases, default27LongR1Bases, default27LongR2Bases, null, null,
                        default27LongR1Qualities, default27LongR2Qualities, default27LongR1Qualities, default27LongR2Qualities, null, null},


                //                  123456789.123456--789.123456789.123456789.123456789.123456789.123456789.123456789
                //                                     MMMMMMMMMMMMMMMMMMsssssssss>
                //                          <sssssMMiiMMMMMMMMMMMMMMMMMM
                //                                  ^^
                //   with an insertion of two bases here should become
                //                    123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //                                     MMMMMMMMMMMMMMMMMSsssssssss>
                //                          <sssssSSSSSMMMMMMMMMMMMMMMMM
                {27, false, 18, 15, "18M9S", "5S2M2I18M", false, true, 18, 18, "17M10S", "10S17M",
                        default27LongR1Bases, default27LongR2Bases, default27LongR1Bases, default27LongR2Bases, null, null,
                        default27LongR1Qualities, default27LongR2Qualities, default27LongR1Qualities, default27LongR2Qualities, null, null},


                //       123456789.123456789.123456789.1234--56789.123456789.123456789.123456789.123456789
                //            <sssssssMMMMMMMMMMMMMMMMMMMM
                //                        MMMMMMMMMMMMMMMMMiiMsssssss>
                //                                         ^^
                //         with an insertion of two bases here should become
                //       123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //            <sssssssSSSSMMMMMMMMMMMMMMMM
                //                        MMMMMMMMMMMMMMMMSSSSsssssss>
                {27, false, 14, 18, "7S20M", "17M2I1M7S", true, false, 18, 18, "11S16M", "16M11S",
                        default27LongR1Bases, default27LongR2Bases, default27LongR1Bases, default27LongR2Bases, null, null,
                        default27LongR1Qualities, default27LongR2Qualities, default27LongR1Qualities, default27LongR2Qualities, null, null},


                //       123456789.123456789.123456789.123--456789.123456789.123456789.123456789.123456789
                //            <sssssssMMMMMMMMMMMMMMMMMMMM
                //                        MMMMMMMMMMMMMMMMiiMMsssssss>
                //                                        ^^
                //        with an insertion of two bases here should become
                //       123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //            <sssssssSSSSMMMMMMMMMMMMMMMM
                //                        MMMMMMMMMMMMMMMMSSSSsssssss>
                {27, false, 14, 18, "7S20M", "16M2I2M7S", true, false, 18, 18, "11S16M", "16M11S",
                        default27LongR1Bases, default27LongR2Bases, default27LongR1Bases, default27LongR2Bases, null, null,
                        default27LongR1Qualities, default27LongR2Qualities, default27LongR1Qualities, default27LongR2Qualities, null, null},


                //       123456789.123456789.123456789.123-456789.123456789.123456789.123456789.123456789
                //            <sssssssMMMMMMMMMMMMMMMMMMMM
                //                        MMMMMMMMMMMMMMMMiMMMsssssss>
                //                                        ^
                //        with an insertion of one base here should become
                //       123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //            <ssssssSSSSSMMMMMMMMMMMMMMMM
                //                        MMMMMMMMMMMMMMMMSSSSsssssss>
                {27, false, 14, 18, "7S20M", "16M1I3M7S", true, false, 18, 18, "11S16M", "16M11S",
                        default27LongR1Bases, default27LongR2Bases, default27LongR1Bases, default27LongR2Bases, null, null,
                        default27LongR1Qualities, default27LongR2Qualities, default27LongR1Qualities, default27LongR2Qualities, null, null},

                //       123456789.123456789.123456789.1234--56789.123456789.123456789.123456789.123456789
                //            <sssssssMMMMMMMMMMMMMMMMMMMM
                //                        MMMMMMMMMMMMMMMMMiiMMssssss>
                //                                         ^^
                //         with an insertion of two bases here should become
                //       123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //            <sssssssSSSSMMMMMMMMMMMMMMMM
                //                        MMMMMMMMMMMMMMMMSSSSSssssss>
                {27, false, 14, 18, "7S20M", "17M2I2M6S", true, false, 18, 18, "11S16M", "16M11S",
                        default27LongR1Bases, default27LongR2Bases, default27LongR1Bases, default27LongR2Bases, null, null,
                        default27LongR1Qualities, default27LongR2Qualities, default27LongR1Qualities, default27LongR2Qualities, null, null},


                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //                  MMMMMMMMMMMMMMMMMMsssssssss>
                //     <ssssssssMMM-MMMMMMMMMMMMMMMM
                //          should become
                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //                  MMMMMMMMMMMMMMMMHHHHHHHHHHH>
                //      <HHHHHHHHHHHMMMMMMMMMMMMMMMM
                {27, true, 18, 14, "18M9S", "8S3M1D16M", false, true, 18, 18, "16M11H", "11H16M",
                        default27LongR1Bases, default27LongR2Bases,
                        default27LongR1Bases.substring(0, 16), default27LongR2Bases.substring(11, 27),
                        default27LongR1Bases.substring(16, 27), SequenceUtil.reverseComplement(default27LongR2Bases.substring(0, 11)),
                        default27LongR1Qualities, default27LongR2Qualities,
                        default27LongR1Qualities.substring(0, 16), default27LongR2Qualities.substring(11, 27),
                        default27LongR1Qualities.substring(16, 27), StringUtil.reverseString(default27LongR2Qualities.substring(0, 11)),
                },


                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //                  MMMMMMMMMMMMMMMMMMsssssssss>
                //   <ssssssssMMM---MMMMMMMMMMMMMMMM
                //          should become
                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //                  MMMMMMMMMMMMMMMMHHHHHHHHHHH>
                //      <HHHHHHHHHHHMMMMMMMMMMMMMMMM
                {27, true, 18, 12, "18M9S", "8S3M3D16M", false, true, 18, 18, "16M11H", "11H16M",
                        default27LongR1Bases, default27LongR2Bases,
                        default27LongR1Bases.substring(0, 16), default27LongR2Bases.substring(11, 27),
                        default27LongR1Bases.substring(16, 27), SequenceUtil.reverseComplement(default27LongR2Bases.substring(0, 11)),
                        default27LongR1Qualities, default27LongR2Qualities,
                        default27LongR1Qualities.substring(0, 16), default27LongR2Qualities.substring(11, 27),
                        default27LongR1Qualities.substring(16, 27), StringUtil.reverseString(default27LongR2Qualities.substring(0, 11))
                },


                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //     <ssssssssSSSSMMMMMMMMMMMMMMM
                //                  MMMMMMMMMMMMMMMSSSsssssssss>
                //          should become
                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //     <HHHHHHHHHHHHMMMMMMMMMMMMMMM
                //                  MMMMMMMMMMMMMMMHHHHHHHHHHHH>
                {27, true, 14, 18, "8S19M", "15M1D3M9S", true, false, 18, 18, "12H15M", "15M12H",
                        default27LongR1Bases, default27LongR2Bases,
                        default27LongR1Bases.substring(12, 27), default27LongR2Bases.substring(0, 15),
                        SequenceUtil.reverseComplement(default27LongR1Bases.substring(0, 12)), default27LongR2Bases.substring(15, 27),
                        default27LongR1Qualities, default27LongR2Qualities,
                        default27LongR1Qualities.substring(12, 27), default27LongR2Qualities.substring(0, 15),
                        StringUtil.reverseString(default27LongR1Qualities.substring(0, 12)), default27LongR2Qualities.substring(15, 27)},


                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //     <ssssssssMMMMMMMMMMMMMMMMMMM
                //                  MMMMMMMMMMMMMMM---MMMsssssssss>
                //          should become
                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //     <HHHHHHHHHHHHMMMMMMMMMMMMMMM
                //                  MMMMMMMMMMMMMMMHHHHHHHHHHHH>
                {27, true, 14, 18, "8S19M", "15M3D3M9S", true, false, 18, 18, "12H15M", "15M12H",
                        default27LongR1Bases, default27LongR2Bases,
                        default27LongR1Bases.substring(12, 27), default27LongR2Bases.substring(0, 15),
                        SequenceUtil.reverseComplement(default27LongR1Bases.substring(0, 12)), default27LongR2Bases.substring(15, 27),
                        default27LongR1Qualities, default27LongR2Qualities,
                        default27LongR1Qualities.substring(12, 27), default27LongR2Qualities.substring(0, 15),
                        StringUtil.reverseString(default27LongR1Qualities.substring(0, 12)), default27LongR2Qualities.substring(15, 27)},


                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //     <ssssssssMMMMMMMMMMMMMMMMMMM
                //                  MMMMMMMMMMMMMM---MMMMsssssssss>
                //          should become
                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //     <HHHHHHHHHHHHMMMMMMMMMMMMMMM
                //                  MMMMMMMMMMMMMMSHHHHHHHHHHHH>
                {27, true, 14, 18, "8S19M", "14M3D4M9S", true, false, 18, 18, "12H15M", "14M1S12H",
                        default27LongR1Bases, default27LongR2Bases,
                        default27LongR1Bases.substring(12, 27), default27LongR2Bases.substring(0, 15),
                        SequenceUtil.reverseComplement(default27LongR1Bases.substring(0, 12)), default27LongR2Bases.substring(15, 27),
                        default27LongR1Qualities, default27LongR2Qualities,
                        default27LongR1Qualities.substring(12, 27), default27LongR2Qualities.substring(0, 15),
                        StringUtil.reverseString(default27LongR1Qualities.substring(0, 12)), default27LongR2Qualities.substring(15, 27)},


                // 123456789.123456789.123456789.12-3456789.123456789.123456789.123456789.123456789
                //     <ssssssssMMMMMMMMMMMMMMMMMMM
                //      MMMMMMMMMMMMMMMMMMMMMMMMMMM-> ## deletion right at the end
                //          should remain
                // 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //     <ssssssssMMMMMMMMMMMMMMMMMMM
                //      MMMMMMMMMMMMMMMMMMMMMMMMMMM->  ## deletion right at the end
                {27, true, 14, 6, "8S19M", "27M1D", true, false, 14, 6, "8S19M", "27M1D",
                        default27LongR1Bases, default27LongR2Bases, default27LongR1Bases, default27LongR2Bases, null, null,
                        default27LongR1Qualities, default27LongR2Qualities, default27LongR1Qualities, default27LongR2Qualities, null, null},


                //                  123456789.123456--789.123456789.123456789.123456789.123456789.123456789.123456789
                //                                     MMMMMMMMMMMMMMMMMMsssssssss>
                //                          <sssssMMiiMMMMMMMMMMMMMMMMMM
                //                                  ^^
                //   with an insertion of two bases here should become
                //                    123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //                                     MMMMMMMMMMMMMMMMMHHHHHHHHHH>
                //                          <HHHHHHHHHHMMMMMMMMMMMMMMMMM
                {27, true, 18, 15, "18M9S", "5S2M2I18M", false, true, 18, 18, "17M10H", "10H17M",
                        default27LongR1Bases, default27LongR2Bases,
                        default27LongR1Bases.substring(0, 17), default27LongR2Bases.substring(10, 27),
                        default27LongR1Bases.substring(17, 27), SequenceUtil.reverseComplement(default27LongR2Bases.substring(0, 10)),
                        default27LongR1Qualities, default27LongR2Qualities,
                        default27LongR1Qualities.substring(0, 17), default27LongR2Qualities.substring(10, 27),
                        default27LongR1Qualities.substring(17, 27), StringUtil.reverseString(default27LongR2Qualities.substring(0, 10))},


                //       123456789.123456789.123456789.1234--56789.123456789.123456789.123456789.123456789
                //            <sssssssMMMMMMMMMMMMMMMMMMMM
                //                        MMMMMMMMMMMMMMMMMiiMsssssss>
                //                                         ^^
                //         with an insertion of two bases here should become
                //       123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //            <HHHHHHHHHHHMMMMMMMMMMMMMMMM
                //                        MMMMMMMMMMMMMMMMHHHHHHHHHHH>
                {27, true, 14, 18, "7S20M", "17M2I1M7S", true, false, 18, 18, "11H16M", "16M11H",
                        default27LongR1Bases, default27LongR2Bases,
                        default27LongR1Bases.substring(11, 27), default27LongR2Bases.substring(0, 16),
                        SequenceUtil.reverseComplement(default27LongR1Bases.substring(0, 11)), default27LongR2Bases.substring(16, 27),
                        default27LongR1Qualities, default27LongR2Qualities,
                        default27LongR1Qualities.substring(11, 27), default27LongR2Qualities.substring(0, 16),
                        StringUtil.reverseString(default27LongR1Qualities.substring(0, 11)), default27LongR2Qualities.substring(16, 27)},


                //       123456789.123456789.123456789.123--456789.123456789.123456789.123456789.123456789
                //            <sssssssMMMMMMMMMMMMMMMMMMMM
                //                        MMMMMMMMMMMMMMMMiiMMsssssss>
                //                                        ^^
                //        with an insertion of two bases here should become
                //       123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //            <sssssssSSSSMMMMMMMMMMMMMMMM
                //                        MMMMMMMMMMMMMMMMSSSSsssssss>
                {27, true, 14, 18, "7S20M", "16M2I2M7S", true, false, 18, 18, "11H16M", "16M11H",
                        default27LongR1Bases, default27LongR2Bases,
                        default27LongR1Bases.substring(11, 27), default27LongR2Bases.substring(0, 16),
                        SequenceUtil.reverseComplement(default27LongR1Bases.substring(0, 11)), default27LongR2Bases.substring(16, 27),
                        default27LongR1Qualities, default27LongR2Qualities,
                        default27LongR1Qualities.substring(11, 27), default27LongR2Qualities.substring(0, 16),
                        StringUtil.reverseString(default27LongR1Qualities.substring(0, 11)), default27LongR2Qualities.substring(16, 27)},


                //       123456789.123456789.123456789.123-456789.123456789.123456789.123456789.123456789
                //            <sssssssMMMMMMMMMMMMMMMMMMMM
                //                        MMMMMMMMMMMMMMMMiMMMsssssss>
                //                                        ^
                //        with an insertion of one base here should become
                //       123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //            <HHHHHHHHHHHMMMMMMMMMMMMMMMM
                //                        MMMMMMMMMMMMMMMMHHHHHHHHHHH>
                {27, true, 14, 18, "7S20M", "16M1I3M7S", true, false, 18, 18, "11H16M", "16M11H",
                        default27LongR1Bases, default27LongR2Bases,
                        default27LongR1Bases.substring(11, 27), default27LongR2Bases.substring(0, 16),
                        SequenceUtil.reverseComplement(default27LongR1Bases.substring(0, 11)), default27LongR2Bases.substring(16, 27),
                        default27LongR1Qualities, default27LongR2Qualities,
                        default27LongR1Qualities.substring(11, 27), default27LongR2Qualities.substring(0, 16),
                        StringUtil.reverseString(default27LongR1Qualities.substring(0, 11)), default27LongR2Qualities.substring(16, 27)},


                //       123456789.123456789.123456789.1234--56789.123456789.123456789.123456789.123456789
                //            <sssssssMMMMMMMMMMMMMMMMMMMM
                //                        MMMMMMMMMMMMMMMMMiiMMssssss>
                //                                         ^^
                //         with an insertion of two bases here should become
                //       123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                //            <HHHHHHHHHHHMMMMMMMMMMMMMMMM
                //                        MMMMMMMMMMMMMMMMHHHHHHHHHHH>
                {27, true, 14, 18, "7S20M", "17M2I2M6S", true, false, 18, 18, "11H16M", "16M11H",
                        default27LongR1Bases, default27LongR2Bases,
                        default27LongR1Bases.substring(11, 27), default27LongR2Bases.substring(0, 16),
                        SequenceUtil.reverseComplement(default27LongR1Bases.substring(0, 11)), default27LongR2Bases.substring(16, 27),
                        default27LongR1Qualities, default27LongR2Qualities,
                        default27LongR1Qualities.substring(11, 27), default27LongR2Qualities.substring(0, 16),
                        StringUtil.reverseString(default27LongR1Qualities.substring(0, 11)), default27LongR2Qualities.substring(16, 27)},
        };
    }

    @DataProvider
    public Iterator<Object[]> overlapReadDataWithSwaps() {
        List<Object[]> tests = new LinkedList<>();

        for (Object[] inputs : overlapReadData()) {
            tests.add(inputs);
            final Object[] swappedInputs = new Object[inputs.length];
            swappedInputs[0] = inputs[0]; // read length
            swappedInputs[1] = inputs[1]; // hard_clip

            for (int i = 2; i < inputs.length; i += 2) {
                swappedInputs[i] = inputs[i + 1];
                swappedInputs[i + 1] = inputs[i];
            }
            tests.add(swappedInputs);
        }
        return tests.iterator();
    }

    @Test(dataProvider = "overlapReadDataWithSwaps")
    public void testOverlappedReadClipping(
            final int readLength, final boolean hardClipOverlappingReads,
            final int start1, final int start2,
            final String cigar1, final String cigar2,
            final boolean strand1, final boolean strand2,
            final int r1ExpectedAlignmentStart, final int r2ExpectedAlignmentStart,
            final String expectedR1Cigar, final String expectedR2Cigar,
            final String read1Bases, final String read2Bases,
            final String expectedR1Bases, final String expectedR2Bases,
            final String expectedR1ClippedBases, final String expectedR2ClippedBases,
            final String read1Qualities, final String read2Qualities,
            final String expectedR1Qualities, final String expectedR2Qualities,
            final String expectedR1ClippedQualities, final String expectedR2ClippedQualities) {

        final SAMRecordSetBuilder set = new SAMRecordSetBuilder();
        set.setReadLength(readLength);
        final List<SAMRecord> recs = set.addPair("q1", 0, start1, start2, false, false, cigar1, cigar2, strand1, strand2, 30);
        final SAMRecord r1 = recs.get(0);
        final SAMRecord r2 = recs.get(1);

        r1.setReadBases(StringUtil.stringToBytes(read1Bases));
        r2.setReadBases(StringUtil.stringToBytes(read2Bases));

        r1.setBaseQualities(SAMUtils.fastqToPhred(read1Qualities));
        r2.setBaseQualities(SAMUtils.fastqToPhred(read2Qualities));

        AbstractAlignmentMerger.clipForOverlappingReads(r1, r2, hardClipOverlappingReads);
        Assert.assertEquals(r1.getAlignmentStart(), r1ExpectedAlignmentStart, "r1 POS");
        Assert.assertEquals(r1.getCigarString(), expectedR1Cigar, "r1 CIGAR");
        Assert.assertEquals(r2.getAlignmentStart(), r2ExpectedAlignmentStart, "r2 POS");
        Assert.assertEquals(r2.getCigarString(), expectedR2Cigar, "r2 CIGAR");
        Assert.assertEquals(r1.getReadString(), expectedR1Bases, "r1 BASES");
        Assert.assertEquals(r2.getReadString(), expectedR2Bases, "r1 BASES");
        Assert.assertEquals(SAMUtils.phredToFastq(r1.getBaseQualities()), expectedR1Qualities, "r1 QUAL");
        Assert.assertEquals(SAMUtils.phredToFastq(r2.getBaseQualities()), expectedR2Qualities, "r2 QUAL");

        Assert.assertEquals(r1.getAttribute(AbstractAlignmentMerger.HARD_CLIPPED_BASES_TAG), expectedR1ClippedBases, "r1 CLIPPED BASES");
        Assert.assertEquals(r2.getAttribute(AbstractAlignmentMerger.HARD_CLIPPED_BASES_TAG), expectedR2ClippedBases, "r2 CLIPPED BASES");

        Assert.assertEquals(r1.getAttribute(AbstractAlignmentMerger.HARD_CLIPPED_BASE_QUALITIES_TAG), expectedR1ClippedQualities, "r1 CLIPPED QUALS");
        Assert.assertEquals(r2.getAttribute(AbstractAlignmentMerger.HARD_CLIPPED_BASE_QUALITIES_TAG), expectedR2ClippedQualities, "r2 CLIPPED QUALS");
    }


    @Test
    public void testUnmapBacterialContamination() throws IOException {
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.queryname);
        final SAMFileHeader header = builder.getHeader();
        final SAMFileHeader.SortOrder sortOrder = header.getSortOrder();
        final SAMFileHeader newHeader = SAMRecordSetBuilder.makeDefaultHeader(sortOrder, 100000, true);
        builder.setHeader(newHeader);

        final File reference = File.createTempFile("reference", ".fasta");
        reference.deleteOnExit();

        builder.writeRandomReference(reference.toPath());

        builder.addPair("overlappingpair", 0, 500, 500, false, false, "20S20M60S", "20S20M60M", true, false, 45);
        builder.addPair("overlappingpairFirstAligned", 0, 500, 500, false, true, "20S20M60S", null, true, false, 45);
        builder.addPair("overlappingpairSecondAligned", 0, 500, 500, true, false, null, "20S20M60S", true, false, 45);
        builder.addPair("overlappingpairFirstAlignedB", 0, 500, 500, false, true, "20S20M60S", null, false, true, 45);
        builder.addPair("overlappingpairSecondAlignedB", 0, 500, 500, true, false, null, "20S20M60S", false, true, 45);

//        builder.addFrag("frag",1,500,false,false,"20S20M60S",null, 45);
//        builder.addFrag("frag2",1,500,true,false,"20S20M60S",null, 45);
//        builder.addFrag("frag3",1,500,false,false,"20S20M60S",null, 45);
//        builder.addFrag("frag4",1,500,true,false,"20S20M60S",null, 45);

        final PicardHtsPath file = PicardHtsPath.fromPath(newTempSamFile("aligned").toPath());

        try (SAMFileWriter writer = new SAMFileWriterFactory().makeWriter(builder.getHeader(), true, file.toPath(), (Path) null)) {
            builder.getRecords().forEach(writer::addAlignment);
        }

        final RevertSam revertSam = new RevertSam();

        revertSam.INPUT = file;
        final File fileUnaligned = newTempSamFile("unaligned");
        revertSam.OUTPUT = fileUnaligned;

        revertSam.SANITIZE = false;
        revertSam.REMOVE_ALIGNMENT_INFORMATION = true;
        revertSam.REMOVE_DUPLICATE_INFORMATION = true;

        revertSam.SORT_ORDER = SAMFileHeader.SortOrder.queryname;

        Assert.assertEquals(revertSam.doWork(), 0);

        MergeBamAlignment mergeBamAlignment = new MergeBamAlignment();

        mergeBamAlignment.ALIGNED_BAM = Collections.singletonList(file.toPath().toFile()); // TODO update to use Path when MergeBamAlignment is updated to use Path
        mergeBamAlignment.UNMAPPED_BAM = fileUnaligned;
        mergeBamAlignment.UNMAP_CONTAMINANT_READS = true;

        //yuck!
        final RequiredReferenceArgumentCollection requiredReferenceArgumentCollection = new RequiredReferenceArgumentCollection();
        requiredReferenceArgumentCollection.REFERENCE_SEQUENCE = new PicardHtsPath(reference.getAbsolutePath());
        mergeBamAlignment.referenceSequence = requiredReferenceArgumentCollection;

        final File fileMerged = newTempSamFile("merged");

        mergeBamAlignment.OUTPUT = fileMerged;

        // merge file with itself.
        Assert.assertEquals(mergeBamAlignment.doWork(), 0);

        // check that all reads have been unmapped due to bacterial contamination as needed.
        try (SamReader mergedReader = SamReaderFactory.makeDefault().open(fileMerged)) {
            for (SAMRecord mergedRecord : mergedReader) {
                Assert.assertTrue(mergedRecord.getReadUnmappedFlag(), mergedRecord.toString());
                Assert.assertTrue(!mergedRecord.getReadPairedFlag() || mergedRecord.getMateUnmappedFlag(), mergedRecord.toString());
            }
        }
    }

    @Override
    public String getCommandLineProgramName() {
        return this.getClass().getSimpleName();
    }


    private static File newTempSamFile(final String filename) throws IOException {
        final File file = File.createTempFile(filename, ".sam");
        file.deleteOnExit();
        return file;
    }

    @DataProvider(name = "readPositionIgnoringSoftClips")
    public Object[][] readPositionIgnoringSoftClips() {
        return new Object[][]{
                {"26S58M62S", 3688, 3827, 0}, // This is from the read that made us aware of a bug
                {"26S58M62S", 3688, 3665, 4},
                {"26S58M62S", 3688, 3660, 0}, // Before soft clip
                {"10S100M2S", 5, 10, 16},
                {"10S100M2S", 5, 3, 9},
                {"10S100M2S", 10, 12, 13},
                {"10S100M2S", 5, 107, 0}
        };
    }

    @Test(dataProvider = "readPositionIgnoringSoftClips")
    public void testGetReadPositionIgnoringSoftClips(final String cigarString, final int startPosition, final int queryPosition, final int expectedReadPosititon) {
        final SAMFileHeader newHeader = SAMRecordSetBuilder.makeDefaultHeader(SAMFileHeader.SortOrder.queryname, 100000, false);
        final SAMRecord rec = new SAMRecord(newHeader);

        rec.setCigarString(cigarString);
        rec.setAlignmentStart(startPosition);

        final int readPosition = AbstractAlignmentMerger.getReadPositionAtReferencePositionIgnoreSoftClips(rec, queryPosition);

        Assert.assertEquals(readPosition, expectedReadPosititon);
    }


    @DataProvider
    public Object[][] referencePositionAndReadPositions() {
        return new Object[][]{
                {1, 0, false}, {1, 0, true},
                {2, 0, false}, {2, 0, true},
                {3, 0, false}, {3, 0, true},
                {4, 1, false}, {4, 1, true},
                {5, 2, false}, {5, 2, true},
                {6, 3, false}, {6, 3, true},
                {7, 4, false}, {7, 4, true},
                {8, 5, false}, {8, 5, true},
                {9, 6, false}, {9, 6, true},
                {10, 7, false}, {10, 7, true},
                {11, 8, false}, {11, 8, true},
                {12, 8, false}, {12, 8, true},
                {13, 8, false}, {13, 8, true},
                {14, 8, false}, {14, 8, true},
                {15, 8, false}, {15, 8, true},
                {16, 9, false}, {16, 9, true},
                {17, 10, false}, {17, 10, true},
                {18, 11, false}, {18, 11, true},
                {19, 12, false}, {19, 12, true},
                {20, 17, false}, {20, 17, true},
                {21, 18, false}, {21, 18, true},
                {22, 19, false}, {22, 19, true},
                {23, 20, false}, {23, 20, true},
                {24, 20, false}, {24, 20, true},
                {25, 20, false}, {25, 20, true},
                {26, 20, false}, {26, 20, true},
                {27, 20, false}, {27, 20, true},
                {28, 21, false}, {28, 21, true},
                {29, 22, false}, {29, 22, true},
                {30, 23, false}, {30, 23, true},
                {31, 24, false}, {31, 24, true},
                {32, 0, false}, {32, 0, true},
        };
    }

    @Test(dataProvider = "referencePositionAndReadPositions")
    public void testGetReadPositionAtReferencePositionIgnoreSoftClips(final int refPos, final int readPos, final boolean negativeStrand) {
        final SAMRecordSetBuilder set = new SAMRecordSetBuilder();
        //REF   123456789.123456789----.123456789.123456789.
        //         ....||||----||||....||||----....
        //         SSSSMMMMDDDDMMMMIIIIMMMMDDDDSSSS>
        //READ     12345678----9.123456789.----12345678
        final SAMRecord samRecord = set.addFrag("test", 0, 8, negativeStrand, false, "4S4M4D4M4I4M4D4S", "null", 30);

        Assert.assertEquals(AbstractAlignmentMerger.getReadPositionAtReferencePositionIgnoreSoftClips(samRecord, refPos), readPos);
    }
}

