package net.sf.picard.illumina.parser;

import net.sf.picard.PicardException;

import java.io.File;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static net.sf.picard.util.IlluminaUtilTest.iToB;

public class QSeqTdUtil {

    private static Map<Integer, List<Integer>> laneToReadNos = new HashMap<Integer, List<Integer>>();
    private static Map<String, File> fnPrefixToQSeqFile = new HashMap<String, File>();
    private static Map<String, Integer> fnPrefixToNumReads = new HashMap<String, Integer>();
    private static Map<String, List<QseqReadData>> fnPrefixToQSeqReadData = new HashMap<String, List<QseqReadData>>();


    public static final String s_1_1 = "s_1_1";
    public static final String s_1_2 = "s_1_2";
    public static final String s_1_1_0001 = "s_1_1_0001";
    public static final String s_1_1_0002 = "s_1_1_0002";
    public static final String s_1_1_0003 = "s_1_1_0003";
    public static final String s_1_2_0001 = "s_1_2_0001";
    public static final String s_1_2_0002 = "s_1_2_0002";
    public static final String s_1_2_0003 = "s_1_2_0003";

    public static final String s_5_1 = "s_5_1";
    public static final String s_5_1_0001 = "s_5_1_0001";

    public static final String s_6_1 = "s_6_1";
    public static final String s_6_1_0001 = "s_6_1_0001";
    public static final String s_6_2_0001 = "s_6_2_0001";
    public static final String s_6_3_0001 = "s_6_3_0001";

    public static final String s_8_1_0001 = "s_8_1_0001";
    public static final String s_8_2_0001 = "s_8_2_0001";

    public static final byte A = (byte)65;
    public static final byte C = (byte)67;
    public static final byte G = (byte)71;
    public static final byte T = (byte)84;
    public static final byte P = (byte)46; //dot
    public static final String NO_TILE_PREFIX_STRING = "^s_(\\d)_(\\d)$";
    public static final Pattern NO_TILE_PREFIX_PATTERN = Pattern.compile(NO_TILE_PREFIX_STRING);
    public static final String PREFIX_PATTERN_STRING = "^s_(\\d)_(\\d)_(\\d{4})$";
    public static final Pattern PREFIX_PATTERN = Pattern.compile(PREFIX_PATTERN_STRING);
    private static final int LANE = 0;
    private static final int READ = 1;
    private static final int TILE = 2;

    public static final File PRIMARY_TESTDATA_DIR = new File("testdata/net/sf/picard/illumina/IlluminaTests/BasecallsDir");

    public static int getLaneFileSize(int lane) {
        switch(lane) {
            case 7:
                return 200;
            default:
                return 20;
        }
    }

    public static String makeFnPrefix(int lane, int end, int tile) {
        String prefix = "s_" + lane + "_" + end + "_";
        for(int i = 0; i < 4 - String.valueOf(tile).length(); i++) {
            prefix += "0";
        }
        prefix += tile;
        return prefix;
    }

    private static int[] varsFromPrefixAndRegex(final String prefix, final Pattern pattern, final String errMsg, int tokensExpected) {
        Matcher matcher = pattern.matcher(prefix);
        if(!matcher.find()) {
            throw new PicardException(errMsg);
        }

        int [] out = new int[tokensExpected];
        for(int i = 0; i < tokensExpected; i++) {
            out[i] = Integer.parseInt(matcher.group(i+1));
        }
        return out;
    }

    public static int[] varsFromNoTilePrefix(final String prefix) {
        return varsFromPrefixAndRegex(prefix, NO_TILE_PREFIX_PATTERN, "Couldn't find NO_TILE_PREFIX_PATTERN: " + NO_TILE_PREFIX_STRING + " in " + prefix, 2);
    }

    public static int[] varsFromPrefix(final String prefix) {
        return varsFromPrefixAndRegex(prefix, PREFIX_PATTERN, "Couldn't find PREFIX_PATTERN: " + PREFIX_PATTERN_STRING + " in " + prefix, 3);
    }

    public static void addPrefToFileMap(final String pref) {
        fnPrefixToQSeqFile.put(pref, new File(PRIMARY_TESTDATA_DIR, pref + "_qseq.txt"));
    }

    public static void addFileSize(final String pref, int numReads) {
        fnPrefixToNumReads.put(pref, numReads);
    }

    public static int getFileSize(final String pref) {
        return fnPrefixToNumReads.get(pref);
    }

    static {
        //lane 1
        laneToReadNos.put(1, Arrays.asList(0, 9, 19));

        String pref = makeFnPrefix(1,1,1);
        addPrefToFileMap(pref);
        addFileSize(pref, 20);
        fnPrefixToQSeqReadData.put(pref, Arrays.asList(
        makeQSeqReadData(new byte[][]{new byte[]{G, P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P, C, P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P, T,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,  T, P,P,P,P,P,P,P,P,P,P}},
                         new int [][]{new int[] {28,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 28,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 27,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 27,4,4,4,4,4,4,4,4,4,4}},
                         false, 1793, 1420, 1, 1),

        makeQSeqReadData(new byte[][]{new byte[]{G, P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,  C, P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P, C,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,  C, P,P,P,P,P,P,P,P,P,P}},
                         new int [][]{new int[] {25,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  27,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 8,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  27,4,4,4,4,4,4,4,4,4,4}},
                         false, 1793, 1718, 1, 1),

        makeQSeqReadData(new byte[][]{new byte[]{A, P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,  C, P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P, G,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,  G, P,P,P,P,P,P,P,P,P,P}},
                         new int [][]{new int[] {7, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  26,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 16,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 26,4,4,4,4,4,4,4,4,4,4}},
                         false, 1793, 1011, 1, 1)
        ));

        pref = makeFnPrefix(1,1,2);
        addPrefToFileMap(pref);
        addFileSize(pref, 20);
        fnPrefixToQSeqReadData.put(pref, Arrays.asList(
        makeQSeqReadData(new byte[][]{new byte[]{G, P, T, C, C, P,P,P,P,P,  P,P,A,P,T, C, C, C, C,G,  T, G, G, A, T, G, A, A,A, A,  T, G, C, T, G, T, A, C, A, T,      G, T, T, C, T, T, G, T, T, G,    A, T, A, A, C, A, A,T, G, G, G,   C, A, A, G, G, T, C, T, A,   G, G, A, C, A, G}},
                         new int [][]{new int[] {15,4, 26,34,27,4,4,4,4,4,  4,4,8,4,24,34,34,31,7,29, 30,34,31,16,28,32,22,7,15,15, 28,33,34,32,33,24,21,32,29,27,     32,30,31,34,32,32,33,31,32,32,   20,29,22,32,33,16,6,21,29,22,29,  33,24,13,24,33,22,33,20,8,   22,26,20,33,27,29}},
                         false, 1793, 1893, 1, 2),

        makeQSeqReadData(new byte[][]{new byte[]{C, P,G, A,G, P,P,P,P,P,  P,P,T,P,G, G,G, T, A, G,  G, A, A, A, C, A, C, A, G, C,   C, T, T, G, C, T,C, C, A, C,   A, G, C,G, C, A, C, T, G, T,   C, A, G, C, A, A,G, A, C,G,    C, T, C, C, T, T, T, C, T, T,   C,G,G, A, G, A}},
                         new int [][]{new int[] {30,4,24,7,24,4,4,4,4,4,  4,4,4,4,22,8,24,13,13,27, 25,25,29,18,29,16,26,22,22,27,  33,15,20,30,30,8,33,29,15,28,  16,17,7,17,32,21,28,19,25,19,  24,20,21,30,14,4,20,12,4,12,   26,16,26,21,20,19,17,26,14,12,  5,7,17,21,11,7}},
                         false, 1793, 666, 1, 2),

        makeQSeqReadData(new byte[][]{new byte[]{C, P,T, C, C, P,P,P,P,P,   P,P,T, P,T, A,T, T, T, A,  C, T, A, T, T, T, T, C, T, G,   A, T, T, T, T, T, A, A, A, A,   T, G, A, C, A, G, T, G, G, C,   A, A,T, T, A, C, C, A, T, T,  T, A, T, A, C, T, G, T, G, T,   T, A, T, T, T, G}},
                         new int [][]{new int[] {30,4,23,34,30,4,4,4,4,4,   4,4,24,4,21,7,28,30,25,7,  28,31,33,32,32,31,31,33,31,32,  27,27,27,27,28,31,28,26,20,30,  31,24,10,30,32,31,25,30,24,31,  20,5,16,28,29,32,30,10,25,21, 21,21,23,28,32,24,30,22,30,24,  30,28,27,24,28,21}},
                         false, 1794, 490, 1, 2)
        ));

        pref = makeFnPrefix(1,1,3);
        addPrefToFileMap(pref);
        addFileSize(pref, 20);
        fnPrefixToQSeqReadData.put(pref, Arrays.asList(
        makeQSeqReadData(new byte[][]{new byte[]{G, P,T, A,G, P,P,P,P,P,    P,P,C,P,G, G, A, C,T, C,     T, G, G, C, G, T, C,A, C, C,    T, T,T, G, G,C, G, C, T, G,     A,G, C,G, C, C, C, C, A,G,     G, C, C, C, G,C, C, A, G, C,     C, C, G, C, C,C, G, C, C, C,     A, C,T,G, C, C}},
                         new int [][]{new int[] {27,4,26,7,13,4,4,4,4,4,    4,4,7,4,13,21,13,7,26,33,    31,20,20,10,24,29,9,23,31,7,    29,9,31,26,7,13,26,31,22,20,    6,20,6,26,26,31,26,26,9,19,    27,20,27,17,9,31,26,14,29,27,    31,20,29,20,7,16,21,26,16,26,    22,8,6,11,13,13}},
                         false, 1793, 1282, 1, 3),

        makeQSeqReadData(new byte[][]{new byte[]{G, P,G, G, G, P,P,P,P,P,    P,P,T, P,C, C, A, T, T, G,     G, T, T, T, C, T, G, A, A, A,     G, T, A, T, T, C, A,C, A, T,     C, A, T, T, T, G, G, G, A, T,     A, C, C, A, G, A, T, A,G, C,     T, C, A, A, T, A, C, T, C, T,     C, T, G, A, G, T}},
                         new int [][]{new int[] {23,4,29,33,27,4,4,4,4,4,    4,4,25,4,30,33,23,33,33,33,    33,30,33,33,34,32,34,30,31,21,    31,16,26,32,33,31,9,31,27,29,    33,24,29,30,31,31,33,30,24,27,    19,33,33,23,31,23,17,8,27,34,    31,34,32,25,30,23,33,31,34,31,    34,31,30,14,28,6}},
                         false, 1793, 1709, 1, 3),

        makeQSeqReadData(new byte[][]{new byte[]{C, P,T, A, T, P,P,P,P,P,    P,P,C, P,P,C, G, G, G, T,     A, C, C, A, C, A, G, T, T, G,     A, G, G, A, C, T, G, A, C, A,     T, T, C, T, G, A, A, C, C, C,     T, G, A, T, G, T, T, T, C, T,     A, A, A, G, A, A, A, C, G, A,    C, A, G, T, A, T}},
                         new int [][]{new int[] {30,4,21,31,23,4,4,4,4,4,    4,4,30,4,4,30,31,30,31,21,    32,32,32,30,30,21,27,29,29,31,    21,14,22,30,32,30,30,21,29,27,    27,23,28,31,17,20,17,26,29,31,    23,19,27,24,29,20,23,31,30,22,    13,12,21,26,22,23,26,27,19,6,    24,12,27,25,21,23}},
                         false, 1793, 456, 1, 3)
        ));

        pref = makeFnPrefix(1,2,1);
        addPrefToFileMap(pref);
        addFileSize(pref, 20);
        fnPrefixToQSeqReadData.put(pref, Arrays.asList(
        makeQSeqReadData(new byte[][]{new byte[]{C, P,P,P,A,P,P,P,P,P,   P,P,P,P,P,P,P,P,P,P,  P,A, P,P,P,P,P,P,P,P,  P,P,P,P,P,P,P,P,P,P,  P,P,P,C, P,P,P,P,P,P,  P,P,P,P,P,P,P,P,P,P,  P,P,P,P,P,T, P,P,P,P,  P,P,P,P,P,P}},
                         new int [][]{new int[] {30,4,4,4,9,4,4,4,4,4,   4,4,4,4,4,4,4,4,4,4,  4,11,4,4,4,4,4,4,4,4,  4,4,4,4,4,4,4,4,4,4,  4,4,4,27,4,4,4,4,4,4,  4,4,4,4,4,4,4,4,4,4,  4,4,4,4,4,13,4,4,4,4,  4,4,4,4,4,4}},
                         false, 1793, 1420, 1, 1),

        makeQSeqReadData(new byte[][]{new byte[]{A,P,P,P,A,P,P,P,P,P,    P,P,P,P,P,P,P,P,P,P,    P,A, P,P,P,P,P,P,P,P,    P,P,P,P,P,P,P,P,P,P,    P,P,P,A,P,P,P,P,P,P,    P,P,P,P,P,P,P,P,P,P,    P,P,P,P,P,G,P,P,P,P,    P,P,P,P,P,P}},
                         new int [][]{new int[] {7,4,4,4,4,4,4,4,4,4,    4,4,4,4,4,4,4,4,4,4,    4,14,4,4,4,4,4,4,4,4,    4,4,4,4,4,4,4,4,4,4,    4,4,4,6,4,4,4,4,4,4,    4,4,4,4,4,4,4,4,4,4,    4,4,4,4,4,7,4,4,4,4,    4,4,4,4,4,4}},
                         false, 1793, 1718, 1, 1),

        makeQSeqReadData(new byte[][]{new byte[]{A, P,P,P,T,P,P,P,P,P,    P,P,P,P,P,P,P,P,P,P,    P,T,P,P,P,P,P,P,P,P,    P,P,P,P,P,P,P,P,P,P,    P,P,P,C, P,P,P,P,P,P,    P,P,P,P,P,P,P,P,P,P,    P,P,P,P,P,C,P,P,P,P,    P,P,P,P,P,P}},
                         new int [][]{new int[] {18,4,4,4,4,4,4,4,4,4,    4,4,4,4,4,4,4,4,4,4,    4,7,4,4,4,4,4,4,4,4,    4,4,4,4,4,4,4,4,4,4,    4,4,4,30,4,4,4,4,4,4,    4,4,4,4,4,4,4,4,4,4,    4,4,4,4,4,9,4,4,4,4,    4,4,4,4,4,4}},
                         false, 1793, 1011, 1, 1)
        ));

        pref = makeFnPrefix(1,2,2);
        addPrefToFileMap(pref);
        addFileSize(pref, 20);
        fnPrefixToQSeqReadData.put(pref, Arrays.asList(
        makeQSeqReadData(new byte[][]{new byte[]{G, G, G, C, A, G, G, G, G, G,     C,T, T, T, T, T, C, T, C, A,     C,T, C, T, C, T, C, T, C, A,     T, A,T, C, T, T, C, T, A,G,     G, G, T, A, A, C, T, A,C, A,     T, G, A, A,C, A, C,A, C, G,     C, T, T, C, T, C, T, C, C, C,     C,T, T, C,C,G}},
                         new int [][]{new int[] {32,33,33,33,19,28,33,33,33,25,    8,28,33,27,32,32,34,30,34,29,    8,30,34,31,34,30,34,33,33,18,    26,8,27,34,31,33,34,24,7,29,    30,32,18,31,12,31,28,8,31,21,    31,28,19,6,30,19,7,19,17,20,    34,24,28,33,27,34,31,34,33,31,    7,15,22,7,7,25}},
                         false, 1793, 1893, 1, 2),

        makeQSeqReadData(new byte[][]{new byte[]{A, T, T, T, T, T, G, C, T, T,     T, C, C, T, T, G, G, C, C, C,     C, C, A, C, C, A, A, T, T, T,     A, T, A, C, A, T, C, T, C, C,     A, T, T, T, T, C, C, G, A,C,     C, T, C, T, G, G,A, C, T, A,     A, C, T, G, C,T,T, G, C, T,    C, A, G, C, A,C}},
                         new int [][]{new int[] {21,26,33,32,34,30,32,33,31,31,    29,31,23,26,31,29,31,33,33,33,    33,32,25,31,33,23,28,30,26,16,    29,21,31,32,32,27,33,19,32,33,    25,24,19,27,26,28,11,15,6,27,    31,27,32,23,22,8,22,33,28,22,    18,31,26,18,8,6,14,21,24,7,    27,11,21,22,6,4}},
                         false, 1793, 1909, 1, 2),

        makeQSeqReadData(new byte[][]{new byte[]{A, T, A, A, A, C, T, T, G, A,     A, A, A, T, A, A, T, T, T, T,     C, T, A, T, G, A, T, A, C, A,     G, C, T, T, T, C, A, G, G, T,     A, G, A,A, A, A, A, T, G, A,     A, T, T, T, T, C, G, T, C, G,     T, G, T, T, T, A, A, C, A, A,    T, G, T, T, G, T}},
                         new int [][]{new int[] {32,33,33,34,33,34,32,32,34,31,    26,20,31,31,22,22,33,33,33,30,    34,34,33,33,33,30,31,33,34,33,    33,34,32,33,31,31,18,27,32,21,    29,29,8,24,24,30,15,27,33,31,    21,30,29,31,31,33,33,21,32,24,    29,33,20,32,27,31,21,27,4,10,    29,31,27,27,25,15}},
                         false, 1794, 490, 1, 2)
        ));

        pref = makeFnPrefix(1,2,3);
        addPrefToFileMap(pref);
        addFileSize(pref, 20);
        fnPrefixToQSeqReadData.put(pref, Arrays.asList(
        makeQSeqReadData(new byte[][]{new byte[]{A,A,G, C, A, C, A, C,C, A,     A, G, C, T, T, C, C, C, G, C,     A, G, C, T, G, C,G, G,C, G,    C, C, C, C, C, G,C, T, C, C,     A,A,G,G, C, C, A,G, C,A,     G, C, T, G, C, C, C, C, C, T,    C, T, G, C, C, C,C, C, A,C,    C, C, C, C, T, C}},
                         new int [][]{new int[] {7,9,20,22,24,13,16,7,13,14,    14,29,33,31,12,16,22,10,28,13,    13,14,31,24,28,9,28,7,13,7,    13,26,27,13,16,6,29,26,27,25,    6,6,6,22,28,26,6,13,6,25,    28,28,25,20,20,16,24,18,23,9,    22,13,22,26,25,4,22,25,6,6,    26,31,26,21,15,24}},
                         false, 1793, 1282, 1, 3),

        makeQSeqReadData(new byte[][]{new byte[]{A, G, C, C,A, C, T, G, T, C,     C, C, G, T, G, T, A, T, A, A,     C, T, T, G, G, C, A, T, T, A,     G, A, G, C, A, C, C, A,G, G,     T, C, T, G, T, T, G, G, A, T,     G, G, T,G, G, T, G, G, C, A,     G, G, C, C, C, A,A, T,T,A,     A,T, T, T, T,T}},
                         new int [][]{new int[] {11,30,31,8,22,31,33,31,21,33,    34,34,34,26,34,30,31,32,26,10,    31,32,26,33,33,33,28,26,33,32,    31,21,31,33,21,33,31,7,29,33,    23,34,32,32,21,27,33,30,21,24,    33,30,7,33,28,12,25,31,29,12,    22,31,23,20,25,4,13,8,8,18,    5,22,18,12,4,8}},
                         false, 1793, 1709, 1, 3),

        makeQSeqReadData(new byte[][]{new byte[]{T, C, C, A, T, C, C, A, C, T,     T, C, C, C, T, G, A, G, C, C,     T, C, A, G, A, A, A, A, G, G,     G, C, A, A, G, G, C, A, T, G,     G, C, T, C, A, C, A, T, A, C,     T, C, T, C, A,G, C, C, A, C,     G, G, C, C, T,G,G, C, C,T,     G, C, T, G, C, C}},
                         new int [][]{new int[] {33,33,33,27,32,33,33,31,33,33,    33,33,33,33,28,31,31,32,33,33,    30,33,20,31,22,22,22,28,26,26,    32,33,32,24,32,26,33,30,28,29,    27,33,33,31,21,31,32,32,30,33,    29,33,24,29,9,29,32,32,24,32,    20,18,30,29,7,4,23,24,7,13,    24,29,26,27,25,7}},
                         false, 1793, 456, 1, 3)
        ));

        laneToReadNos.put(5, Arrays.asList(0, 9, 19));
        pref = makeFnPrefix(5,1,1);
        addPrefToFileMap(pref);
        addFileSize(pref, 20);
        fnPrefixToQSeqReadData.put(pref, Arrays.asList(
        makeQSeqReadData(new byte[][]{new byte[]{G, A, C, T, T, T, G, G, G, A,     A, G, G, G, T, C, A, T, T, A,     C, T, G, C, C, C, T, T, G, T,     A, G, A, A, A, G, A, A, C, A,     C, C, T, C, A, T, G, T, T, C,     C, T, T, A, T, C, G, A, G, A,     G, C, G, G, C, C, G, C, T, G,     C, T, G, A,T,C}},
                         new int [][]{new int[] {23,27,32,34,34,34,34,31,34,33,    19,28,32,31,28,34,34,33,34,34,    33,23,18,32,34,34,33,32,32,33,    34,31,34,34,34,34,34,33,34,34,    34,33,33,34,34,30,30,30,28,33,    33,31,33,34,32,33,31,32,27,32,    22,33,19,20,27,30,19,23,20,24,    23,14,25,5,8,13}},
                         true, 0, 357, 5, 1),

        makeQSeqReadData(new byte[][]{new byte[]{C, A, A, A, A, C,A, A, C, T,     C, A, G, T, T, T, G, T, T, C,     C, A, A, A, A, C, A, A, T, G,     T, G, A, G, T, T, C, C, C, A,     G, A, T, T, T, A, G, C, C, T,     T, G, T, C, T, T, A, A, T, A,     T, A, T, T, G, A, C, C, T, T,     A, G, T, T, C, C}},
                         new int [][]{new int[] {33,32,32,33,21,9,29,30,31,17,    31,33,33,33,33,34,22,28,34,32,    27,32,33,34,34,31,19,33,33,32,    33,34,29,34,33,34,34,31,25,32,    31,33,33,29,33,33,32,30,33,33,    34,34,24,23,33,28,32,33,33,33,    34,33,33,33,28,26,25,10,17,28,    22,25,30,19,14,10}},
                         false, 0, 1484, 5, 1),

        makeQSeqReadData(new byte[][]{new byte[]{C, A, C, A, C, A, C, A, C, A,     C, A, C, A, C, A, C, A, C, A,     C, A, C, C, A, C, C, T, T, T,     T, G, G, C, T, T, A, T, C, T,     G, C, A, C, G, C, G, G, C, C,    G, C, G, T, G, C, C, C,T, A,     C, C, C, T, A, C, C, C, C, A,    T,G, G,G, A, T}},
                         new int [][]{new int[] {33,31,33,33,30,28,18,33,28,32,    33,33,24,33,31,33,33,31,33,33,    33,33,33,33,33,27,30,24,30,32,    32,22,27,32,31,32,33,30,33,33,    15,32,32,30,31,19,10,21,5,12,    20,22,13,11,22,20,16,6,15,11,    10,10,13,14,11,20,18,18,10,5,    5,16,6,11,20,18}},
                         true, 0, 1250, 5, 1)
        ));

        laneToReadNos.put(6, Arrays.asList(0, 4, 9));
        pref = makeFnPrefix(6,1,1);
        addPrefToFileMap(pref);
        addFileSize(pref, 10);
        fnPrefixToQSeqReadData.put(pref, Arrays.asList(
        makeQSeqReadData(new byte[][]{new byte[]{T, A, G, A,G,A,T,G,G,C,    P,C,T,P,P,P,P,P,P,P,    P,P,T,P,P,P,P,P,P,C,    P,P,P,P,P,P,P,P,G,P,    P,P,T,C,C,A,G,A,C,C,    G,C,C,C,A,T,T,C,T,C,    T,G,C,C,T,G,C,C}},
                         new int [][]{new int[] {30,32,30,2,2,2,2,2,2,2,    2,2,2,2,2,2,2,2,2,2,    2,2,2,2,2,2,2,2,2,2,    2,2,2,2,2,2,2,2,2,2,    2,2,2,2,2,2,2,2,2,2,    2,2,2,2,2,2,2,2,2,2,    2,2,2,2,2,2,2,2}},
                         false, 4, 1969, 6, 1),

        makeQSeqReadData(new byte[][]{new byte[]{T, A, A, A, C, A,G,A,T,G,    P,T,T,P,P,P,P,P,P,P,    P,C,A,P,P,P,P,P,P,T,    P,P,P,P,P,P,P,P,A,P,    P,P,C,C,A,A,T,C,C,C,    T,A,A,T,C,T,C,C,A,G,    T,A,A,T,C,C,G,G}},
                         new int [][]{new int[] {30,34,34,32,31,2,2,2,2,2,    2,2,2,2,2,2,2,2,2,2,    2,2,2,2,2,2,2,2,2,2,    2,2,2,2,2,2,2,2,2,2,    2,2,2,2,2,2,2,2,2,2,    2,2,2,2,2,2,2,2,2,2,    2,2,2,2,2,2,2,2}},
                         false, 4, 1248, 6, 1),

        makeQSeqReadData(new byte[][]{new byte[]{C, A, A, C, T, C, T,T,G,T,    P,G,T,P,P,P,P,P,P,P,    P,G,T,P,P,P,P,P,P,A,    P,P,P,P,P,P,P,P,G,P,    P,P,A,A,T,A,T,A,T,T,    C,T,G,A,A,A,C,T,C,A,    G,C,A,A,T,G,T,T}},
                         new int [][]{new int[] {33,33,33,34,33,33,2,2,2,2,    2,2,2,2,2,2,2,2,2,2,    2,2,2,2,2,2,2,2,2,2,    2,2,2,2,2,2,2,2,2,2,    2,2,2,2,2,2,2,2,2,2,    2,2,2,2,2,2,2,2,2,2,    2,2,2,2,2,2,2,2}},
                         false, 4, 151, 6, 1)
        ));

        pref = makeFnPrefix(6,2,1);
        addPrefToFileMap(pref);
        addFileSize(pref, 10);
        fnPrefixToQSeqReadData.put(pref, Arrays.asList(
        makeQSeqReadData(new byte[][]{new byte[]{C, C, A,C, C, C, A,C}},
                         new int [][]{new int[] {31,26,6,26,30,29,2,2}},
                         false, 4, 1969, 6, 1),

        makeQSeqReadData(new byte[][]{new byte[]{G, C, A, C, C, C, G,A}},
                         new int [][]{new int[] {20,31,24,31,10,31,2,2}},
                         false, 4, 1248, 6, 1),

        makeQSeqReadData(new byte[][]{new byte[]{T, C, G, G, A, A, T, G}},
                         new int [][]{new int[] {33,34,25,32,34,34,31,30}},
                         false, 4, 151, 6, 1)
        ));

        pref = makeFnPrefix(6,3,1);
        addPrefToFileMap(pref);
        addFileSize(pref, 10);
        fnPrefixToQSeqReadData.put(pref, Arrays.asList(
        makeQSeqReadData(new byte[][]{new byte[]{C, C, T, C,T, A, A, T, C, C,     C, A, G,C, A, C, T,A, T, C,     C, G,A, G, A, C, C, A, A,A,     T,C, A, G, G,C, A, A, A, T,     C, A, C, T, T, G,A, A,G, T,     C, A, G, G, A, G, T, T, C, G,     A, G, A, C, C, A, G, C}},
                         new int [][]{new int[] {29,29,22,9,19,17,11,31,13,28,    32,13,8,24,28,26,6,13,33,16,    23,8,24,25,21,33,29,26,8,10,    7,33,21,12,7,16,24,20,18,19,    28,23,27,32,26,8,31,7,13,21,    33,31,27,29,13,29,16,20,21,26,    24,29,22,33,26,33,21,28}},
                         false, 4, 1969, 6, 1),

        makeQSeqReadData(new byte[][]{new byte[]{A,A, T, A, T, T, C, T, T, T,     T, A, A, G, G, T, C, T, C, T,     G, G, T, T, T, T,C, C, T, A,     G, G, C, A, G, A, G, G, A, C,     C, C, T, G, C, G, G, C, C, T,     T, C, C, G, C, A, G, T, G, T,     T, T, G, T, G, T, C, C}},
                         new int [][]{new int[] {9,21,17,27,23,16,31,33,33,33,    29,23,34,34,27,33,34,32,34,21,    26,31,33,34,34,9,34,34,34,27,    34,27,31,28,31,33,34,34,33,33,    34,33,33,33,34,17,31,33,34,21,    14,31,34,34,34,32,34,32,34,30,    32,31,12,13,25,21,32,34}},
                         false, 4, 1248, 6, 1),

        makeQSeqReadData(new byte[][]{new byte[]{T, A, A, C, T, T, T, C, A, G,     A, G, G, C, C, C, T, T, C, A,     G, G, A, G, G, C, C, C, T, G,     G, C, C, T, G, T,C, A, A, G,     T, A, C, C, T, T, T, A, C, A,     G, T, G, A, T, G, G, G, T, A,     T, A, G,A,C,T,T,T}},
                         new int [][]{new int[] {33,34,34,34,33,34,34,33,31,34,    34,33,34,34,31,29,19,29,33,34,    31,33,34,31,31,30,33,34,19,15,    18,25,24,30,18,6,27,33,33,32,    31,33,31,27,24,22,33,32,27,23,    21,14,32,33,30,31,31,33,28,26,    12,28,2,2,2,2,2,2}},
                         false, 4, 151, 6, 1)
        ));

        laneToReadNos.put(8, Arrays.asList(0, 9, 19));
        pref = makeFnPrefix(8,1,1);
        addPrefToFileMap(pref);
        addFileSize(pref, 20);
        fnPrefixToQSeqReadData.put(pref, Arrays.asList(
        makeQSeqReadData(new byte[][]{new byte[]{G, P,P,P}},
                         new int [][]{new int[] {28,4,4,4}},
                         false, 1793, 1420, 8, 1),

        makeQSeqReadData(new byte[][]{new byte[]{G, P,P,P}},
                         new int [][]{new int[] {25,4,4,4}},
                         false, 1793, 1718, 8, 1),

        makeQSeqReadData(new byte[][]{new byte[]{A, P,P,P}},
                         new int [][]{new int[] {7,4,4,4}},
                         false, 1793, 1011, 8, 1)
        ));

        pref = makeFnPrefix(8,2,1);
        addPrefToFileMap(pref);
        addFileSize(pref, 20);
        fnPrefixToQSeqReadData.put(pref, Arrays.asList(
        makeQSeqReadData(new byte[][]{new byte[]{C, P,P,P}},
                         new int [][]{new int[] {30,4,4,4}},
                         false, 1793, 1420, 8, 2),

        makeQSeqReadData(new byte[][]{new byte[]{A,P,P,P}},
                         new int [][]{new int[] {7,4,4,4}},
                         false, 1793, 1718, 8, 2),

        makeQSeqReadData(new byte[][]{new byte[]{A, P,P,P}},
                         new int [][]{new int[] {18,4,4,4}},
                         false, 1793, 1011, 8, 2)
        ));
    }

    private static QseqReadData makeQSeqReadData(final byte [][] bases, final int [][] qualities, final boolean pf,
                                                final int xCoord, final int yCoord, final int lane, final int tile) {
        final int [] outputLengths = new int[bases.length];
        for(int i = 0; i < outputLengths.length; i++) {
            outputLengths[i] = bases[i].length;
        }
        final QseqReadData qseqRd = new QseqReadData(outputLengths);
        for(int i = 0; i < bases.length; i++) {
            System.arraycopy(bases[i],           0, qseqRd.getBases()[i],  0, bases[i].length);
            System.arraycopy(iToB(qualities[i]), 0, qseqRd.getQualities()[i], 0, qualities[i].length);
        }

        qseqRd.setOrCheckPf(pf);
        qseqRd.setOrCheckXCoordinate(xCoord);
        qseqRd.setOrCheckYCoordinate(yCoord);
        qseqRd.setOrCheckLane(lane);
        qseqRd.setOrCheckTile(tile);

        return qseqRd;
    }

    private static void inMemorySplit(final byte[][] source, final byte[][] dest, final int dstOffset) {
        int srcIndex = 0;
        int dstIndex = 0;
        int srcPos = 0;
        int dstPos = dstOffset;

        while(dstPos > dest[dstIndex].length) {
            dstPos -= dest[dstIndex].length;
            ++dstIndex;
        }

        while(srcIndex < source.length) {
            if(dest[dstIndex].length > 0) {
                final int ln = Math.min(source[srcIndex].length - srcPos, dest[dstIndex].length - dstPos);
                System.arraycopy(source[srcIndex], srcPos, dest[dstIndex], dstPos, ln);

                srcPos += ln;
                dstPos += ln;
            }

            if(srcPos >= source[srcIndex].length) {
                srcPos = 0;
                ++srcIndex;
            }

            if(dstPos >= dest[dstIndex].length) {
                dstPos = 0;
                ++dstIndex;
            }
        }
    }

    private static void fillArrays(byte [][] array, byte value, int fillLength) {
        int arrIndex = 0;
        while(fillLength > 0) {
            final int filledThisRound = Math.min(array[arrIndex].length, fillLength);
            Arrays.fill(array[arrIndex], 0, filledThisRound, value);
            fillLength -= filledThisRound;
            ++arrIndex;
        }
    }

    private static void backFillArrays(byte [][] array, byte value, int fillLength) {
        int arrIndex = array.length - 1;
        while(fillLength > 0) {
            final int filledThisRound = Math.min(array[arrIndex].length, fillLength);
            Arrays.fill(array[arrIndex], array[arrIndex].length - filledThisRound, array[arrIndex].length-1, value);
            fillLength -= filledThisRound;
            --arrIndex;
        }
    }

    private static Map<Integer,QseqReadData> splitReadIntoLengths(final int thisReadLength, final int [] outputLengths, final Map<Integer, QseqReadData> templates, final int writeOffset) {
        int totalLength = 0;
        for(int i = 0; i < outputLengths.length; i++) {
            totalLength += outputLengths[i];
        }

        final Map<Integer,QseqReadData> outData = new HashMap<Integer, QseqReadData>();
        for(final Map.Entry<Integer, QseqReadData> readToTemplate : templates.entrySet()) {
            final QseqReadData template = readToTemplate.getValue();
            final QseqReadData qseqRd = new QseqReadData(outputLengths);
            fillArrays(qseqRd.getBases(),     (byte)0, writeOffset);
            fillArrays(qseqRd.getQualities(), (byte)0, writeOffset);
            backFillArrays(qseqRd.getBases(),     (byte)0, totalLength - thisReadLength - writeOffset);
            backFillArrays(qseqRd.getQualities(), (byte)0, totalLength - thisReadLength - writeOffset);

            inMemorySplit(template.getBases(),     qseqRd.getBases(),     writeOffset);
            inMemorySplit(template.getQualities(), qseqRd.getQualities(), writeOffset);
            qseqRd.setOrCheckPf(template.isPf());
            qseqRd.setOrCheckXCoordinate(template.getXCoordinate());
            qseqRd.setOrCheckYCoordinate(template.getYCoordinate());
            qseqRd.setOrCheckLane(template.getLane());
            qseqRd.setOrCheckTile(template.getTile());
            outData.put(readToTemplate.getKey(), qseqRd);
        }

        return outData;
    }

    private static int totalLength(byte [][] arrs) {
        int total = 0;
        for(int i = 0; i < arrs.length; i++) {
            for(byte[] arr : arrs) {
                total += arr.length;
            }
        }

        return total;
    }

    private static void flattenWLimit(byte [][] src, byte [] dest, int srcStart) {
        int srcArrIndex = 0;
        while(srcStart >= src[srcArrIndex].length) {
            srcStart -= src[srcArrIndex++].length;
        }

        int written = 0;


        while(written < dest.length) {

            final int toWrite = Math.min(dest.length, src[srcArrIndex].length);
            System.arraycopy(src[srcArrIndex], srcStart, dest, 0, toWrite);
            written += toWrite;

            ++srcArrIndex; //if there is another iteration it is because src[srcArrIndex] was too short
            srcStart = 0;
        }
    }


    private static QseqReadData [] splitQSeq(final QseqReadData qseqReadData, int splitStart, int splitLength) {
        int availableLength = totalLength(qseqReadData.getBases());
        int totalQseqReads = 1;
        if(splitStart > 0) {
            ++totalQseqReads;
        }

        int lastStart = splitStart + splitLength;
        if(lastStart < availableLength) {
            ++totalQseqReads;
        }

        QseqReadData [] outData = new QseqReadData[totalQseqReads];

        QseqReadData qrd;
        int outIndex = 0;
        if(splitStart > 0) {
            qrd = new QseqReadData(new int[]{splitStart});
            flattenWLimit(qseqReadData.getBases(), qrd.getBases()[0], 0);
            flattenWLimit(qseqReadData.getQualities(), qrd.getQualities()[0], 0);
            outData[outIndex++] = qrd;
        }


        qrd = new QseqReadData(new int[]{splitLength});
        flattenWLimit(qseqReadData.getBases(), qrd.getBases()[0], splitStart);
        flattenWLimit(qseqReadData.getQualities(), qrd.getQualities()[0],splitStart);
        outData[outIndex++] = qrd;

        if(splitStart + splitLength < availableLength) {
            qrd = new QseqReadData(new int[]{availableLength - lastStart});
            flattenWLimit(qseqReadData.getBases(), qrd.getBases()[0], lastStart);
            flattenWLimit(qseqReadData.getQualities(), qrd.getQualities()[0], lastStart);
            outData[outIndex++] = qrd;
        }

        for(int i = 0; i < outData.length; i++) {
            outData[i].setOrCheckPf(qseqReadData.isPf());
            outData[i].setOrCheckXCoordinate(qseqReadData.getXCoordinate());
            outData[i].setOrCheckYCoordinate(qseqReadData.getYCoordinate());
            outData[i].setOrCheckLane(qseqReadData.getLane());
            outData[i].setOrCheckTile(qseqReadData.getTile());
        }

        return outData;
    }

    private static int flatten(byte [][] bytes, byte [] out, int index) {
        int outIndex = index;
        for(int i = 0; i < bytes.length; i++) {
            for(int j = 0; j < bytes[i].length; j++) {
                out[outIndex++] = bytes[i][j];
            }
        }
        return outIndex;
    }

    private static QseqReadData combineReads(final QseqReadData qrd1, final QseqReadData qrd2) {
        int totalSize = 0;
        byte [][] bases1 = qrd1.getBases();
        byte [][] qualities1 = qrd1.getQualities();

        byte [][] bases2 = qrd2.getBases();
        byte [][] qualities2 = qrd2.getQualities();

        for(byte [] bases : bases1) {
            totalSize += bases.length;
        }

        for(byte [] bases : bases2) {
            totalSize += bases.length;
        }

        final QseqReadData outQrd = new QseqReadData(new int[]{totalSize});
        byte [] outBases = outQrd.getBases()[0];
        byte [] outQuals = outQrd.getQualities()[0];
        int bases2Start = flatten(bases1, outBases, 0);
        flatten(bases2, outBases, bases2Start);
        int quals2Start = flatten(qualities1, outQuals, 0); //should be the same as bases2Start
        flatten(qualities2, outQuals, quals2Start);

        outQrd.setOrCheckLane(qrd1.getLane());
        outQrd.setOrCheckTile(qrd1.getTile());
        outQrd.setOrCheckPf(qrd1.isPf());
        outQrd.setOrCheckXCoordinate(qrd1.getXCoordinate());
        outQrd.setOrCheckYCoordinate(qrd1.getYCoordinate());
        return outQrd;
    }

    private static Map<Integer, QseqReadData> combineReads(Map<Integer, QseqReadData> map1, Map<Integer, QseqReadData> map2) {
        final Map<Integer, QseqReadData> outMap = new HashMap<Integer, QseqReadData>();
        if(map1.size() != map2.size()) {
            throw new PicardException("Map1 and Map2 are not of the same size!");
        }

        for(Integer key : map1.keySet()) {
            final QseqReadData read1 = map1.get(key);
            final QseqReadData read2 = map2.get(key);
            if(read1 == null) {
                throw new PicardException("Null value in map 1 for key " + key);
            }
            if(read2 == null) {
                throw new PicardException("Null value in map 2 for key " + key);
            }

            outMap.put(key, combineReads(read1, read2));
        }

        return outMap;
    }

    private static <T,V> Map<T,V> toMap(List<T> keys, List<V> values) {
        if(keys.size() != values.size())
            throw new PicardException("Key list is not the same size as the value list!");

        Map<T, V> hmap = new HashMap<T, V>();
        for(int i = 0; i < keys.size(); i++) {
            hmap.put(keys.get(i), values.get(i));
        }

        return hmap;
    }

    public static Map<Integer, QseqReadData> addToKeys(Map<Integer, QseqReadData> map, int addend) {
        Map<Integer, QseqReadData> outMap = new HashMap<Integer, QseqReadData>();
        for(final Map.Entry<Integer, QseqReadData> entry : map.entrySet()) {
            outMap.put(entry.getKey() + addend, entry.getValue());
        }
        return outMap;
    }
    
    public static Map<Integer, QseqReadData> getReadData(int lane, int end, int tile) {
        final Map<Integer, QseqReadData> qMap = new HashMap<Integer, QseqReadData>();
        return toMap(getReadNos(lane, 1), fnPrefixToQSeqReadData.get(makeFnPrefix(lane, end, tile)));
    }

    public static Map<Integer, QseqReadData> getReadData(final String fileNamePrefix) {
       int [] vars = varsFromPrefix(fileNamePrefix);
       return getReadData(vars[LANE], vars[READ], vars[TILE]);
    }

    public static Map<Integer, QseqReadData> getTiledReadData(final String fileNamePrefix_noTile, final List<Integer> tiles) {
        int [] vars = varsFromNoTilePrefix(fileNamePrefix_noTile);
        return getTiledReadData(vars[LANE], vars[READ], tiles);
    }

    public static Map<Integer, QseqReadData> getTiledReadData(int lane, int end, final List<Integer> tiles) {
        final Map<Integer, QseqReadData> outMap = new HashMap<Integer, QseqReadData>();
        int totalAddend = 0;
        for(final Integer tile : tiles) {
            outMap.putAll(addToKeys(getReadData(lane, end, tile), totalAddend));
            totalAddend += getLaneFileSize(lane);
        }
        return outMap;
    }

    //final int lane, final int [] writeLengths, final List<IlluminaFileMap> readNumberToTile, final Map<Integer, QseqReadData> testAgainst
    public static Map<Integer, QseqReadData> getTiledReadData(final List<String> fileNames) {
        final Map<Integer, QseqReadData> outMap = new HashMap<Integer, QseqReadData>();
        int totalAddend = 0;
        for(final String fn : fileNames) {
            int [] vars = varsFromPrefix(fn);
            outMap.putAll(addToKeys(getReadData(vars[LANE], vars[READ], vars[TILE]), totalAddend));
            totalAddend += getLaneFileSize(vars[LANE]);
        }
        return outMap;
    }

    public static Map<Integer, QseqReadData> getTiledReadData(final List<String> fileNames, int totalReadLength, int [] outputLengths, int offset) {
        return splitReadIntoLengths(totalReadLength, outputLengths, getTiledReadData(fileNames), offset);
    }


    //Note: we want to use this like getReadNos(1,1,1) if we had 3 tiles for lane 1
    public static List<Integer> getReadNos(int lane, int files) {
        final List<Integer> readNos = new ArrayList<Integer>();
        int addend = 0;
        int laneSize = getLaneFileSize(lane);
        final List<Integer> laneNos = laneToReadNos.get(lane);
        for(int i = 0; i < files; i++) {
            for(Integer rn : laneNos) {
                readNos.add(rn + addend);
            }

            addend += laneSize;
        }

        return readNos;
    }

    public static Map<Integer,QseqReadData> getSplitOffsetReadData(final String fileNamePrefix, int totalLength, int [] outputLengths, int offsets) {
        return splitReadIntoLengths(totalLength, outputLengths, getReadData(fileNamePrefix), offsets);
    }

    public static Map<Integer, QseqReadData> getSplitOffsetReadData(final String fileName_preTile, List<Integer> tiles, int totalLength, int [] outputLengths, int offset) {
        return splitReadIntoLengths(totalLength, outputLengths, getTiledReadData(fileName_preTile, tiles), offset);
    }

    public static Map<Integer, QseqReadData> combineReads(final Map<Integer, QseqReadData> ... ends) {
        Map<Integer, QseqReadData> outList = new HashMap<Integer, QseqReadData>();
        outList.putAll(ends[0]);
        for(int i = 1; i < ends.length; i++) {
            outList = combineReads(outList, ends[i]);
        }

        return outList;
    }

    public static Map<Integer, QseqReadData> combineReads(int totalLength, int [] outputLengths, int offsets, final Map<Integer, QseqReadData> ... ends) {
        return splitReadIntoLengths(totalLength, outputLengths, combineReads(ends), offsets);
    }

    public static List<File> getQseqs(final String ... prefixes) {
        final List<File> files = new ArrayList<File>();
        for(final String prefix : prefixes) {
            files.add(fnPrefixToQSeqFile.get(prefix));
        }

        return files;
    }

    public static Map<Integer, ClusterData> qseqDataToClusterMap(final Map<Integer, QseqReadData> read1, Map<Integer, QseqReadData> read2, final int barcodeCycle, final int barcodeLength, final IlluminaDataType ... dataTypes) {
        if(barcodeCycle < 1 || barcodeLength <= 0) {
            throw new PicardException("This testdata method is to be used only with valid barcode data (i.e. barcodeCycle > 0 and barcodeLength > 0");
        }

        final Map<Integer, QseqReadData> combinedReads;

        if(read2 != null) {
            combinedReads = combineReads(read1, read2);
        } else {
            combinedReads = read1;
        }

        final Map<Integer, ClusterData> readNoToClusterData = new HashMap<Integer, ClusterData>();
        for(final Integer key : combinedReads.keySet()) {
            final QseqReadData [] qrds = splitQSeq(combinedReads.get(key), barcodeCycle - 1, barcodeLength);

            QseqReadData curRead1 = null;
            QseqReadData curRead2 = null;
            QseqReadData curBarcode;
            if(barcodeCycle > 1) {
                curRead1 = qrds[0];
                curBarcode = qrds[1];
                if(qrds.length > 2) {
                    curRead2 = qrds[2];
                }
            } else {
                curBarcode = qrds[0];
                if(qrds.length > 1) {
                    curRead2 = qrds[1];
                }
            }

            readNoToClusterData.put(key, qseqDataToClusterData(
                (curRead1   == null) ? null : curRead1,
                (curRead2   == null) ? null : curRead2,
                curBarcode,
                dataTypes)
            );
        }

        return readNoToClusterData;
    }

    public static Map<Integer, ClusterData> qseqDataToClusterMap(final Map<Integer, QseqReadData> read1, Map<Integer, QseqReadData> read2, Map<Integer, QseqReadData> barcode, final IlluminaDataType ... dataTypes) {
        final Set<Integer> keySet;
        if(read1 != null) {
            keySet = read1.keySet();
        } else if(read2 != null) {
            keySet = read2.keySet();
        } else if(barcode != null) {
            keySet = barcode.keySet();
        } else {
            throw new PicardException("Read to QseqReadData maps are all null!");
        }

        final Map<Integer, ClusterData> readNoToClusterData = new HashMap<Integer, ClusterData>();
        for(final Integer key : keySet) {
            readNoToClusterData.put(key, qseqDataToClusterData(
                (read1   == null) ? null : read1.get(key),
                (read2   == null) ? null : read2.get(key),
                (barcode == null) ? null : barcode.get(key),
                dataTypes)
            );
        }

        return readNoToClusterData;
    }

    public static ClusterData qseqDataToClusterData(final QseqReadData read1, final QseqReadData read2, final QseqReadData barcode,  final IlluminaDataType ... dataTypes) {
        int arrSize = ((read1 != null)   ? 1 : 0) +
                      ((read2 != null)   ? 1 : 0) +
                      ((barcode != null) ? 1 : 0);
        int readIndex = 0;
        ReadData [] reads = new ReadData[arrSize];

        ReadData rd1 = makeRead(ReadType.Template, read1,   dataTypes);
        if(rd1 != null) {
            reads[readIndex++] = rd1;
        }
        ReadData rdBarcode = makeRead(ReadType.Barcode,  barcode, dataTypes);
        if(rdBarcode != null) {
            reads[readIndex++] = rdBarcode;
        }
        ReadData rd2 = makeRead(ReadType.Template, read2,   dataTypes);
        if(rd2 != null) {
            reads[readIndex++] = rd2;
        }

        final ClusterData cd = new ClusterData(reads);
        for(IlluminaDataType idt : dataTypes) {
            switch(idt) {
                case Position:
                    cd.setTile(read1.getTile());
                    cd.setLane(read1.getLane());
                    cd.setX(read1.getXCoordinate());
                    cd.setY(read1.getYCoordinate());
                    break;

                case PF:
                    cd.setPf(read1.isPf());
                    break;
            }
        }

        return cd;
    }

    public static ReadData makeRead(final ReadType rt, final QseqReadData read, final IlluminaDataType ... dataTypes) {
        if(read != null) {
            ReadData rd = new ReadData(rt);

            for(IlluminaDataType idt : dataTypes) {
                switch(idt) {
                    case BaseCalls:
                        rd.setBases(read.getBases()[0]);
                        break;

                    case QualityScores:
                        rd.setQualities(read.getQualities()[0]);
                        break;
                }
            }

            return rd;
        }

        return null;
    }
}
