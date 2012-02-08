package net.sf.picard.illumina.parser;

import static net.sf.picard.util.CollectionUtil.makeList;

import java.io.File;
import java.util.*;

//Illumina Dir Test Data
public class BinTdUtil {
    public static final File IntensitiesDir = new File("testdata/net/sf/picard/illumina/CompleteIlluminaDir/Intensities/");
    public static final File basecallDir = new File("testdata/net/sf/picard/illumina/CompleteIlluminaDir/Intensities/BaseCalls");

    public static final String ltStr(final int lane, final int tile) {
        return "s_" + lane + "_" + tile;
    }

    public static final byte A = (byte)65;
    public static final byte C = (byte)67;
    public static final byte G = (byte)71;
    public static final byte T = (byte)84;
    public static final byte P = (byte)46; //dot
    public static final Map<String, List<ClusterData>> goldData     = new HashMap<String, List<ClusterData>>();
    public static final Map<String, List<Integer>>     goldIndices  = new HashMap<String, List<Integer>>();
    public static final Map<String, Integer>           goldSizes    = new HashMap<String, Integer>();
    static {
        int lane = 1;
        int tile = 1101;
        goldIndices.put(ltStr(lane, tile), makeList(0, 2, 10, 18, 19));
        goldData.put(ltStr(lane, tile), makeList(
            makeCd(lane, tile, 1075, 1001, true,
                new byte[]{ P,  G,  T,  G,  T,  G,  A,  A,  G,  T,  T,  C,  A,  A,  G,  A,  A,  G,  C,  A,  G,  A,  C,  A,  A,  P,  G,  G,  T,  C,  G,  C,  A,  A,  A,  C,  T,  G,  A,  P,  P,  A,  A,  A,  A,  A,  C,  A,  C,  T,  C,  C,  A,  A,  A,  T,  A,  A },
                new byte[]{ 2, 16, 25, 35, 35, 35, 35, 30, 33, 37, 37, 39, 37, 35, 38, 37, 39, 37, 39, 40, 36, 39, 29, 31, 27,  2, 16, 28, 33, 35, 37, 37, 37, 27, 31, 31, 35, 30, 33,  2,  2, 17, 17, 34, 37, 37, 27, 34, 35, 30, 36, 39, 38, 10, 17, 27, 27, 34 },
                "AGGTCGCA"),
            makeCd(lane, tile, 1184, 1010, true,
                new byte[]{ P,  G,  T,  G,  A,  A,  C,  A,  A,  A,  C,  T,  G,  G,  T,  T,  T,  T,  A,  T,  C,  T,  G,  G,  T,  P,  G,  G,  T,  C,  G,  C,  A,  C,  C,  C,  T,  A,  T,  C,  A,  G,  T,  T,  T,  G,  T,  T,  G,  T,  C,  A,  C,  T,  T,  C,  C,  A },
                new byte[]{ 2, 16, 28, 35, 35, 35, 37, 37, 39, 39, 39, 39, 39, 41, 40, 41, 41, 41, 36, 39, 39, 40, 40, 41, 38,  2, 16, 28, 35, 35, 37, 37, 37, 33, 33, 34, 37, 32, 35, 35, 36, 39, 39, 39, 39, 39, 40, 41, 40, 41, 41, 41, 39, 39, 38, 39, 39, 41 },
                "AGGTCGCA"),
            makeCd(lane, tile, 1175, 1044, true,
                new byte[]{ C,  G,  G,  T,  C,  A,  G,  C,  A,  A,  A,  G,  G,  C,  T,  A,  T,  T,  C,  T,  C,  A,  T,  C,  T,  P,  G,  G,  T,  C,  G,  C,  A,  G,  T,  C,  C,  C,  A,  G,  G,  G,  A,  G,  G,  A,  G,  T,  C,  C,  C,  T,  G,  T,  G,  G,  G,  A },
                new byte[]{ 34, 34, 34, 37, 37, 37, 37, 37, 39, 39, 39, 39, 39, 40, 41, 41, 40, 41, 41, 40, 40, 34, 39, 37, 39,  2, 16, 28, 32, 35, 37, 37, 37, 31, 33, 34, 35, 37, 35, 35, 37, 39, 35, 37, 37, 30, 37, 37, 38, 41, 41, 41, 41, 40, 40, 40, 40, 31 },
                "AGGTCGCA"),
            makeCd(1, 1101, 1006, 1090, false,
                new byte[]{ P,  A,  C,  C,  P,  C,  T,  C,  A,  G,  G,  A,  G,  C,  A,  G,  A,  G,  T,  G,  P,  P,  P,  P,  P,  P,  P,  G,  T,  C,  G,  C,  P,  P,  P,  P,  P,  P,  P,  P,  P,  P,  P,  P,  P,  C,  A,  C,  C,  A,  P,  P,  P,  P,  P,  T,  C,  P },
                new byte[]{ 2,  7,  7, 17,  2, 17, 27, 12, 12, 27, 29, 30, 26, 30, 27, 30, 30, 30,  2,  2,  2,  2,  2,  2,  2,  2,  2, 16, 19, 32, 35, 37,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2 },
                null),
            makeCd(lane, tile, 1018, 1110, false,
                new byte[]{ P,  G,  A,  G,  G,  T,  T,  T,  T,  C,  T,  C,  C,  A,  G,  C,  A,  C,  C,  C,  A,  P,  P,  P,  C,  P,  G,  G,  T,  C,  G,  C,  A,  C,  P,  P,  P,  P,  G,  P,  P,  G,  C,  C,  A,  G,  A,  C,  C,  A,  G,  G,  T,  A,  G,  T,  T,  C },
                new byte[]{ 2, 15, 15, 28, 31, 27, 30, 30, 30, 30, 30, 30, 30, 31, 30, 31, 30, 30, 31, 30, 30,  2,  2,  2,  2,  2, 16, 28, 33, 35, 35, 37, 35, 33,  2,  2,  2,  2, 17,  2,  2, 17, 17, 32, 37, 37, 38, 40, 41, 41, 38, 41, 33, 34, 38, 37, 39, 40 },
                "AGGTCGCA")
        ));
        goldSizes.put(ltStr(lane, tile), 20);

        tile = 1201;
        goldIndices.put(ltStr(lane, tile), makeList(0, 1, 18, 19));
        goldData.put(ltStr(lane, tile), makeList(
            makeCd(lane, tile, 1134, 1002, true,
                new byte[]{ P,  A,  C,  C,  T,  G,  C,  G,  T,  G,  T,  C,  A,  G,  C,  A,  A,  C,  A,  T,  C,  C,  G,  C,  C,  P,  G,  G,  T,  C,  G,  C,  A,  A,  P,  G,  G,  G,  C,  P,  P,  C,  A,  T,  C,  A,  C,  A,  G,  G,  A,  G,  C,  C,  T,  G,  C,  C },
                    new byte[]{ 2, 16, 28, 32, 33, 35, 35, 35, 39, 39, 39, 37, 30, 37, 35, 27, 26, 32, 37, 38, 31, 36, 32, 39, 40,  2, 16, 28, 33, 35, 37, 37, 35, 31,  2, 16, 19, 32, 35,  2,  2, 18, 27, 34, 37, 37, 35, 39, 38, 39, 40, 40, 40, 38, 34, 32, 34, 38 },
                    "AGGTCGCA"),
            makeCd(lane, tile, 1224, 1012, false,
                new byte[]{ P,  A,  T,  A,  A,  T,  G,  A,  C,  C,  T,  G,  G,  G,  G,  C,  T,  A,  C,  T,  G,  A,  A,  T,  C,  P,  G,  G,  T,  C,  G,  C,  A,  G,  C,  C,  T,  T,  T,  A,  C,  T,  T,  A,  T,  G,  A,  T,  C,  A,  C,  A,  G,  T,  T,  T,  A,  T },
                    new byte[]{ 2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 19, 28, 32, 35, 35, 37, 37,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2 },
                    "AGGTCGCA"),
            makeCd(lane, tile, 1136, 1112, true,
                new byte[]{ A,  G,  G,  A,  A,  T,  G,  G,  C,  T,  T,  G,  C,  T,  T,  A,  A,  G,  A,  C,  T,  T,  G,  C,  C,  P,  G,  G,  T,  C,  G,  C,  A,  T,  A,  T,  C,  C,  A,  A,  A,  A,  T,  C,  A,  A,  A,  A,  T,  G,  A,  A,  A,  T,  G,  C,  A,  G },
                    new byte[]{ 34, 34, 34, 37, 37, 37, 37, 37, 38, 39, 39, 16, 32, 36, 37, 36, 39, 35, 34, 38, 40, 39, 34, 37, 39,  2, 16, 28, 33, 35, 35, 37, 37, 31, 34, 34, 37, 37, 37, 37, 37, 39, 39, 39, 39, 39, 40, 41, 40, 40, 40, 41, 40, 41, 41, 40, 36, 39 },
                    "AGGTCGCA"),
            makeCd(lane, tile, 1201, 1117, true,
                new byte[]{ C,  C,  T,  C,  A,  G,  T,  T,  T,  C,  G,  G,  G,  A,  G,  A,  T,  C,  A,  T,  C,  C,  A,  C,  A,  P,  G,  G,  T,  C,  G,  C,  A,  C,  T,  C,  A,  A,  T,  A,  C,  T,  A,  A,  A,  T,  C,  T,  C,  A,  A,  G,  A,  T,  T,  G,  A,  T },
                    new byte[]{ 27, 31, 31, 35, 33, 35, 35, 35, 39, 39, 39, 39, 35, 31, 30, 33, 35, 32, 34, 32, 39, 29, 27, 37, 24,  2, 16, 28, 33, 35, 35, 37, 37, 30, 31, 31, 37, 35, 35, 36, 29, 37, 35, 35, 39, 33, 31, 37, 38, 31, 25, 25, 10, 32, 40, 35, 33, 39 },
                    "AGGTCGCA")
        ));
        goldSizes.put(ltStr(lane, tile), 20);

        tile = 2101;
        goldIndices.put(ltStr(lane, tile), makeList(7, 15, 16, 19));
        goldData.put(ltStr(lane, tile), makeList(
            makeCd(lane, tile, 1218, 1037, true,
                new byte[]{ C,  C,  G,  G,  A,  C,  C,  C,  T,  G,  G,  G,  C,  T,  C,  A,  G,  C,  C,  C,  T,  G,  A,  T,  G,  P,  G,  G,  T,  C,  G,  C,  A,  T,  G,  C,  T,  G,  G,  A,  C,  T,  G,  T,  T,  T,  G,  T,  G,  C,  A,  G,  G,  C,  G,  G,  C,  A },
                    new byte[]{ 31, 34, 34, 37, 37, 35, 37, 37, 39, 38, 39, 39, 39, 41, 41, 41, 41, 41, 40, 41, 40, 41, 38, 40, 40,  2, 16, 28, 33, 35, 35, 37, 37, 34, 31, 31, 37, 37, 37, 37, 37, 38, 39, 39, 37, 37, 40, 38, 40, 41, 40, 40, 40, 41, 40, 41, 39, 33 },
                    "AGGTCGCA"),
            makeCd(lane, tile, 1069, 1101, false,
                new byte[]{ G,  P,  T,  T,  P,  P,  P,  T,  T,  A,  A,  A,  C,  A,  T,  G,  G,  T,  G,  C,  T,  T,  A,  G,  T,  P,  G,  G,  T,  C,  G,  C,  A,  A,  A,  T,  T,  C,  C,  C,  C,  C,  T,  C,  C,  C,  C,  C,  A,  T,  C,  T,  A,  T,  A,  A,  T,  T },
                    new byte[]{ 27,  2, 15, 28,  2,  2,  2, 17, 18, 25, 29, 30, 31, 30, 30, 31, 31, 31, 30, 31, 31, 30, 30, 31, 30,  2, 16, 28, 33, 35, 37, 37, 37, 31, 33, 31, 37, 35, 35, 36, 37, 37, 39, 39, 39, 39, 41, 41, 40, 40, 38, 40, 40, 38, 34, 33, 37, 40 },
                    "AGGTCGCA"),
            makeCd(lane, tile, 1209, 1106, true,
                new byte[]{ C,  A,  A,  A,  T,  C,  T,  G,  A,  T,  T,  G,  C,  A,  T,  T,  A,  T,  T,  C,  A,  C,  C,  T,  G,  P,  G,  G,  T,  C,  G,  C,  A,  A,  G,  C,  T,  T,  G,  T,  T,  T,  G,  G,  G,  G,  G,  T,  C,  T,  G,  C,  T,  G,  T,  A,  G,  T },
                    new byte[]{ 34, 34, 34, 37, 37, 35, 37, 37, 39, 39, 39, 39, 39, 40, 40, 41, 40, 40, 41, 38, 40, 41, 40, 41, 41,  2, 16, 28, 33, 35, 37, 37, 37, 33, 33, 34, 37, 37, 37, 35, 37, 39, 39, 39, 39, 38, 41, 37, 39, 40, 41, 41, 41, 41, 39, 40, 40, 40 },
                    "AGGTCGCA"),
            makeCd(lane, tile, 1103, 1112, false,
                new byte[]{ A,  C,  C,  A,  P,  P,  P,  T,  G,  C,  C,  C,  C,  G,  G,  T,  G,  G,  G,  G,  G,  C,  T,  G,  C,  P,  G,  G,  T,  C,  G,  C,  A,  A,  G,  G,  A,  G,  G,  G,  G,  G,  G,  G,  G,  G,  G,  T,  G,  A,  G,  G,  G,  A,  A,  C,  C,  G },
                    new byte[]{ 26, 27, 27, 31,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 19, 28, 33, 35, 37, 37, 37, 25, 30,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2 },
                    "AGGTCGCA")
        ));
        goldSizes.put(ltStr(lane, tile), 20);

        lane = 2;

        tile = 1101;
        goldIndices.put(ltStr(lane, tile), makeList(0, 8, 9, 10));
        goldData.put(ltStr(lane, tile), makeList(
            makeCd(lane, tile, 1178, 1003, false,
                new byte[]{ P,  C,  A,  A,  P,  P,  P,  C,  A,  C,  C,  T,  C,  C,  C,  A,  G,  A,  G,  T,  C,  G,  T,  G,  G,  A,  G,  G,  T,  C,  G,  C,  A,  A,  G,  P,  G,  A,  C,  C,  A,  G,  A,  A,  G,  C,  T,  C,  T,  C,  T,  G,  A,  G,  T,  T,  A,  G },
                    new byte[]{ 2, 15, 26, 31,  2,  2,  2, 17, 17, 28, 28, 29, 30, 30, 31, 30, 30, 29, 30, 30, 30, 30, 31,  2,  2, 33, 33, 34, 37, 37, 37, 37, 35, 31, 33,  2, 19, 32, 35, 35, 37, 39, 39, 38, 39, 39, 40, 40, 41, 38, 40, 41, 34, 39, 39, 40, 40, 41 },
                    null),
            makeCd(lane, tile, 1215, 1071, true,
                new byte[]{ A,  G,  G,  G,  A,  T,  T,  T,  C,  A,  G,  T,  T,  G,  G,  T,  G,  G,  T,  G,  G,  G,  G,  C,  A,  A,  G,  G,  T,  C,  G,  C,  A,  G,  A,  A,  C,  A,  T,  A,  T,  T,  T,  G,  T,  C,  T,  C,  A,  G,  A,  G,  A,  G,  A,  C,  T,  T },
                    new byte[]{ 30, 31, 31, 32, 35, 35, 35, 35, 32, 23, 35, 30, 34, 38, 36, 27, 32, 38, 27, 34, 37, 38, 36, 32, 37, 33, 31, 34, 37, 37, 37, 37, 37, 31, 31, 31, 35, 35, 35, 33, 25, 35, 30, 39, 37, 32, 34, 37, 37, 39, 40, 40, 29, 38, 28, 35, 31, 37 },
                    null),
            makeCd(lane, tile, 1189, 1093, false,
                new byte[]{ C,  C,  A,  C,  C,  G,  T,  G,  A,  G,  A,  A,  T,  G,  C,  G,  C,  C,  A,  T,  C,  T,  G,  C,  A,  A,  G,  G,  T,  C,  G,  C,  A,  A,  A,  A,  A,  C,  A,  G,  A,  A,  C,  C,  C,  A,  G,  T,  T,  A,  G,  C,  A,  G,  T,  A,  T,  G },
                    new byte[]{ 28, 16, 22, 27, 16, 23, 32, 32, 10, 17, 28, 17, 17, 30, 19, 32, 21, 30,  8, 27, 10, 31, 27, 28, 31, 33, 33, 33, 37, 37, 37, 37, 37, 16, 23, 10, 32, 30, 10, 32, 32, 18, 10, 27, 34, 16, 10, 17, 27, 11, 11, 32, 34, 39,  2,  2,  2,  2 },
                    null),
            makeCd(lane, tile, 1171, 1103, false,
                new byte[]{ G,  T,  P,  G,  P,  P,  P,  G,  G,  A,  G,  A,  T,  G,  T,  A,  C,  T,  G,  G,  T,  G,  A,  A,  A,  A,  G,  G,  T,  C,  G,  C,  A,  T,  C,  T,  T,  T,  G,  G,  C,  T,  T,  G,  C,  T,  A,  A,  A,  T,  T,  T,  T,  A,  T,  T,  T,  A },
                    new byte[]{ 2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 31, 31, 34, 35, 35, 37, 37, 37, 16, 10, 25, 32, 10, 17, 10, 10, 32, 33, 31, 30, 30, 25, 30,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2 },
                    null)
        ));
        goldSizes.put(ltStr(lane, tile), 20);

        tile = 1201;
        goldIndices.put(ltStr(lane, tile), makeList(1, 3, 4, 17));
        goldData.put(ltStr(lane, tile), makeList(
            makeCd(lane, tile, 1228, 1012, true,
                new byte[]{ P,  T,  G,  G,  T,  C,  A,  T,  C,  T,  G,  C,  A,  G,  G,  T,  T,  T,  C,  T,  G,  A,  G,  A,  T,  P,  G,  G,  T,  C,  G,  C,  A,  T,  T,  T,  C,  T,  T,  A,  T,  A,  G,  A,  A,  A,  C,  A,  T,  C,  T,  T,  T,  T,  A,  T,  T,  T },
                    new byte[]{ 2, 16, 25, 33, 35, 33, 19, 35, 35, 35, 37, 35, 37, 38, 37, 35, 35, 38, 36, 37, 40, 34, 39, 33, 37,  2, 16, 25, 28, 28, 33, 31, 23, 31, 33, 31, 35, 37, 37, 35, 35, 33, 37, 37, 39, 39, 41, 40, 41, 37, 39, 41, 40, 40, 37, 39, 39, 34 },
                    null),
            makeCd(lane, tile, 1229, 1034, true,
                new byte[]{ C,  C,  T,  G,  T,  T,  G,  T,  A,  C,  G,  T,  C,  C,  C,  A,  G,  T,  A,  T,  G,  G,  A,  G,  C,  P,  G,  G,  T,  C,  G,  C,  A,  C,  T,  C,  T,  A,  T,  C,  C,  C,  A,  G,  A,  C,  C,  C,  T,  T,  C,  T,  C,  A,  G,  G,  C,  A },
                    new byte[]{ 34, 34, 34, 37, 37, 37, 37, 37, 38, 39, 39, 39, 39, 40, 38, 39, 39, 41, 40, 40, 41, 41, 24, 31, 34,  2, 16, 28, 32, 35, 35, 37, 37, 31, 31, 31, 35, 35, 37, 35, 35, 27, 34, 35, 39, 35, 39, 35, 37, 39, 38, 40, 36, 39, 38, 39, 39, 40 },
                    null),
            makeCd(lane, tile, 1194, 1053, true,
                new byte[]{ A,  G,  A,  T,  C,  T,  C,  A,  T,  A,  T,  C,  G,  T,  C,  G,  C,  T,  C,  G,  T,  C,  A,  T,  G,  P,  A,  T,  T,  A,  G,  A,  T,  T,  G,  T,  C,  G,  A,  T,  T,  A,  T,  C,  G,  C,  A,  C,  T,  G,  G,  T,  G,  C,  G,  A,  A,  T },
                    new byte[]{ 31, 31, 31, 37, 37, 35, 35, 37, 35, 37, 39, 39, 39, 38, 40, 38, 40, 40, 41, 41, 41, 36, 34, 36, 36,  2, 16, 16, 33, 28, 35, 32, 33, 31, 31, 31, 37, 37, 35, 35, 37, 37, 39, 39, 37, 30, 39, 33, 36, 36, 39, 36, 39, 39, 38, 30, 30, 30 },
                    null),
            makeCd(lane, tile, 1173, 1158, true,
                new byte[]{ A,  C,  G,  C,  P,  G,  G,  C,  A,  A,  T,  G,  A,  T,  G,  G,  C,  G,  T,  C,  T,  C,  G,  C,  A,  P,  G,  G,  T,  C,  G,  C,  A,  T,  G,  C,  A,  T,  C,  G,  A,  C,  T,  G,  C,  C,  C,  T,  A,  T,  T,  T,  T,  G,  T,  C,  T,  A },
                    new byte[]{ 31, 30, 31, 35,  2, 17, 32, 32, 37, 37, 37, 37, 37, 40, 40, 40, 40, 40, 40, 40, 40, 40, 38, 38, 34,  2, 16, 28, 35, 35, 37, 37, 37, 33, 31, 31, 30, 35, 35, 35, 35, 30, 37, 39, 37, 37, 38, 37, 38, 39, 38, 34, 29, 36, 38, 39, 41, 38 },
                    null)

        ));
        goldSizes.put(ltStr(lane, tile), 20);

        tile = 2101;
        goldIndices.put(ltStr(lane, tile), makeList(0, 11, 13, 19));
        goldData.put(ltStr(lane, tile), makeList(
            makeCd(lane, tile, 1139, 1012, true,
                new byte[]{ P,  G,  T,  A,  T,  C,  A,  T,  T,  T,  T,  G,  T,  T,  A,  A,  A,  A,  G,  T,  C,  C,  T,  A,  C,  P,  G,  G,  T,  C,  G,  C,  A,  C,  A,  A,  A,  A,  T,  P,  A,  T,  G,  C,  T,  C,  A,  A,  G,  A,  G,  T,  G,  C,  C,  C,  A,  G },
                    new byte[]{ 2, 16, 16, 32, 33, 35, 26, 35, 37, 37, 35, 30, 34, 38, 37, 37, 37, 33, 37, 36, 39, 40, 37, 35, 31,  2, 16, 28, 35, 35, 37, 37, 37, 27, 30, 26, 35, 35, 35,  2, 17, 17, 32, 27, 30, 34, 38, 39, 33, 38, 38, 40, 39, 34, 29, 31, 37, 36 },
                    null),
            makeCd(lane, tile, 1177, 1057, true,
                new byte[]{ A,  A,  G,  A,  G,  A,  A,  C,  A,  C,  G,  T,  T,  A,  T,  A,  G,  G,  A,  C,  A,  T,  T,  T,  T,  P,  G,  G,  T,  C,  G,  C,  A,  A,  G,  G,  C,  C,  A,  A,  G,  T,  G,  C,  C,  C,  C,  T,  C,  T,  T,  G,  G,  T,  G,  T,  C,  T },
                    new byte[]{ 31, 27, 30, 35, 35, 30, 35, 35, 35, 37, 39, 37, 37, 38, 39, 36, 39, 39, 34, 36, 38, 38, 39, 40, 40,  2, 16, 28, 35, 35, 37, 37, 37, 31, 31, 31, 35, 35, 32, 37, 35, 39, 39, 32, 37, 35, 30, 36, 38, 36, 38, 39, 41, 36, 38, 36, 24, 34 },
                    null),
            makeCd(lane, tile, 1070, 1062, false,
                new byte[]{ A,  P,  T,  A,  P,  P,  P,  C,  A,  G,  G,  G,  A,  T,  A,  A,  C,  T,  C,  A,  A,  P,  P,  G,  C,  P,  G,  G,  T,  C,  G,  C,  A,  A,  A,  C,  T,  C,  A,  G,  A,  C,  T,  A,  C,  A,  G,  G,  T,  G,  C,  A,  G,  A,  A,  A,  A,  C },
                    new byte[]{ 27,  2, 15, 31,  2,  2,  2, 19, 17, 31, 30, 31, 30, 30, 30, 30, 30, 30, 30, 30, 30,  2,  2,  2,  2,  2, 19, 28, 35, 35, 37, 37, 35, 31, 31, 34, 35, 35, 33, 37, 37, 37, 39, 39, 39, 39, 41, 40, 37, 39, 38, 38, 38, 39, 36, 39, 38, 38 },
                    null),
            makeCd(lane, tile, 1246, 1085, true,
                new byte[]{ T,  C,  T,  T,  C,  C,  C,  A,  T,  G,  A,  G,  G,  G,  C,  A,  C,  A,  G,  T,  T,  T,  G,  A,  C,  P,  G,  G,  T,  C,  G,  C,  A,  A,  G,  T,  C,  G,  A,  A,  T,  T,  G,  T,  A,  A,  T,  T,  C,  C,  A,  T,  T,  T,  G,  C,  C,  C },
                    new byte[]{ 34, 34, 34, 37, 37, 37, 37, 37, 39, 39, 39, 39, 39, 40, 39, 39, 40, 40, 39, 38, 39, 40, 41, 39, 40,  2, 16, 28, 33, 35, 35, 37, 37, 34, 34, 34, 37, 37, 37, 37, 37, 39, 39, 39, 39, 39, 40, 41, 41, 39, 40, 41, 41, 41, 41, 40, 41, 40 },
                    null)
        ));
        goldSizes.put(ltStr(lane, tile), 20);

        lane = 3;
        tile = 1101;
        goldIndices.put(ltStr(lane, tile), makeList(0, 8, 9, 10));
        goldData.put(ltStr(lane, tile), makeList(
            makeCd(3, 1101, 1224, 1003, true,
                    new byte[]{ P,  C,  C,  A,  A,  T,  T,  A,  C,  A,  T,  T,  T,  T,  C,  T,  C,  T,  C,  T,  A,  T,  C,  A,  C,  A,  G,  G,  T,  C,  G,  C,  A,  T,  T,  T,  T,  C,  T,  T,  G,  G,  T,  T,  G,  G,  G,  T,  C,  T,  C,  G,  G,  G,  C,  T,  T,  G },
                    new byte[]{ 2, 16, 25, 33, 35, 35, 37, 37, 39, 39, 39, 39, 39, 40, 38, 40, 39, 40, 40, 40, 39, 39, 40, 40, 40, 31, 33, 34, 35, 37, 37, 37, 37, 34, 34, 34, 37, 37, 37, 37, 37, 39, 39, 39, 39, 39, 40, 39, 40, 41, 41, 41, 41, 40, 38, 38, 36, 39 },
                    null),
            makeCd(3, 1101, 1195, 1128, false,
                    new byte[]{ T,  C,  P,  T,  P,  P,  C,  T,  C,  T,  G,  C,  G,  C,  T,  C,  C,  C,  T,  C,  T,  C,  G,  C,  T,  A,  G,  G,  T,  C,  G,  C,  A,  A,  C,  C,  A,  G,  G,  A,  G,  G,  A,  A,  A,  G,  T,  G,  T,  G,  G,  G,  G,  C,  C,  A,  G,  A },
                    new byte[]{ 27, 26,  2, 17,  2,  2, 17, 17, 31, 30, 24, 30, 28, 30, 30, 30, 27, 30, 30, 27, 31, 22, 28, 30, 30, 34, 34, 31, 37, 37, 37, 37, 37, 31, 30, 33, 33, 35, 35, 30, 35, 37, 34, 35, 37, 39, 40, 38, 33, 38, 38, 36, 38, 32, 30, 36, 34, 37 },
                    null),
            makeCd(3, 1101, 1232, 1155, true,
                    new byte[]{ A,  A,  A,  C,  C,  C,  A,  A,  T,  G,  A,  T,  C,  A,  G,  G,  T,  A,  T,  G,  T,  A,  C,  C,  C,  A,  G,  G,  T,  C,  G,  C,  A,  C,  T,  T,  A,  C,  C,  T,  A,  G,  T,  T,  T,  C,  C,  A,  G,  C,  A,  G,  T,  C,  T,  G,  C,  T },
                    new byte[]{ 34, 34, 34, 37, 37, 37, 37, 37, 39, 39, 33, 39, 37, 40, 40, 41, 36, 39, 38, 39, 39, 40, 41, 40, 40, 34, 34, 34, 35, 37, 37, 37, 37, 34, 34, 34, 37, 37, 37, 37, 35, 39, 39, 37, 39, 39, 41, 41, 40, 41, 40, 41, 41, 40, 41, 41, 34, 39 },
                    null),
            makeCd(3, 1101, 1240, 1185, false,
                    new byte[]{ C,  C,  A,  G,  C,  A,  G,  G,  C,  G,  G,  G,  G,  G,  C,  A,  G,  G,  G,  G,  G,  G,  C,  A,  G,  A,  G,  G,  T,  C,  G,  C,  A,  A,  C,  A,  G,  T,  A,  G,  G,  C,  A,  C,  T,  C,  A,  C,  T,  A,  C,  A,  T,  G,  C,  G,  G,  C },
                    new byte[]{ 2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 31, 33, 34, 35, 37, 37, 37, 37, 33, 10, 10, 32, 33, 25, 35, 35, 30, 35, 26, 17, 32, 32, 36, 33, 30, 10, 27, 34, 19, 10, 27, 31, 37 },
                    null)
        ));
        goldSizes.put(ltStr(lane, tile), 20);

        tile = 1201;
        goldIndices.put(ltStr(lane, tile), makeList(2, 16, 17, 19));
        goldData.put(ltStr(lane, tile), makeList(
            makeCd(3, 1201, 1208, 1031, true,
                    new byte[]{ C,  T,  G,  C,  P,  C,  A,  G,  G,  C,  A,  T,  T,  C,  A,  C,  A,  A,  T,  G,  G,  A,  G,  G,  C,  P,  G,  G,  T,  C,  G,  C,  A,  T,  G,  C,  T,  G,  G,  G,  A,  T,  T,  A,  C,  A,  G,  G,  C,  G,  T,  G,  G,  G,  G,  A,  C,  T },
                    new byte[]{ 31, 31, 31, 35,  2, 17, 28, 32, 37, 33, 30, 37, 37, 27, 34, 36, 27, 25, 34, 39, 36, 32, 37, 38, 27,  2, 16, 28, 32, 35, 35, 37, 35, 33, 31, 30, 35, 37, 32, 37, 30, 37, 39, 30, 35, 33, 36, 38, 33, 37, 39, 30,  2,  2,  2,  2,  2,  2 },
                    null),
            makeCd(3, 1201, 1242, 1201, true,
                    new byte[]{ A,  T,  C,  T,  G,  C,  C,  G,  C,  A,  C,  C,  T,  C,  T,  G,  A,  C,  T,  T,  T,  G,  T,  A,  C,  P,  G,  G,  T,  C,  G,  C,  A,  G,  G,  A,  A,  A,  G,  G,  G,  T,  A,  A,  T,  A,  A,  G,  T,  A,  T,  G,  A,  C,  T,  T,  G,  A },
                    new byte[]{ 31, 31, 27, 35, 35, 35, 35, 35, 32, 32, 37, 37, 37, 33, 18, 30, 27, 36, 37, 34, 39, 27, 34, 27, 24,  2, 19, 28, 33, 35, 37, 37, 35, 26, 30, 30, 35, 32, 35, 23, 35, 10, 32, 33, 37, 35, 39, 40, 19, 34, 37, 41, 39, 40, 38, 39, 39, 39 },
                    null),
            makeCd(3, 1201, 1220, 1217, true,
                    new byte[]{ A,  G,  A,  C,  C,  G,  G,  C,  G,  G,  A,  G,  A,  T,  G,  T,  G,  A,  A,  C,  G,  T,  G,  G,  G,  P,  G,  G,  T,  C,  G,  C,  A,  T,  G,  G,  G,  A,  A,  G,  T,  G,  C,  T,  G,  C,  A,  G,  G,  C,  T,  C,  A,  C,  T,  G,  C,  G },
                    new byte[]{ 30, 30, 30, 35, 35, 35, 35, 35, 35, 35, 25, 35, 30, 26, 33, 30, 35, 19, 33, 23, 35, 14, 33, 34,  2,  2, 16, 28, 33, 35, 35, 37, 37, 31, 34, 34, 37, 37, 37, 37, 35, 37, 39, 38, 35, 39, 40, 39, 40, 38, 40, 41, 38, 38, 38, 39, 40, 40 },
                    null),
            makeCd(3, 1201, 1364, 1001, false,
                    new byte[]{ P,  T,  T,  T,  A,  A,  A,  C,  C,  C,  C,  A,  A,  C,  A,  A,  T,  T,  A,  A,  T,  T,  T,  T,  A,  P,  G,  G,  T,  C,  G,  C,  A,  T,  G,  P,  P,  G,  A,  A,  G,  T,  G,  C,  A,  G,  T,  T,  A,  G,  A,  T,  C,  C,  T,  T,  C,  A },
                    new byte[]{ 2, 10, 16, 22, 32, 32, 32, 31, 32, 10, 30, 28, 32, 34, 27, 32, 32, 32, 34, 22, 22, 31, 33, 30, 30,  2, 16, 28, 33, 35, 37, 37, 37, 28, 33,  2,  2, 19, 10, 28, 32, 17, 32, 10, 17, 17, 27, 19, 34, 37, 10, 18, 32, 32, 34, 19, 10, 32 },
                    null)
        ));
        goldSizes.put(ltStr(lane, tile), 20);

        tile = 2101;
        goldIndices.put(ltStr(lane, tile), makeList(0, 1, 2, 3));
        goldData.put(ltStr(lane, tile), makeList(
            makeCd(3, 2101, 1125, 1003, true,
                new byte[]{ P,  G,  T,  C,  T,  G,  G,  C,  C,  T,  T,  G,  C,  T,  G,  C,  T,  G,  T,  C,  C,  C,  T,  G,  G,  A,  G,  G,  T,  C,  G,  C,  A,  C,  C,  C,  A,  G,  A,  T,  T,  T,  T,  T,  T,  A,  T,  C,  A,  C,  A,  G,  A,  G,  T,  T,  G,  A },
                new byte[]{ 2, 16, 28, 35, 35, 35, 37, 37, 39, 39, 39, 39, 39, 41, 41, 38, 40, 40, 40, 40, 41, 40, 40, 41, 41, 31, 34, 34, 35, 37, 37, 37, 37, 31, 31, 31, 37, 37, 35, 36, 37, 39, 39, 39, 39, 37, 40, 39, 38, 40, 38, 39, 38, 40, 39, 39, 40, 39 },
                null),
            makeCd(3, 2101, 1099, 1006, true,
                new byte[]{ P,  T,  G,  C,  T,  G,  T,  T,  T,  C,  A,  C,  A,  G,  C,  C,  T,  C,  T,  T,  C,  A,  A,  T,  C,  A,  G,  G,  T,  C,  G,  C,  A,  A,  T,  T,  T,  T,  A,  C,  T,  A,  T,  T,  C,  T,  T,  T,  C,  A,  G,  G,  C,  T,  T,  T,  C,  A },
                new byte[]{ 2, 16, 28, 35, 35, 37, 35, 35, 39, 37, 39, 39, 38, 41, 40, 40, 40, 40, 40, 41, 40, 39, 40, 40, 40, 34, 34, 34, 37, 37, 37, 37, 37, 31, 31, 31, 37, 37, 36, 35, 37, 37, 39, 39, 39, 39, 41, 41, 41, 41, 41, 40, 40, 41, 40, 41, 39, 37 },
                null),
            makeCd(3, 2101, 1198, 1014, true,
                new byte[]{ P,  C,  C,  A,  G,  A,  T,  G,  T,  T,  A,  C,  A,  T,  G,  G,  T,  G,  A,  G,  C,  C,  A,  G,  A,  A,  G,  G,  T,  C,  G,  C,  A,  T,  C,  A,  G,  T,  C,  T,  C,  T,  C,  A,  C,  T,  G,  T,  G,  C,  T,  G,  T,  G,  T,  C,  C,  T },
                new byte[]{ 2, 16, 28, 35, 35, 35, 37, 37, 39, 39, 39, 39, 39, 41, 41, 41, 38, 40, 39, 40, 40, 41, 40, 39, 39, 34, 34, 34, 37, 37, 37, 37, 37, 34, 34, 34, 37, 37, 37, 37, 37, 39, 37, 39, 39, 39, 41, 41, 41, 41, 41, 41, 40, 41, 41, 41, 41, 41 },
                null),
            makeCd(3, 2101, 1078, 1021, true,
                new byte[]{ G,  C,  A,  A,  A,  P,  G,  A,  T,  G,  A,  G,  A,  T,  C,  G,  C,  T,  G,  G,  G,  C,  C,  T,  G,  A,  G,  G,  T,  C,  G,  C,  A,  G,  C,  A,  C,  T,  C,  G,  G,  C,  C,  T,  C,  G,  A,  G,  T,  A,  T,  C,  C,  T,  T,  T,  A,  G },
                new byte[]{ 34, 34, 34, 37, 37,  2, 17, 32, 37, 39, 39, 39, 39, 38, 37, 40, 40, 40, 40, 40, 40, 39, 40, 40, 40, 31, 33, 34, 35, 37, 37, 37, 37, 34, 34, 34, 37, 37, 37, 37, 37, 39, 37, 39, 38, 39, 40, 41, 41, 40, 38, 40, 41, 39, 40, 40, 29, 37 },
                null)
        ));
        goldSizes.put(ltStr(lane, tile), 20);

        //Using the same data as Lane 2 except breaking it up differently, this will cause the EAMSS filtering to mess with the data a little, rather than 25T8B25T it will become 19T8B8B19T
        lane = 4;
        tile = 1101;
        goldIndices.put(ltStr(lane, tile), makeList(0, 8, 9, 10));
        goldData.put(ltStr(lane, tile), makeList(
                makeCd(lane, tile, 1178, 1003, false,
                        new byte[]{ P,  C,  A,  A,  P,  P,  P,  C,  A,  C,  C,  T,  C,  C,  C,  A,  G,  A,  G,  T,  C,  G,  T,  G,  G,  A,  G,  G,  T,  C,  G,  C,  A,  A,  G,  P,  G,  A,  C,  C,  A,  G,  A,  A,  G,  C,  T,  C,  T,  C,  T,  G,  A,  G,  T,  T,  A,  G },
                        new byte[]{ 2, 15, 26, 31,  2,  2,  2, 17, 17, 28, 28, 29, 30, 30, 31, 30, 30, 29, 30, 30, 30, 30, 31,  7,  7, 33, 33, 34, 37, 37, 37, 37, 35, 31, 33,  2, 19, 32, 35, 35, 37, 39, 39, 38, 39, 39, 40, 40, 41, 38, 40, 41, 34, 39, 39, 40, 40, 41 },
                        null),
                makeCd(lane, tile, 1215, 1071, true,
                        new byte[]{ A,  G,  G,  G,  A,  T,  T,  T,  C,  A,  G,  T,  T,  G,  G,  T,  G,  G,  T,  G,  G,  G,  G,  C,  A,  A,  G,  G,  T,  C,  G,  C,  A,  G,  A,  A,  C,  A,  T,  A,  T,  T,  T,  G,  T,  C,  T,  C,  A,  G,  A,  G,  A,  G,  A,  C,  T,  T },
                        new byte[]{ 30, 31, 31, 32, 35, 35, 35, 35, 32, 23, 35, 30, 34, 38, 36, 27, 32, 38, 27, 34, 37, 38, 36, 32, 37, 33, 31, 34, 37, 37, 37, 37, 37, 31, 31, 31, 35, 35, 35, 33, 25, 35, 30, 39, 37, 32, 34, 37, 37, 39, 40, 40, 29, 38, 28, 35, 31, 37 },
                        null),
                makeCd(lane, tile, 1189, 1093, false,
                        new byte[]{ C,  C,  A,  C,  C,  G,  T,  G,  A,  G,  A,  A,  T,  G,  C,  G,  C,  C,   A,  T,  C,  T,  G,  C,  A,  A,  G,  G,  T,  C,  G,  C,  A,  A,  A,  A,  A,  C,  A,  G,  A,  A,  C,  C,  C,  A,  G,  T,  T,  A,  G,  C,  A,  G,  T,  A,  T,  G },
                        new byte[]{ 28, 16, 22, 27, 16, 23, 32, 32, 10, 17, 28, 17, 17, 30, 19, 32, 21, 30,  2, 27, 10, 31, 27, 28, 31, 33, 33, 33, 37, 37, 37, 37, 37, 16, 23, 10, 32, 30,  2, 32, 32, 18, 10, 27, 34, 16, 10, 17, 27, 11, 11, 32, 34, 39,  2,  2,  2,  2 },
                        null),
                makeCd(lane, tile, 1171, 1103, false,
                        new byte[]{ G,  T,  P,  G,  P,  P,  P,  G,  G,  A,  G,  A,  T,  G,  T,  A,  C,  T,  G,  G,  T,  G,  A,  A,  A,  A,  G,  G,  T,  C,  G,  C,  A,  T,  C,  T,  T,  T,  G,  G,  C,  T,  T,  G,  C,  T,  A,  A,  A,  T,  T,  T,  T,  A,  T,  T,  T,  A },
                        new byte[]{ 2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 13, 27,  8,  8, 13, 21, 31, 31, 34, 35, 35, 37, 37, 37, 16, 10, 25, 32,  2,  2, 10, 10, 32, 33, 31, 30, 30, 25, 30,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2 },
                        null)
        ));
        goldSizes.put(ltStr(lane, tile), 20);

        tile = 1201;
        goldIndices.put(ltStr(lane, tile), makeList(1, 3, 4, 17));
        goldData.put(ltStr(lane, tile), makeList(
                makeCd(lane, tile, 1228, 1012, true,
                        new byte[]{ P,  T,  G,  G,  T,  C,  A,  T,  C,  T,  G,  C,  A,  G,  G,  T,  T,  T,  C,  T,  G,  A,  G,  A,  T,  P,  G,  G,  T,  C,  G,  C,  A,  T,  T,  T,  C,  T,  T,  A,  T,  A,  G,  A,  A,  A,  C,  A,  T,  C,  T,  T,  T,  T,  A,  T,  T,  T },
                        new byte[]{ 2, 16, 25, 33, 35, 33, 19, 35, 35, 35, 37, 35, 37, 38, 37, 35, 35, 38, 36, 37, 40, 34, 39, 33, 37,  2,  2,  2,  2, 28, 33, 31, 23, 31, 33, 31, 35, 37, 37, 35, 35, 33, 37, 37, 39, 39, 41, 40, 41, 37, 39, 41, 40, 40, 37, 39, 39, 34 },
                        null),
                makeCd(lane, tile, 1229, 1034, true,
                        new byte[]{ C,  C,  T,  G,  T,  T,  G,  T,  A,  C,  G,  T,  C,  C,  C,  A,  G,  T,  A,  T,  G,  G,  A,  G,  C,  P,  G,  G,  T,  C,  G,  C,  A,  C,  T,  C,  T,  A,  T,  C,  C,  C,  A,  G,  A,  C,  C,  C,  T,  T,  C,  T,  C,  A,  G,  G,  C,  A },
                        new byte[]{ 34, 34, 34, 37, 37, 37, 37, 37, 38, 39, 39, 39, 39, 40, 38, 39, 39, 41, 40, 40, 41, 41, 24, 31, 34,  2, 16, 28, 32, 35, 35, 37, 37, 31, 31, 31, 35, 35, 37, 35, 35, 27, 34, 35, 39, 35, 39, 35, 37, 39, 38, 40, 36, 39, 38, 39, 39, 40 },
                        null),
                makeCd(lane, tile, 1194, 1053, true,
                        new byte[]{ A,  G,  A,  T,  C,  T,  C,  A,  T,  A,  T,  C,  G,  T,  C,  G,  C,  T,  C,  G,  T,  C,  A,  T,  G,  P,  A,  T,  T,  A,  G,  A,  T,  T,  G,  T,  C,  G,  A,  T,  T,  A,  T,  C,  G,  C,  A,  C,  T,  G,  G,  T,  G,  C,  G,  A,  A,  T },
                        new byte[]{ 31, 31, 31, 37, 37, 35, 35, 37, 35, 37, 39, 39, 39, 38, 40, 38, 40, 40, 41, 41, 41, 36, 34, 36, 36,  2, 16, 16, 33, 28, 35, 32, 33, 31, 31, 31, 37, 37, 35, 35, 37, 37, 39, 39, 37, 30, 39, 33, 36, 36, 39, 36, 39, 39, 38, 30, 30, 30 },
                        null),
                makeCd(lane, tile, 1173, 1158, true,
                        new byte[]{ A,  C,  G,  C,  P,  G,  G,  C,  A,  A,  T,  G,  A,  T,  G,  G,  C,  G,  T,  C,  T,  C,  G,  C,  A,  P,  G,  G,  T,  C,  G,  C,  A,  T,  G,  C,  A,  T,  C,  G,  A,  C,  T,  G,  C,  C,  C,  T,  A,  T,  T,  T,  T,  G,  T,  C,  T,  A },
                        new byte[]{ 31, 30, 31, 35,  2, 17, 32, 32, 37, 37, 37, 37, 37, 40, 40, 40, 40, 40, 40, 40, 40, 40, 38, 38, 34,  2, 16, 28, 35, 35, 37, 37, 37, 33, 31, 31, 30, 35, 35, 35, 35, 30, 37, 39, 37, 37, 38, 37, 38, 39, 38, 34, 29, 36, 38, 39, 41, 38 },
                        null)

        ));
        goldSizes.put(ltStr(lane, tile), 20);

        tile = 2101;
        goldIndices.put(ltStr(lane, tile), makeList(0, 11, 13, 19));
        goldData.put(ltStr(lane, tile), makeList(
                makeCd(lane, tile, 1139, 1012, true,
                        new byte[]{ P,  G,  T,  A,  T,  C,  A,  T,  T,  T,  T,  G,  T,  T,  A,  A,  A,  A,  G,  T,  C,  C,  T,  A,  C,  P,  G,  G,  T,  C,  G,  C,  A,  C,  A,  A,  A,  A,  T,  P,  A,  T,  G,  C,  T,  C,  A,  A,  G,  A,  G,  T,  G,  C,  C,  C,  A,  G },
                        new byte[]{ 2, 16, 16, 32, 33, 35, 26, 35, 37, 37, 35, 30, 34, 38, 37, 37, 37, 33, 37, 36, 39, 40, 37, 35, 31,  2, 16, 28, 35, 35, 37, 37, 37, 27, 30, 26, 35, 35, 35,  2, 17, 17, 32, 27, 30, 34, 38, 39, 33, 38, 38, 40, 39, 34, 29, 31, 37, 36 },
                        null),
                makeCd(lane, tile, 1177, 1057, true,
                        new byte[]{ A,  A,  G,  A,  G,  A,  A,  C,  A,  C,  G,  T,  T,  A,  T,  A,  G,  G,  A,  C,  A,  T,  T,  T,  T,  P,  G,  G,  T,  C,  G,  C,  A,  A,  G,  G,  C,  C,  A,  A,  G,  T,  G,  C,  C,  C,  C,  T,  C,  T,  T,  G,  G,  T,  G,  T,  C,  T },
                        new byte[]{ 31, 27, 30, 35, 35, 30, 35, 35, 35, 37, 39, 37, 37, 38, 39, 36, 39, 39, 34, 36, 38, 38, 39, 40, 40,  2, 16, 28, 35, 35, 37, 37, 37, 31, 31, 31, 35, 35, 32, 37, 35, 39, 39, 32, 37, 35, 30, 36, 38, 36, 38, 39, 41, 36, 38, 36, 24, 34 },
                        null),
                makeCd(lane, tile, 1070, 1062, false,
                        new byte[]{ A,  P,  T,  A,  P,  P,  P,  C,  A,  G,  G,  G,  A,  T,  A,  A,  C,  T,  C,  A,  A,  P,  P,  G,  C,  P,  G,  G,  T,  C,  G,  C,  A,  A,  A,  C,  T,  C,  A,  G,  A,  C,  T,  A,  C,  A,  G,  G,  T,  G,  C,  A,  G,  A,  A,  A,  A,  C },
                        new byte[]{ 27,  2, 15, 31,  2,  2,  2, 19, 17, 31, 30, 31, 30, 30, 30, 30, 30, 30, 30, 30, 30,  2,  2,  2,  2,  2, 2,  2,  2, 35, 37, 37, 35, 31, 31, 34, 35, 35, 33, 37, 37, 37, 39, 39, 39, 39, 41, 40, 37, 39, 38, 38, 38, 39, 36, 39, 38, 38 },
                        null),
                makeCd(lane, tile, 1246, 1085, true,
                        new byte[]{ T,  C,  T,  T,  C,  C,  C,  A,  T,  G,  A,  G,  G,  G,  C,  A,  C,  A,  G,  T,  T,  T,  G,  A,  C,  P,  G,  G,  T,  C,  G,  C,  A,  A,  G,  T,  C,  G,  A,  A,  T,  T,  G,  T,  A,  A,  T,  T,  C,  C,  A,  T,  T,  T,  G,  C,  C,  C },
                        new byte[]{ 34, 34, 34, 37, 37, 37, 37, 37, 39, 39, 39, 39, 39, 40, 39, 39, 40, 40, 39, 38, 39, 40, 41, 39, 40,  2, 16, 28, 33, 35, 35, 37, 37, 34, 34, 34, 37, 37, 37, 37, 37, 39, 39, 39, 39, 39, 40, 41, 41, 39, 40, 41, 41, 41, 41, 40, 41, 40 },
                        null)
        ));
        goldSizes.put(ltStr(lane, tile), 20);
    }

    public static Map<Integer, ClusterData> clusterData(final int lane, final List<Integer> tiles, final String readStructure, final IlluminaDataType ... dataTypes) {
        final List<Integer> sortedTiles = new ArrayList<Integer>(tiles);
        Collections.sort(sortedTiles);

        final Map<Integer, ClusterData> data = new HashMap<Integer, ClusterData>();
        int offset = 0;
        for(final int tile : sortedTiles) {
            final String key = ltStr(lane,tile);
            final List<ClusterData> cds = goldData.get(key);
            final List<Integer> readNos = goldIndices.get(key);
            final int size = goldSizes.get(key);
            
            for(int i = 0; i < cds.size(); i++) {
                data.put(offset + readNos.get(i), selectiveCopyCd(cds.get(i), readStructure, dataTypes));
            }

            offset += size;
        }
        return data;
    }

    public static ReadData [] copyReadData(final ReadStructure rs, final IlluminaDataType [] dts, final ClusterData toCopy) {
        boolean doBases = false;
        boolean doQuals = false;
        boolean doInts  = false;
        boolean doNoise = false;

        for(final IlluminaDataType dt : dts) {
            switch(dt) {
                case BaseCalls:
                    doBases = true;
                    break;

                case QualityScores:
                    doQuals = true;
                    break;

                case RawIntensities:
                    doInts = true;
                    break;

                case Noise:
                    doNoise = true;
                    break;
            }
        }

        if(!doBases && !doQuals && !doInts && !doNoise)
            return null;

        final ReadData rdToCopy = toCopy.getRead(0); //Only gonna be one read in this
        final ReadData [] rds = new ReadData[rs.descriptors.size()];

        int index = 0;
        int baseIndex = 0;
        for(final ReadDescriptor readDesc : rs.descriptors) {
            final ReadData curRead = new ReadData(readDesc.type);
            if(doBases) {
                final byte [] bases = Arrays.copyOfRange(rdToCopy.getBases(), baseIndex, baseIndex + readDesc.length);
                curRead.setBases(bases);
            }

            if(doQuals) {
                final byte [] quals = Arrays.copyOfRange(rdToCopy.getQualities(), baseIndex, baseIndex + readDesc.length);
                curRead.setQualities(quals);
            }

            if(doInts) {
                final FourChannelIntensityData fcid = copyIntensities(rdToCopy.getRawIntensities(), baseIndex, readDesc.length);
                curRead.setRawIntensities(fcid);
            }

            if(doNoise) {
                final FourChannelIntensityData fcid = copyIntensities(rdToCopy.getNoise(), baseIndex, readDesc.length);
                curRead.setNoise(fcid);
            }

            baseIndex += readDesc.length;
            rds[index++] = curRead;
        }

        return rds;
    }

    private static FourChannelIntensityData copyIntensities(final FourChannelIntensityData toCopy, final int start, final int length) {
        final FourChannelIntensityData fcid = new FourChannelIntensityData(length);

        System.arraycopy(toCopy.getA(), start, fcid.getA(), 0, length);
        System.arraycopy(toCopy.getC(), start, fcid.getC(), 0, length);
        System.arraycopy(toCopy.getG(), start, fcid.getG(), 0, length);
        System.arraycopy(toCopy.getT(), start, fcid.getT(), 0, length);
        return fcid;
    }

    public static ClusterData selectiveCopyCd(final ClusterData toCopy, final String readStructure, final IlluminaDataType ... dataTypes) {
        final ReadStructure rs = new ReadStructure(readStructure);
        final ReadData [] rd = copyReadData(rs, dataTypes, toCopy);
        final ClusterData cd = new ClusterData(rd);

        for(final IlluminaDataType idt : dataTypes) {
            switch(idt) {
                case Position:
                    cd.setTile(toCopy.getTile());
                    cd.setLane(toCopy.getLane());
                    cd.setX(toCopy.getX());
                    cd.setY(toCopy.getY());
                    break;

                case PF:
                    cd.setPf(toCopy.isPf());
                    break;

                case Barcodes:
                    cd.setMatchedBarcode(toCopy.getMatchedBarcode());
                    break;

                default:
                    break;
            }
        }

        return cd;
    }

    public static ClusterData makeCd(final int lane, final int tile, final int xCoord, final int yCoord, final boolean pf, final byte [] bases, final byte[] qualities, final String matchedBarcode) {
        final ReadData rd = new ReadData();
        rd.setBases(Arrays.copyOf(bases, bases.length));
        rd.setQualities(Arrays.copyOf(qualities, bases.length));
        rd.setReadType(ReadType.T); //This will be ignored, as the cluster will be chopped up by ReadStructure

        final ClusterData cd = new ClusterData(rd);
        cd.setLane(lane);
        cd.setTile(tile);
        cd.setX(xCoord);
        cd.setY(yCoord);
        cd.setPf(pf);
        cd.setMatchedBarcode(matchedBarcode);

        return cd;
    }

}
