package picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.math3.util.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static htsjdk.samtools.SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS;

@CommandLineProgramProperties(
        summary = StabilizeQualityScores.USAGE_SUMMARY + StabilizeQualityScores.USAGE_DETAILS,
        oneLineSummary = StabilizeQualityScores.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
public class StabilizeQualityScores extends CommandLineProgram {

    static final String BIN_THRESHOLD_STRATEGY_DESCRIPTION = "Bins qualities >= N as QN, and qualities below it as Q2.";
    static final String BIN_BINS_STRATEGY_DESCRIPTION = "Bins qualities to the closest bin in the list.";

    static final String STABILIZE_DROP_STRATEGY_DESCRIPTION = "Stabilize by dropping runs of <=N quals down to the lowest value seen in the oscillation";
    static final String STABILIZE_AVERAGE_STRATEGY_DESCRIPTION = "Stabilize by setting runs of <=N quals to the average value seen in the oscillation";
    static final String STABILIZE_MERGE_STRATEGY_DESCRIPTION = "Stabilize by setting runs of <=N quals to the qual immediately before the oscillation";

    static final String USAGE_SUMMARY = "Stabilize quality scores.";
    static final String USAGE_DETAILS = "This tool bins and stabilizes the quality scores of a SAM or BAM file " +
            "to eliminate short runs of zigzaggy qualities.\n" +
            "Binning strategies available:\n\n" +
            "THRESHOLD N:\n " + BIN_THRESHOLD_STRATEGY_DESCRIPTION + "\n" +
            "BINS Q1 Q2 Q3 ...:\n " + BIN_BINS_STRATEGY_DESCRIPTION + "\n" +

            "\n" +

            "Stabilization strategies available:\n\n" +
            "DROP N:\n " + STABILIZE_DROP_STRATEGY_DESCRIPTION + "\n" +
            "AVERAGE N:\n " + STABILIZE_AVERAGE_STRATEGY_DESCRIPTION + "\n" +
            "MERGE N:\n " + STABILIZE_MERGE_STRATEGY_DESCRIPTION + "\n" +

            "<h3>Usage examples:</h3>\n" +
            "<h4>Stabilize file, using a Q30 threshold and dropping runs of oscillations of length &lt; 4</h4>\n"+
            "\n"+
            "java -jar picard.jar StabilizeQualityScores \\\n" +
            "      I=input.bam \\\n" +
            "      O=stabilized.bam \\\n" +
            "      THRESHOLD=30\n"+
            "      DROP=4\n";

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to downsample.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The stabilized SAM or BAM file to write.")
    public File OUTPUT;

    @Argument(shortName = "BANKS", optional = true, doc = "For use with THRESHOLD. Sets the above-threshold value.")
    public Integer THRESHOLD_UP = null;

    @Argument(shortName = "T", doc = "Any qual <T becomes Q2. Any qual >=T becomes T, or the value of THRESHOLD_UP if specified.", mutex = {"BINS", "CUTOFF"})
    public Integer THRESHOLD = 20;

    @Argument(shortName = "B", doc = "All qualities round to nearest bin in this list.", mutex = {"THRESHOLD", "CUTOFF"})
    public List<Integer> BINS = null;

    @Argument(shortName = "C", doc = "Find the cutoff where the end of the read drops consistently at or below this qual. Set the end of the read to Q2 and the rest to Q-average.", mutex = {"THRESHOLD", "BINS", "DROP", "AVERAGE", "MERGE"})
    public Integer CUTOFF = null;

    @Argument(shortName = "D", doc = "Drop runs of <= N quals to the lowest qual in the run.", mutex = {"AVERAGE", "MERGE", "CUTOFF"})
    public Integer DROP = null;

    @Argument(shortName = "A", doc = "Average runs of <= N quals to the binned average quality of the run.", mutex = {"DROP", "MERGE", "CUTOFF"})
    public Integer AVERAGE = 4;

    @Argument(shortName = "M", doc = "Merge runs of <= N quals to the qual immediately previous to the run.", mutex = {"AVERAGE", "DROP", "CUTOFF"})
    public Integer MERGE = null;

    abstract static class Binner {
        public abstract int[] bin(int[] quals);
    }

    static class ThresholdBinner extends Binner {
        protected Integer threshold;
        protected Integer thresholdUp;
        ThresholdBinner(Integer threshold, Integer thresholdUp) {
            this.threshold = threshold;
            this.thresholdUp = thresholdUp != null ? thresholdUp : threshold;
        }

        public int[] bin(int[] quals) {
            int[] iquals = new int[quals.length];
            for (int i = 0; i < quals.length; ++i) {
                iquals[i] = quals[i] >= this.threshold ? this.thresholdUp : 2;
            }
            return iquals;
        }
    }

    static class NearestBinner extends Binner {
        protected List<Integer> bins;
        NearestBinner(List<Integer> bins) { Collections.sort(bins); this.bins = bins; }

        public int[] bin(int[] quals) {
            int[] iquals = new int[quals.length];
            for (int i = 0; i < quals.length; ++i) {
                final int q = quals[i];
                iquals[i] = bins.stream()
                        .min(Comparator.comparingInt(n -> Math.abs(n - q)))
                        .orElseThrow(() -> new NoSuchElementException("No value present"));
            }
            return iquals;
        }
    }

    protected Binner makeBinner(Integer threshold, Integer thresholdUp, List<Integer> bins, Integer cutoff) {
        if(cutoff != null) {
            return null;
        }
        else if(threshold != null) {
            return new ThresholdBinner(threshold, thresholdUp);
        } else {
            return new NearestBinner(bins);
        }
    }

    static class RLEElem {
        int qual; int count;
        RLEElem(int q, int c) {
            qual = q; count = c;
        }

        @Override
        public boolean equals(Object obj) {
            if(obj.getClass() == this.getClass()) {
                RLEElem thatRLE = (RLEElem) obj;
                return thatRLE.qual == this.qual && thatRLE.count == this.count;
            }
            return false;
        }
    }

    static ArrayList<RLEElem> toRLE(int[] iquals) {
        ArrayList<RLEElem> rle = new ArrayList<>();

        RLEElem last = null;

        for (int iqual: iquals) {
            if( last != null && last.qual == iqual) {
                last.count++;
            } else {
                last = new RLEElem(iqual, 1);
                rle.add(last);
            }
        }
        return rle;
    }

    static byte[] toQuals(List<RLEElem> rles) {
        ArrayList<Byte> byteList = new ArrayList<>();

        for(RLEElem rle: rles) {
            byteList.addAll( Collections.nCopies(rle.count, (byte) rle.qual) );
        }

        byte[] bytes = new byte[byteList.size()];

        for(int i = 0; i < byteList.size(); ++i) {
            bytes[i] = byteList.get(i);
        }
        return bytes;
    }

    abstract static class Stabilizer {
        public Integer minRunLength;
        Stabilizer(Integer rl) { this.minRunLength = rl; }

        //Return a list of prevRle + whatever your deoscillation does.
        //prevRle may be null.
        public abstract List<RLEElem> stabilize(List<Integer> oquals, RLEElem prevRle, List<RLEElem> deoscillate);
    }

    /* merges RLEs with the previous qual in the list */
    static class MergeStabilizer extends Stabilizer {
        MergeStabilizer(Integer rl) { super(rl); }

        public List<RLEElem> stabilize(List<Integer> oquals, RLEElem prevRle, List<RLEElem> deoscillate) {
            int thisCount = deoscillate.stream().mapToInt(e -> e.count).sum();

            RLEElem newRLE = new RLEElem(0, 0);
            if(prevRle == null) {
                newRLE.qual = deoscillate.get(0).qual;
            } else {
                newRLE.qual = prevRle.qual;
                newRLE.count = prevRle.count;
            }
            newRLE.count += thisCount;

            List<RLEElem> outList = Arrays.asList(newRLE);
            return outList;
        }
    }

    /* sets the qual for the RLE to the minimum value seen */
    static class DropStabilizer extends Stabilizer {
        DropStabilizer(Integer rl) { super(rl); }

        //Return a list of prevRle + the deoscillation's minimum qual.
        //prevRle may be null.
        public List<RLEElem> stabilize(List<Integer> oquals, RLEElem prevRle, List<RLEElem> deoscillate) {
            int thisCount = deoscillate.stream().mapToInt(e -> e.count).sum();
            int minQual = deoscillate.stream().mapToInt(e -> e.qual).min().getAsInt();
            List<RLEElem> outList = new ArrayList<>();
            if(prevRle != null) {
                outList.add(prevRle);
            }
            outList.add(new RLEElem(minQual, thisCount));
            return outList;
        }
    }

    /* sets the qual for the RLE to the binned average value */
    static class AverageStabilizer extends Stabilizer {
        protected Binner binner;
        AverageStabilizer(Integer rl, Binner binner) { super(rl); this.binner = binner; }

        //Return a list of prevRle + whatever your deoscillation does.
        //prevRle may be null.
        public List<RLEElem> stabilize(List<Integer> oquals, RLEElem prevRle, List<RLEElem> deoscillate) {
            int thisCount = deoscillate.stream().mapToInt(e -> e.count).sum();
            int avgQual = (int) Math.round(oquals.stream().mapToInt(i->i).average().getAsDouble());
            int binnedAvgQual = binner.bin(new int[] {avgQual})[0]; //lazy much?
            List<RLEElem> outList = new ArrayList<>();
            if(prevRle != null) {
                outList.add(prevRle);
            }

            outList.add(new RLEElem(binnedAvgQual, thisCount));
            return outList;
        }
    }

    protected Stabilizer makeStabilizer(Binner binner, Integer drop, Integer average, Integer merge, Integer cutoff) {
        if(cutoff != null) {
            return null;
        } else if(merge != null) {
            return new MergeStabilizer(merge);
        } else if(average != null) {
            return new AverageStabilizer(average, binner);
        } else {
            return new DropStabilizer(drop);
        }
    }

    private int[] byte2intQuals(byte[] bquals) {
        int[] iquals = new int[bquals.length];
        for (int i = 0; i < bquals.length; ++i) {
            iquals[i] = Byte.toUnsignedInt(bquals[i]);
        }
        return iquals;
    }

    public static List<Pair<Boolean, List<RLEElem>>> groupify(List<RLEElem> rles, int minRunLength) {
        ArrayList<Pair<Boolean, List<RLEElem>>> groups = new ArrayList<>();

        Pair<Boolean, List<RLEElem>> last = null;

        for (RLEElem rle: rles) {
            boolean shouldStabilize = rle.count <= minRunLength;
            if( last != null && last.getFirst() == shouldStabilize) {
                last.getSecond().add(rle);
            } else {
                ArrayList<RLEElem> newElem = new ArrayList<>();
                newElem.add(rle);
                last = new Pair<>(shouldStabilize, newElem);
                groups.add(last);
            }
        }
        return groups;
    }

    protected LinkedList<RLEElem> binAndStabilize(Binner binner, Stabilizer stabilizer, int[] iquals) {
        LinkedList<RLEElem> outputRle = new LinkedList<>();
        //bin
        int[] binnedQuals = binner.bin(iquals);

        //rle
        ArrayList<RLEElem> rle = toRLE(binnedQuals);

        List<Integer> sliceableOQuals = Arrays.stream(iquals).boxed().collect(Collectors.toList());
        int numStabilized = 0;
        //groupify returns pairs of "should we stabilize this run" and "the run"
        for (Pair<Boolean, List<RLEElem>> pair : groupify(rle, stabilizer.minRunLength)) {
            int rleLen = pair.getSecond().stream().mapToInt(e -> e.count).sum();

            if (pair.getFirst()) {
                RLEElem lastRle = outputRle.size() >= 1 ? outputRle.removeLast() : null;
                List<Integer> oquals = sliceableOQuals.subList(numStabilized, numStabilized + rleLen);

                List<RLEElem> addMe = stabilizer.stabilize(oquals, lastRle, pair.getSecond());
                outputRle.addAll(addMe);
            } else {
                outputRle.addAll(pair.getSecond());
            }
            numStabilized += rleLen;
        }
        return outputRle;
    }

    //return the index of the first item in the list that is above the cutoff
    protected static int getIndexAboveCutoff(List<Integer> ls, Integer cutoff) {
        for(int i = 0; i < ls.size(); ++i) {
            if( ls.get(i) > cutoff ) {
                return i;
            }
        }
        return ls.size();
    }

    protected static List<RLEElem> cutoff(int[] iquals, Integer cutoff, boolean readIsReversed) {
        List<Integer> is = Arrays.stream(iquals).boxed().collect(Collectors.toList());
        if(!readIsReversed) {
            Collections.reverse(is);
        }

        int idx = getIndexAboveCutoff(is, cutoff);
        int avg = (int) Math.round(is.subList(idx, iquals.length).stream().mapToInt(i->i).average().orElse(0.0));

        if(idx == 0) {
            return Arrays.asList(new RLEElem(avg, iquals.length));
        }
        else if(idx == iquals.length) {
            return Arrays.asList(new RLEElem(2, iquals.length));
        }
        else if(!readIsReversed) {
            return Arrays.asList(new RLEElem(avg, iquals.length - idx), new RLEElem(2, idx));
        } else {
            return Arrays.asList(new RLEElem(2, idx), new RLEElem(avg, iquals.length - idx));
        }
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final Binner binner = makeBinner(THRESHOLD, THRESHOLD_UP, BINS, CUTOFF);
        final Stabilizer stabilizer = makeStabilizer(binner, DROP, AVERAGE, MERGE, CUTOFF);

        SamReaderFactory readerFactory = SamReaderFactory.makeDefault().setOption(INCLUDE_SOURCE_IN_RECORDS, true).validationStringency(ValidationStringency.SILENT);
        SAMFileWriterFactory writerFactory = new SAMFileWriterFactory().setCreateIndex(true);

        try (final SamReader samReader = readerFactory.open(SamInputResource.of(INPUT)) ) {
            try(final SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(samReader.getFileHeader(), true, OUTPUT)) {
                for (SAMRecord record : samReader) {

                    int[] iquals = byte2intQuals(record.getBaseQualities());
                    List<RLEElem> outputRle;

                    if(CUTOFF != null) {
                        boolean readIsReversed = record.getReadNegativeStrandFlag();
                        outputRle = cutoff(iquals, CUTOFF, readIsReversed);
                    } else {
                        outputRle = binAndStabilize(binner, stabilizer, iquals);
                    }

                    record.setBaseQualities(toQuals(outputRle));
                    writer.addAlignment(record);
                }
            }
        } catch (IOException e) {
        System.out.println( "A problem occurred: " + e.getMessage());
        return 1;
    }

        return 0;
    }
}
