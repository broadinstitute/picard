package picard.util;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.None;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.regex.Pattern;

/**
 * Designs baits for hybrid selection!
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage = "Designs baits or oligos for hybrid selection reactions.",
        usageShort = "Designs baits or oligos for hybrid selection reactions.",
        programGroup = None.class
)
public class BaitDesigner extends CommandLineProgram {
    /**
     * Subclass of Interval for representing Baits, that caches the bait sequence.
     */
    static class Bait extends Interval {
        byte[] bases;

        /** Pass through constructor. */
        public Bait(final String sequence, final int start, final int end, final boolean negative, final String name) {
            super(sequence, start, end, negative, name);
        }

        /** Method that takes in the reference sequence this bait is on and caches the bait's bases. */
        public void addBases(final ReferenceSequence reference, final boolean useStrandInfo) {
            final byte[] tmp = new byte[length()];
            System.arraycopy(reference.getBases(), getStart() - 1, tmp, 0, length());

            if (useStrandInfo && isNegativeStrand()) {
                SequenceUtil.reverseComplement(tmp);
            }

            setBases(tmp);
        }

        public int getMaskedBaseCount() {
            return BaitDesigner.getMaskedBaseCount(bases, 0, bases.length);
        }

        @Override
        public String toString() {
            return "Bait{" +
                    "name=" + getName() +
                    ", bases=" + StringUtil.bytesToString(bases) +
                    '}';
        }

        public void setBases(final byte[] bases) { this.bases = bases;}

        public byte[] getBases() { return bases; }
    }

    /**
     * Set of possible design strategies for bait design.
     */
    public enum DesignStrategy {
        /** Implementation that "constrains" baits to be within the target region when possible. */
        CenteredConstrained {
            List<Bait> design(final BaitDesigner designer, final Interval target, final ReferenceSequence reference) {
                final List<Bait> baits = new LinkedList<Bait>();

                final int baitSize = designer.BAIT_SIZE;
                final int baitOffset = designer.BAIT_OFFSET;

                // Treat targets less than a bait length differently
                if (target.length() <= baitSize) {
                    final int midpoint = target.getStart() + (target.length() / 2);
                    final int baitStart = midpoint - (baitSize / 2);
                    final Bait bait = new Bait(target.getSequence(),
                            baitStart,
                            CoordMath.getEnd(baitStart, baitSize),
                            target.isNegativeStrand(),
                            designer.makeBaitName(target.getName(), 1, 1));
                    bait.addBases(reference, designer.DESIGN_ON_TARGET_STRAND);
                    baits.add(bait);
                } else {
                    // Work out how many baits we need and how to space them
                    final int baitCount = 1 + (int) Math.ceil((target.length() - baitSize) / (double) baitOffset);
                    final int firstBaitStart = target.getStart();
                    final int lastBaitStart = CoordMath.getStart(target.getEnd(), baitSize);
                    final double actualShift = (lastBaitStart - firstBaitStart) / (double) (baitCount - 1);

                    // And then design them
                    int baitIndex = 1;
                    int start = firstBaitStart;
                    while (start <= lastBaitStart) {
                        final int end = CoordMath.getEnd(start, baitSize);
                        final Bait bait = new Bait(target.getSequence(),
                                start,
                                end,
                                target.isNegativeStrand(),
                                designer.makeBaitName(target.getName(), baitIndex, baitCount));
                        bait.addBases(reference, designer.DESIGN_ON_TARGET_STRAND);
                        baits.add(bait);

                        // Recalculate from actualShift to avoid compounding rounding errors
                        start = firstBaitStart + (int) Math.round(actualShift * baitIndex);
                        ++baitIndex;
                    }
                }

                return baits;
            }
        },

        /**
         * Design that places baits at fixed offsets over targets, allowing them to hang off the ends
         * as dictated by the target size and offset.
         */
        FixedOffset {
            List<Bait> design(final BaitDesigner designer, final Interval target, final ReferenceSequence reference) {
                final List<Bait> baits = new LinkedList<Bait>();

                final int baitSize = designer.BAIT_SIZE;
                final int baitOffset = designer.BAIT_OFFSET;
                final int minTargetSize = baitSize + (baitOffset * (designer.MINIMUM_BAITS_PER_TARGET - 1));

                // Redefine the target to be the size of a bait, or the size that a number of baits
                // tiled tiles across the target at the fixed offset
                final Interval t2;
                if (target.length() < minTargetSize) {
                    final int addon = minTargetSize - target.length();
                    final int left = addon / 2;
                    final int right = addon - left;
                    t2 = new Interval(target.getSequence(),
                            Math.max(target.getStart() - left, 1),
                            Math.min(target.getEnd() + right, reference.length()),
                            target.isNegativeStrand(),
                            target.getName());
                } else {
                    t2 = target;
                }

                // Work out how many baits we need and how to space them
                final int baitCount = 1 + (int) Math.ceil((t2.length() - baitSize) / (double) baitOffset);
                final int baitedBases = baitSize + (baitOffset * (baitCount - 1));
                final int firstBaitStart = Math.max(t2.getStart() - ((baitedBases - t2.length()) / 2), 1);

                // And then design them
                final byte[] bases = reference.getBases();
                final int MAX_MASKED = designer.REPEAT_TOLERANCE;

                for (int i = 1; i <= baitCount; ++i) {
                    int start = firstBaitStart + (baitOffset * (i - 1));
                    int end = CoordMath.getEnd(start, baitSize);

                    if (end > reference.length()) break;

                    // If there are too many soft masked bases try shifting it around
                    if (getMaskedBaseCount(bases, start - 1, end) > MAX_MASKED) {
                        final int maxMove = baitOffset * 3 / 4;

                        for (int move = 1; move <= maxMove; move++) {
                            // Move it "backwards?"
                            if (start - move >= 1 && getMaskedBaseCount(bases, start - move - 1, end - move) <= MAX_MASKED) {
                                start = start - move;
                                end = end - move;
                                break;
                            }
                            // Move it "forwards"?
                            if (end + move <= reference.length() && getMaskedBaseCount(bases, start + move - 1, end + move) <= MAX_MASKED) {
                                start = start + move;
                                end = end + move;
                                break;
                            }
                        }
                    }

                    final Bait bait = new Bait(t2.getSequence(),
                            start,
                            end,
                            t2.isNegativeStrand(),
                            designer.makeBaitName(t2.getName(), i, baitCount));
                    bait.addBases(reference, designer.DESIGN_ON_TARGET_STRAND);
                    baits.add(bait);
                }

                return baits;
            }
        },

        /**
         * Ultra simple bait design algorithm that just lays down baits starting at the target start position
         * until either the bait start runs off the end of the target or the bait would run off the sequence
         */
        Simple {
            @Override
            List<Bait> design(final BaitDesigner designer, final Interval target, final ReferenceSequence reference) {
                final List<Bait> baits = new LinkedList<Bait>();
                final int baitSize = designer.BAIT_SIZE;
                final int baitOffset = designer.BAIT_OFFSET;
                final int lastPossibleBaitStart = Math.min(target.getEnd(), reference.length() - baitSize);
                final int baitCount = 1 + (int) Math.floor((lastPossibleBaitStart - target.getStart()) / (double) baitOffset);

                int i = 0;
                for (int start = target.getStart(); start < lastPossibleBaitStart; start += baitOffset) {
                    final Bait bait = new Bait(target.getSequence(),
                            start,
                            CoordMath.getEnd(start, baitSize),
                            target.isNegativeStrand(),
                            designer.makeBaitName(target.getName(), ++i, baitCount));
                    bait.addBases(reference, designer.DESIGN_ON_TARGET_STRAND);
                    baits.add(bait);
                }
                return baits;
            }
        };

        /** Design method that each Design Strategy must implement. */
        abstract List<Bait> design(BaitDesigner designer, Interval target, ReferenceSequence reference);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Options for the Bait Designer
    ///////////////////////////////////////////////////////////////////////////

    @Option(shortName = "T", doc = "The file with design parameters and targets")
    public File TARGETS;

    @Option(doc = "The name of the bait design")
    public String DESIGN_NAME;

    @Option(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc = "The reference sequence fasta file")
    public File REFERENCE_SEQUENCE;

    @Option(doc = "The left amplification primer to prepend to all baits for synthesis")
    public String LEFT_PRIMER = "ATCGCACCAGCGTGT";

    @Option(doc = "The right amplification primer to prepend to all baits for synthesis")
    public String RIGHT_PRIMER = "CACTGCGGCTCCTCA";

    @Option(doc = "The design strategy to use to layout baits across each target")
    public DesignStrategy DESIGN_STRATEGY = DesignStrategy.FixedOffset;

    @Option(doc = "The length of each individual bait to design")
    public int BAIT_SIZE = 120;

    @Option(doc = "The minimum number of baits to design per target.")
    public int MINIMUM_BAITS_PER_TARGET = 2;

    @Option(doc = "The desired offset between the start of one bait and the start of another bait for the same target.")
    public int BAIT_OFFSET = 80;

    @Option(doc = "Pad the input targets by this amount when designing baits. Padding is applied on both sides in this amount.")
    public int PADDING = 0;

    @Option(doc = "Baits that have more than REPEAT_TOLERANCE soft or hard masked bases will not be allowed")
    public int REPEAT_TOLERANCE = 50;

    @Option(doc = "The size of pools or arrays for synthesis. If no pool files are desired, can be set to 0.")
    public int POOL_SIZE = 55000;

    @Option(doc = "If true, fill up the pools with alternating fwd and rc copies of all baits. Equal copies of " +
            "all baits will always be maintained")
    public boolean FILL_POOLS = true;

    @Option(doc = "If true design baits on the strand of the target feature, if false always design on the + strand of " +
            "the genome.")
    public boolean DESIGN_ON_TARGET_STRAND = false;

    @Option(doc = "If true merge targets that are 'close enough' that designing against a merged target would be more efficient.")
    public boolean MERGE_NEARBY_TARGETS = true;

    @Option(doc = "If true also output .design.txt files per pool with one line per bait sequence")
    public boolean OUTPUT_AGILENT_FILES = true;

    @Option(shortName = "O", optional = true,
            doc = "The output directory. If not provided then the DESIGN_NAME will be used as the output directory")
    public File OUTPUT_DIRECTORY;

    // "Output" members that will also get picked up by writeParametersFile()
    int TARGET_TERRITORY;
    int TARGET_COUNT;
    int BAIT_TERRITORY;
    int BAIT_COUNT;
    int BAIT_TARGET_TERRITORY_INTERSECTION;
    int ZERO_BAIT_TARGETS;
    double DESIGN_EFFICIENCY;

    // Utility objects
    private static final Log log = Log.getInstance(BaitDesigner.class);
    private final NumberFormat fmt = NumberFormat.getIntegerInstance();

    /** Takes a target name and a bait index and creates a uniform bait name. */
    String makeBaitName(final String targetName, final int baitIndex, final int totalBaits) {
        final String total = fmt.format(totalBaits);
        String bait = fmt.format(baitIndex);

        // Pad out the bait number to match the longest one for this target
        while (bait.length() < total.length()) bait = "0" + bait;

        return targetName + "_bait#" + bait;
    }

    /** Returns the total of soft or hard masked bases in the interval of bases. */
    public static int getMaskedBaseCount(final byte[] bases, final int from, final int until) {
        int count = 0;
        for (int i = from; i < until; i++) {
            final byte b = bases[i];
            if (b != 'A' && b != 'C' && b != 'G' && b != 'T') ++count;
        }

        return count;
    }

    /** Stock main method. */
    public static void main(final String[] args) {
        new BaitDesigner().instanceMainWithExit(args);
    }

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errors = new ArrayList<String>();

        final Pattern p = Pattern.compile("^[ACGTacgt]*$");
        if (LEFT_PRIMER != null && !p.matcher(LEFT_PRIMER).matches()) {
            errors.add("Left primer " + LEFT_PRIMER + " is not a valid primer sequence.");
        }

        if (RIGHT_PRIMER != null && !p.matcher(RIGHT_PRIMER).matches()) {
            errors.add("Right primer " + RIGHT_PRIMER + " is not a valid primer sequence.");
        }

        if (errors.size() > 0) return errors.toArray(new String[errors.size()]);
        else return null;
    }

    int estimateBaits(final int start, final int end) {
        final int length = end - start + 1;
        return Math.max(MINIMUM_BAITS_PER_TARGET, (int) (Math.ceil(length - BAIT_SIZE) / (double) BAIT_OFFSET) + 1);
    }

    /**
     * Main method that coordinates the checking of inputs, designing of baits and then
     * the writing out of all necessary files.
     *
     * @return
     */
    @Override
    protected int doWork() {
        // Input parameter munging and checking
        if (OUTPUT_DIRECTORY == null) OUTPUT_DIRECTORY = new File(DESIGN_NAME);

        IOUtil.assertFileIsReadable(TARGETS);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);

        if (!OUTPUT_DIRECTORY.exists()) {
            OUTPUT_DIRECTORY.mkdirs();
        }
        IOUtil.assertDirectoryIsWritable(OUTPUT_DIRECTORY);

        // Load up the targets and the reference
        final IntervalList targets;
        final IntervalList originalTargets = IntervalList.fromFile(TARGETS);

        {
            // Apply padding
            final IntervalList padded = new IntervalList(originalTargets.getHeader());
            final SAMSequenceDictionary dict = padded.getHeader().getSequenceDictionary();
            for (final Interval i : originalTargets.getIntervals()) {
                padded.add(new Interval(i.getSequence(),
                        Math.max(i.getStart() - PADDING, 1),
                        Math.min(i.getEnd() + PADDING, dict.getSequence(i.getSequence()).getSequenceLength()),
                        i.isNegativeStrand(),
                        i.getName()));
            }

            log.info("Starting with " + padded.size() + " targets.");
            padded.uniqued();
            log.info("After uniquing " + padded.size() + " targets remain.");

            if (MERGE_NEARBY_TARGETS) {
                final ListIterator<Interval> iterator = padded.getIntervals().listIterator();
                Interval previous = iterator.next();

                targets = new IntervalList(padded.getHeader());

                while (iterator.hasNext()) {
                    final Interval next = iterator.next();
                    if (previous.getSequence().equals(next.getSequence()) &&
                            estimateBaits(previous.getStart(), previous.getEnd()) + estimateBaits(next.getStart(), next.getEnd()) >=
                                    estimateBaits(previous.getStart(), next.getEnd())) {
                        previous = new Interval(previous.getSequence(),
                                previous.getStart(),
                                Math.max(previous.getEnd(), next.getEnd()),
                                previous.isNegativeStrand(),
                                previous.getName());
                    } else {
                        targets.add(previous);
                        previous = next;
                    }
                }

                if (previous != null) targets.add(previous);
                log.info("After collapsing nearby targets " + targets.size() + " targets remain.");
            } else {
                targets = padded;
            }
        }

        final ReferenceSequenceFileWalker referenceWalker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);

        // Check that the reference and the target list have matching headers
        SequenceUtil.assertSequenceDictionariesEqual(referenceWalker.getSequenceDictionary(),
                targets.getHeader().getSequenceDictionary());

        // Design the baits!
        int discardedBaits = 0;
        final IntervalList baits = new IntervalList(targets.getHeader());
        for (final Interval target : targets) {
            final int sequenceIndex = targets.getHeader().getSequenceIndex(target.getSequence());
            final ReferenceSequence reference = referenceWalker.get(sequenceIndex);

            for (final Bait bait : DESIGN_STRATEGY.design(this, target, reference)) {
                if (bait.length() != BAIT_SIZE) {
                    throw new PicardException("Bait designed at wrong length: " + bait);
                }

                if (bait.getMaskedBaseCount() <= REPEAT_TOLERANCE) {
                    baits.add(bait);

                    for (final byte b : bait.getBases()) {
                        final byte upper = StringUtil.toUpperCase(b);
                        if (upper != 'A' && upper != 'C' && upper != 'G' && upper != 'T') {
                            log.warn("Bait contains non-synthesizable bases: " + bait);
                        }
                    }
                } else {
                    log.debug("Discarding bait: " + bait);
                    discardedBaits++;
                }
            }
        }

        calculateStatistics(targets, baits);
        log.info("Designed and kept " + baits.size() + " baits, discarded " + discardedBaits);

        // Write out some files!
        originalTargets.write(new File(OUTPUT_DIRECTORY, DESIGN_NAME + ".targets.interval_list"));
        baits.write(new File(OUTPUT_DIRECTORY, DESIGN_NAME + ".baits.interval_list"));
        writeParametersFile(new File(OUTPUT_DIRECTORY, DESIGN_NAME + ".design_parameters.txt"));
        writeDesignFastaFile(new File(OUTPUT_DIRECTORY, DESIGN_NAME + ".design.fasta"), baits);
        if (POOL_SIZE > 0) writePoolFiles(OUTPUT_DIRECTORY, DESIGN_NAME, baits);

        return 0;
    }

    /** Calculates a few statistics about the bait design that can then be output. */
    void calculateStatistics(final IntervalList targets, final IntervalList baits) {
        this.TARGET_TERRITORY = (int) targets.getUniqueBaseCount();
        this.TARGET_COUNT = targets.size();
        this.BAIT_TERRITORY = (int) baits.getUniqueBaseCount();
        this.BAIT_COUNT = baits.size();
        this.DESIGN_EFFICIENCY = this.TARGET_TERRITORY / (double) this.BAIT_TERRITORY;

        // Figure out the intersection between all targets and all baits
        final IntervalList tmp = new IntervalList(targets.getHeader());
        final OverlapDetector<Interval> detector = new OverlapDetector<Interval>(0, 0);
        detector.addAll(baits.getIntervals(), baits.getIntervals());

        for (final Interval target : targets) {
            final Collection<Interval> overlaps = detector.getOverlaps(target);

            if (overlaps.isEmpty()) {
                this.ZERO_BAIT_TARGETS++;
            } else {
                for (final Interval i : overlaps) tmp.add(target.intersect(i));
            }
        }

        tmp.uniqued();
        this.BAIT_TARGET_TERRITORY_INTERSECTION = (int) tmp.getBaseCount();
    }

    /** Method that writes out all the parameter values that were used in the design using reflection. */
    void writeParametersFile(final File file) {
        try {
            final BufferedWriter out = IOUtil.openFileForBufferedWriting(file);
            for (final Field field : getClass().getDeclaredFields()) {
                if (Modifier.isPrivate(field.getModifiers())) continue;

                final String name = field.getName();

                if (name.toUpperCase().equals(name) && !name.equals("USAGE")) {
                    final Object value = field.get(this);

                    if (value != null) {
                        out.append(name);
                        out.append("=");
                        out.append(value.toString());
                        out.newLine();
                    }
                }
            }
            out.close();
        } catch (Exception e) {
            throw new PicardException("Error writing out parameters file.", e);
        }
    }

    void writeDesignFastaFile(final File file, final IntervalList baits) {
        final BufferedWriter out = IOUtil.openFileForBufferedWriting(file);
        for (final Interval i : baits) {
            writeBaitFasta(out, i, false);
        }
        CloserUtil.close(out);
    }

    /** Writes a Bait out in fasta format to an output BufferedWriter. */
    private void writeBaitFasta(final BufferedWriter out, final Interval i, final boolean rc) {
        try {
            final Bait bait = (Bait) i;
            out.append(">");
            out.append(bait.getName());
            out.newLine();

            final String sequence = getBaitSequence(bait, rc);
            out.append(sequence);
            out.newLine();
        } catch (IOException ioe) {
            throw new PicardException("Error writing out bait information.", ioe);
        }
    }

    /** Gets the bait sequence, with primers, as a String, RC'd as appropriate. */
    private String getBaitSequence(final Bait bait, final boolean rc) {
        String sequence = (LEFT_PRIMER == null ? "" : LEFT_PRIMER) +
                StringUtil.bytesToString(bait.getBases()) +
                (RIGHT_PRIMER == null ? "" : RIGHT_PRIMER);

        if (rc) sequence = SequenceUtil.reverseComplement(sequence);
        return sequence;
    }

    /**
     * Writes out fasta files for each pool and also agilent format files if requested.
     *
     * @param dir      the directory to output files into
     * @param basename the basename of each file
     * @param baits    the set of baits to write out
     */
    void writePoolFiles(final File dir, final String basename, final IntervalList baits) {
        final int copies;
        if (FILL_POOLS && baits.size() < POOL_SIZE) copies = (int) Math.floor(POOL_SIZE / (double) baits.size());
        else copies = 1;

        int written = 0;
        int nextPool = 0;
        BufferedWriter out = null;
        BufferedWriter agilentOut = null;
        final String prefix = DESIGN_NAME.substring(0, Math.min(DESIGN_NAME.length(), 8)) + "_"; // prefix for 15 digit bait id
        final NumberFormat fmt = new DecimalFormat("000000");

        try {
            for (int i = 0; i < copies; ++i) {
                final boolean rc = i % 2 == 1;

                int baitId = 1;
                for (final Interval interval : baits) {
                    final Bait bait = (Bait) interval;

                    if (written++ % POOL_SIZE == 0) {
                        if (out != null) out.close();
                        if (agilentOut != null) agilentOut.close();

                        final String filename = basename + ".pool" + nextPool++ + ".design.";
                        out = IOUtil.openFileForBufferedWriting(new File(dir, filename + "fasta"));
                        if (OUTPUT_AGILENT_FILES) {
                            agilentOut = IOUtil.openFileForBufferedWriting(new File(dir, filename + "txt"));
                        }
                    }

                    writeBaitFasta(out, interval, rc);
                    if (OUTPUT_AGILENT_FILES) {
                        agilentOut.append(prefix).append(fmt.format(baitId++));
                        agilentOut.append("\t");
                        agilentOut.append(getBaitSequence(bait, rc).toUpperCase());
                        agilentOut.newLine();
                    }
                }
            }

            CloserUtil.close(out);
            CloserUtil.close(agilentOut);
        } catch (Exception e) {
            throw new PicardException("Error while writing pool files.", e);
        }
    }
}
