/*
 * The MIT License
 *
 * Copyright (c) 2015-2016 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import picard.sam.markduplicates.MarkDuplicates;
import picard.sam.util.PhysicalLocationInt;
import picard.sam.util.ReadNameParser;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * <h3>Summary</h3>
 * Class to downsample a SAM/BAM file based on the position of the read in a flowcell. As with {@link DownsampleSam}, all the
 * reads with the same queryname are either kept or dropped as a unit.
 * <br/>
 *
 * <h3>Details</h3>
 * The downsampling is <b>not</b> random (and there is no random seed). It is deterministically determined by the position
 * of each read within its tile. Specifically, it draws an ellipse that covers a {@link #FRACTION} of the total tile's
 * area and of all the edges of the tile. It uses this area to determine whether to keep or drop the record. Since reads
 * with the same name have the same position (mates, secondary and supplemental alignments), the decision will be the
 * same for all of them. The main concern of this downsampling method is that due to "optical duplicates" downsampling
 * randomly can create a result that has a different optical duplicate rate, and therefore a different estimated library
 * size (when running MarkDuplicates). This method keeps (physically) close read together, so that (except
 * for reads near the boundary of the circle) optical duplicates are kept or dropped as a group.
 *
 * By default the program expects the read names to have 5 or 7 fields separated by colons (:), and it takes the last two
 * to indicate the x and y coordinates of the reads within the tile whence it was sequenced. See
 * {@link ReadNameParser#DEFAULT_READ_NAME_REGEX} for more detail. By default the program reads the {@link #INPUT} once,
 * learning the size of each tile as its reads arrive, and decides each read against its tile's bounds (single-pass mode).
 * Set {@link #TWO_PASS} to traverse the {@link #INPUT} twice instead: first to find out the size of each of the tiles, and
 * next to perform the downsampling, which is exact. Single-pass needs {@link #READ_COUNT}.
 *
 * Downsampling invalidates the duplicate flag because duplicate reads before downsampling may not all remain duplicated
 * after downsampling. Thus, the default setting also removes the duplicate information.
 *
 * <h3>Example</h3>
 * <pre>
 * java -jar picard.jar PositionBasedDownsampleSam \
 *       I=input.bam \
 *       O=downsampled.bam \
 *       FRACTION=0.1
 * </pre>
 * <h3>Caveats</h3>
 * <ol>
 * <li>
 * This method is <b>technology and read-name dependent</b>. If the read-names do not have coordinate information
 * embedded in them, or if your BAM contains reads from multiple technologies (flowcell versions, sequencing machines).
 * this will not work properly. It has been designed to work with Illumina technology and reads-names. Consider
 * modifying {@link #READ_NAME_REGEX} in other cases.
 * </li>
 * <li>
 * The code has been designed to simulate, as accurately as possible, sequencing less, <b>not</b> for getting an exact
 * downsampled fraction (Use {@link DownsampleSam} for that.) In particular, since the reads may be distributed non-evenly
 * within the lanes/tiles, the resulting downsampling percentage will not be accurately determined by the input argument
 * {@link #FRACTION}.
 * </li>
 * <li>
 * Consider running {@link MarkDuplicates} after downsampling in order to "expose" the duplicates whose representative has
 * been downsampled away.
 * </li>
 * <li>
 * The downsampling assumes a uniform distribution of reads in the flowcell.
 * Input already downsampled with PositionBasedDownsampleSam violates this assumption.
 * To guard against such input,
 * PositionBasedDownsampleSam always places a PG record in the header of its output,
 * and aborts whenever it finds such a PG record in its input.
 * </li>
 * </ol>
 *
 * @author Yossi Farjoun
 */
@CommandLineProgramProperties(
      summary = "<h3>Summary</h3>\n" +
              "Class to downsample a SAM/BAM/CRAM file based on the position of the read in a flowcell. As with DownsampleSam, all the " +
              "reads with the same queryname are either kept or dropped as a unit." +
              "\n\n " +
              "<h3>Details</h3>\n" +
              "The downsampling is _not_ random (and there is no random seed). It is deterministically determined by the position " +
              "of each read within its tile. Specifically, it draws an ellipse that covers a FRACTION of the total tile's " +
              "area and of all the edges of the tile. It uses this area to determine whether to keep or drop the record. Since reads " +
              "with the same name have the same position (mates, secondary and supplemental alignments), the decision will be the " +
              "same for all of them. The main concern of this downsampling method is that due to \"optical duplicates\" downsampling " +
              "randomly can create a result that has a different optical duplicate rate, and therefore a different estimated library " +
              "size (when running MarkDuplicates). This method keeps (physically) close read together, so that (except " +
              "for reads near the boundary of the circle) optical duplicates are kept or dropped as a group. " +
              "By default the program expects the read names to have 5 or 7 fields separated by colons (:) and it takes the last two " +
              "to indicate the x and y coordinates of the reads within the tile whence it was sequenced. See " +
              "DEFAULT_READ_NAME_REGEX for more detail. By default the program reads the INPUT once, learning the size of each " +
              "tile as its reads arrive (single-pass mode, which needs READ_COUNT). TWO_PASS=true traverses the INPUT twice " +
              "instead: first to find out the size of each of the tiles, and next to perform the downsampling, which is exact. " +
              "Downsampling invalidates the duplicate flag because duplicate reads before downsampling " +
              "may not all remain duplicated after downsampling. Thus, the default setting also removes the duplicate information. " +
              "\n\n" +
              "Example\n\n" +
              "java -jar picard.jar PositionBasedDownsampleSam \\\n" +
              "      I=input.bam \\\n" +
              "      O=downsampled.bam \\\n" +
              "      FRACTION=0.1\n" +
              "\n" +

              "Caveats\n\n" +
              "Note 1: " +
              "This method is <b>technology and read-name dependent</b>. If the read-names do not have coordinate information " +
              "embedded in them, or if your BAM contains reads from multiple technologies (flowcell versions, sequencing machines) " +
              "this will not work properly. It has been designed to work with Illumina technology and reads-names. Consider " +
              "modifying {@link #READ_NAME_REGEX} in other cases. " +
              "\n\n" +
              "Note 2: " +
              "The code has been designed to simulate, as accurately as possible, sequencing less, <b>not</b> for getting an exact " +
              "downsampled fraction (Use {@link DownsampleSam} for that.) In particular, since the reads may be distributed non-evenly " +
              "within the lanes/tiles, the resulting downsampling percentage will not be accurately determined by the input argument " +
              "FRACTION. " +
              "\n\n"+
              "Note 3:" +
              "Consider running {@link MarkDuplicates} after downsampling in order to \"expose\" the duplicates whose representative has " +
              "been downsampled away." +
              "\n\n" +
              "Note 4:" +
              "The downsampling assumes a uniform distribution of reads in the flowcell. " +
              "Input already downsampled with PositionBasedDownsampleSam violates this assumption. " +
              "To guard against such input, " +
              "PositionBasedDownsampleSam always places a PG record in the header of its output, " +
              "and aborts whenever it finds such a PG record in its input.",
        oneLineSummary = "Downsample a SAM or BAM file to retain a subset of the reads based on the reads location in each tile in the flowcell.",
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature
public class PositionBasedDownsampleSam extends CommandLineProgram {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM/BAM/CRAM file to downsample.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output, downsampled, SAM/BAM/CRAM file.")
    public File OUTPUT;

    @Argument(shortName = "F", doc = "The (approximate) fraction of reads to be kept, between 0 and 1.")
    public Double FRACTION = null;

    @Argument(doc = "Determines whether the duplicate tag should be reset since the downsampling requires re-marking duplicates.")
    public boolean REMOVE_DUPLICATE_INFORMATION = true;

    @Argument(doc = "Use these regular expressions to parse read names in the input SAM file. Read names are " +
            "parsed to extract three variables: tile/region, x coordinate and y coordinate. The x and y coordinates are used " +
            "to determine the downsample decision. " +
            "Set this option to null to disable optical duplicate detection, e.g. for RNA-seq " +
            "The regular expression should contain three capture groups for the three variables, in order. " +
            "It must match the entire read name. " +
            "Note that if the default regex is specified, a regex match is not actually done, but instead the read name " +
            "is split on colons (:). " +
            "For 5 element names, the 3rd, 4th and 5th elements are assumed to be tile, x and y values. " +
            "For 7 element names (CASAVA 1.8), the 5th, 6th, and 7th elements are assumed to be tile, x and y values. " +
            "Must be left at this default when READ_NAME_FORMAT is set.")
    public String READ_NAME_REGEX = ReadNameParser.DEFAULT_READ_NAME_REGEX;

    @Argument(doc = "A fast, regex-free read name format. ILLUMINA_UMI accepts the default 5 or 7 field Illumina names, " +
            "and also names with one extra trailing UMI field made only of A, C, G, T, N, + and -, e.g. " +
            "FC:1:1101:12345:67890:ACGT-TTGC. Cannot be combined with a READ_NAME_REGEX other than the default.", optional = true)
    public ReadNameFormat READ_NAME_FORMAT = null;

    @Argument(doc = "Use the exact two-pass mode: a first pass finds each tile's bounds, then a second pass decides the reads. " +
            "By default the approximate single-pass mode runs instead, which reads INPUT once and might differ slightly from " +
            "two-pass. Single-pass needs READ_COUNT; without it, or for CRAM input, two-pass runs anyway.", optional = true)
    public boolean TWO_PASS = false;

    @Argument(doc = "The number of reads in INPUT (records, both mates counted), for single-pass. Each tile's count is scaled " +
            "to the whole input with it when widening. Fails if the number of reads seen differs by more than 1%. Required for " +
            "single-pass input from a pipe.", optional = true)
    public Long READ_COUNT = null;

    @Argument(doc = "Stop after processing N reads, mainly for debugging.", optional = true)
    public Long STOP_AFTER = null;

    @Argument(doc = "Allow downsampling again despite this being a bad idea with possibly unexpected results.", optional = true)
    public boolean ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS = false;

    /** Read name formats that {@link #READ_NAME_FORMAT} can select. */
    public enum ReadNameFormat {
        /** Illumina names with an optional trailing UMI field. See {@link #READ_NAME_FORMAT}. */
        ILLUMINA_UMI
    }

    private final Log log = Log.getInstance(PositionBasedDownsampleSam.class);

    private ReadNameParser readNameParser;
    private long total = 0;
    private long kept = 0;
    public static String PG_PROGRAM_NAME = "PositionBasedDownsampleSam";
    private final static double ACCEPTABLE_FUDGE_FACTOR = 0.2;

    /* max-position in tile as a function of tile. We might need to
       look per-readgroup, but at this point I'm making the assumptions that I need to downsample a
       sample where all the readgroups came from the same type of flowcell. */

    CollectionUtil.DefaultingMap.Factory<Coord, Short> defaultingMapFactory = aShort -> new Coord();

    final private Map<Short, Coord> tileCoord = new CollectionUtil.DefaultingMap<>(defaultingMapFactory, true);

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errors = new ArrayList<>();

        if (FRACTION < 0 || FRACTION > 1) {
            errors.add("FRACTION must be a value between 0 and 1, found: " + FRACTION);
        }

        if (READ_NAME_FORMAT != null && !ReadNameParser.DEFAULT_READ_NAME_REGEX.equals(READ_NAME_REGEX)) {
            errors.add("READ_NAME_FORMAT and READ_NAME_REGEX cannot both be set.");
        }

        if (READ_COUNT != null && READ_COUNT < 1) {
            errors.add("READ_COUNT must be at least 1, found: " + READ_COUNT);
        }

        if (!errors.isEmpty())
            return errors.toArray(new String[errors.size()]);

        return super.customCommandLineValidation();
    }

    @Override
    protected int doWork() {

        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        readNameParser = new ReadNameParser(READ_NAME_FORMAT == ReadNameFormat.ILLUMINA_UMI
                ? ReadNameParser.ILLUMINA_UMI_READ_NAME_REGEX
                : READ_NAME_REGEX);

        final boolean isCram = INPUT.getName().toLowerCase().endsWith(".cram");
        if (!TWO_PASS && READ_COUNT != null && !isCram) {
            singlePassSamRecords();
            logFinalStats();
            return 0;
        }
        if (!INPUT.isFile()) {
            throw new PicardException("Two-pass mode needs INPUT to be a regular file, not a pipe. " +
                    "Use single-pass with READ_COUNT.");
        }
        if (!TWO_PASS) {
            log.info("Using two-pass mode: " + (isCram ? "single-pass cannot read CRAM." : "READ_COUNT is not given."));
        }

        log.info("Checking to see if input file has been downsampled with this program before.");
        final SamReader header = openReader();
        checkProgramRecords(header.getFileHeader());
        CloserUtil.close(header);

        log.info("Starting first pass. Examining read distribution in tiles.");
        fillTileMinMaxCoord();
        log.info("First pass done.");

        log.info("Starting second pass. Outputting reads.");
        outputSamRecords();
        log.info("Second pass done.");

        logFinalStats();
        return 0;
    }

    private void logFinalStats() {
        final double finalP = kept / (double) total;
        if (Math.abs(finalP - FRACTION) / (Math.min(finalP, FRACTION) + 1e-10) > ACCEPTABLE_FUDGE_FACTOR) {
            log.warn(String.format("You've requested FRACTION=%g, the resulting downsampling resulted in a rate of %f.", FRACTION, finalP));
        }
        log.info(String.format("Finished! Kept %d out of %d reads (P=%g).", kept, total, finalP));
    }

    /** The input header with this program's PG record added, for the output. */
    private SAMFileHeader outputHeader(final SAMFileHeader inputHeader) {
        final SAMFileHeader header = inputHeader.clone();
        final SAMFileHeader.PgIdGenerator pgIdGenerator = new SAMFileHeader.PgIdGenerator(header);
        final SAMProgramRecord programRecord = new SAMProgramRecord(pgIdGenerator.getNonCollidingId(PG_PROGRAM_NAME));
        programRecord.setProgramName(PG_PROGRAM_NAME);
        programRecord.setCommandLine(getCommandLine());
        programRecord.setProgramVersion(getVersion());
        header.addProgramRecord(programRecord);
        return header;
    }

    /** Async writing: BGZF compression and record handling overlap. The queue fits single-pass bursts. */
    private SAMFileWriter openOutput(final SAMFileHeader inputHeader) {
        return new SAMFileWriterFactory()
                .setUseAsyncIo(true)
                .setAsyncOutputBufferSize(200_000)
                .makeWriter(outputHeader(inputHeader), true, OUTPUT, REFERENCE_SEQUENCE);
    }

    private void outputSamRecords() {

        final ProgressLogger progress = new ProgressLogger(log, (int) 1e7);
        final SamReader in = openReader();
        final SAMFileWriter out = openOutput(in.getFileHeader());

        final CircleSelector selector = new CircleSelector(FRACTION);

        for (final SAMRecord rec : in) {
            if (STOP_AFTER != null && total >= STOP_AFTER) break;

            total++;

            final PhysicalLocationInt pos = getSamRecordLocation(rec);
            final boolean keepRecord = selector.select(pos, tileCoord.get(pos.getTile()));

            if (keepRecord) {
                if (REMOVE_DUPLICATE_INFORMATION) rec.setDuplicateReadFlag(false);
                out.addAlignment(rec);
                kept++;
            }
            progress.record(rec);
        }

        out.close();

        CloserUtil.close(in);
    }

    private void checkProgramRecords(final SAMFileHeader header) {
        for (final SAMProgramRecord pg : header.getProgramRecords()) {
            if (pg.getProgramName() != null && pg.getProgramName().equals(PG_PROGRAM_NAME)) {

                final String outText = "Found previous Program Record that indicates that this file has been downsampled already with this program. Operation not supported! Previous PG: " + pg.toString();

                if (ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS) {
                    log.warn(outText);
                } else {
                    log.error(outText);
                    throw new PicardException(outText);
                }
            }
        }
    }

    // scan all the tiles and find the smallest and largest coordinate (x & y) in that tile. Leaves the bounds unwidened.
    private void fillTileMinMaxCoord() {

        final SamReader in = openReader();

        final ProgressLogger progress = new ProgressLogger(log, (int) 1e7, "Read");

        long total = 0;

        for (final SAMRecord rec : in) {
            if (STOP_AFTER != null && total >= STOP_AFTER) break;

            total++;
            progress.record(rec);
            final PhysicalLocationInt location = getSamRecordLocation(rec);

            //Defaulting map will create a new Coord if it's not there.

            final Coord Pos = tileCoord.get(location.getTile());

            Pos.maxX = Math.max(Pos.maxX, location.getX());
            Pos.minX = Math.min(Pos.minX, location.getX());

            Pos.maxY = Math.max(Pos.maxY, location.getY());
            Pos.minY = Math.min(Pos.minY, location.getY());

            Pos.count++;

        }

        // now that we know what the maximal/minimal numbers were, we should increase/decrease them a little, to account for sampling error
        for (final Coord coord : tileCoord.values()) {

            final int diffX = coord.maxX - coord.minX;
            final int diffY = coord.maxY - coord.minY;

            coord.maxX += diffX / coord.count;
            coord.minX -= diffX / coord.count;

            coord.maxY += diffY / coord.count;
            coord.minY -= diffY / coord.count;
        }

        CloserUtil.close(in);
    }

    // Single-pass tuning: a tile's box is stable once no extreme has moved for max(MIN_READS, FRACTION x reads) reads.
    private static final long SINGLE_PASS_MIN_READS = 20000;
    private static final double SINGLE_PASS_STABLE_FRACTION = 0.1;

    /** Per-tile state for single-pass: the box learned so far, the reads held until it freezes, and the frozen box. */
    private class TilePool {
        final Coord raw = new Coord();   // unwidened, same starting values as pass 1
        long lastChange = 0;             // raw.count at which an extreme last moved
        final List<SAMRecord> buffer = new ArrayList<>();
        final List<PhysicalLocationInt> bufferedPos = new ArrayList<>(); // parsed once, on arrival
        Coord frozen = null;

        void observe(final PhysicalLocationInt pos) {
            boolean changed = false;
            if (pos.getX() > raw.maxX) { raw.maxX = pos.getX(); changed = true; }
            if (pos.getX() < raw.minX) { raw.minX = pos.getX(); changed = true; }
            if (pos.getY() > raw.maxY) { raw.maxY = pos.getY(); changed = true; }
            if (pos.getY() < raw.minY) { raw.minY = pos.getY(); changed = true; }
            raw.count++;
            if (changed) {
                lastChange = raw.count;
            }
        }

        boolean isStable() {
            final double window = Math.max(SINGLE_PASS_MIN_READS, SINGLE_PASS_STABLE_FRACTION * raw.count);
            return raw.count >= SINGLE_PASS_MIN_READS && raw.count - lastChange >= window;
        }

        /** Widens the raw box as pass 1 does, with the count scaled up to the whole input by countScale. */
        void freeze(final double countScale) {
            final long count = Math.max(raw.count, Math.round(raw.count * countScale));
            final int diffX = raw.maxX - raw.minX;
            final int diffY = raw.maxY - raw.minY;
            final Coord w = new Coord();
            w.count = count;
            w.minX = raw.minX - (int) (diffX / count);
            w.maxX = raw.maxX + (int) (diffX / count);
            w.minY = raw.minY - (int) (diffY / count);
            w.maxY = raw.maxY + (int) (diffY / count);
            frozen = w;
        }
    }

    /**
     * Single-pass: reads INPUT once. Each tile's reads are held until its box is stable, then decided against the widened box;
     * later reads of the tile are decided on arrival. INPUT is opened once, so it may be a pipe.
     */
    private void singlePassSamRecords() {
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e7);
        final SamReader in;
        try {
            // htsjdk sniffs the format from the raw bytes, so pass the stream unwrapped.
            in = SamReaderFactory.makeDefault()
                    .referenceSequence(REFERENCE_SEQUENCE)
                    .setUseAsyncIo(true)
                    .open(SamInputResource.of(new FileInputStream(INPUT)));
        } catch (final IOException e) {
            throw new PicardException("Could not open " + INPUT + ": " + e.getMessage(), e);
        }
        checkProgramRecords(in.getFileHeader());

        final SAMFileWriter out = openOutput(in.getFileHeader());
        final CircleSelector selector = new CircleSelector(FRACTION);
        final Map<Short, TilePool> pools = new HashMap<>();

        for (final SAMRecord rec : in) {
            if (STOP_AFTER != null && total >= STOP_AFTER) break;
            total++;

            final PhysicalLocationInt pos = getSamRecordLocation(rec);
            final TilePool pool = pools.computeIfAbsent(pos.getTile(), tile -> new TilePool());
            if (pool.frozen != null) {
                emitIfKept(rec, pos, pool.frozen, selector, out);
            } else {
                pool.observe(pos);
                pool.buffer.add(rec);
                pool.bufferedPos.add(pos);
                if (pool.isStable()) {
                    freezeAndFlush(pool, selector, out);
                }
            }
            progress.record(rec);
        }

        // Pools that never became stable are frozen with what they have.
        for (final TilePool pool : pools.values()) {
            if (pool.frozen == null) {
                freezeAndFlush(pool, selector, out);
            }
        }

        out.close();
        CloserUtil.close(in);

        final double readCountError = Math.abs(total - READ_COUNT) / (double) READ_COUNT;
        if (readCountError > 0.01) {
            throw new PicardException(String.format("READ_COUNT=%d but %d reads were read. The bounds were scaled with READ_COUNT, " +
                    "so the output is not usable.", READ_COUNT, total));
        }
        if (total != READ_COUNT) {
            log.warn(String.format("READ_COUNT=%d but %d reads were read.", READ_COUNT, total));
        }
    }

    private void freezeAndFlush(final TilePool pool, final CircleSelector selector, final SAMFileWriter out) {
        // Scale the tile's count by reads of the whole input per read seen so far.
        pool.freeze(READ_COUNT / (double) Math.max(1L, total));
        for (int i = 0; i < pool.buffer.size(); i++) {
            emitIfKept(pool.buffer.get(i), pool.bufferedPos.get(i), pool.frozen, selector, out);
        }
        pool.buffer.clear();
        pool.bufferedPos.clear();
    }

    private void emitIfKept(final SAMRecord rec, final PhysicalLocationInt pos, final Coord box,
                            final CircleSelector selector, final SAMFileWriter out) {
        if (selector.select(pos, box)) {
            if (REMOVE_DUPLICATE_INFORMATION) rec.setDuplicateReadFlag(false);
            out.addAlignment(rec);
            kept++;
        }
    }

    private SamReader openReader() {
        return SamReaderFactory.makeDefault()
                .referenceSequence(REFERENCE_SEQUENCE)
                .setUseAsyncIo(true) // background-thread BAM reading (htsjdk); output is unchanged
                .open(INPUT);
    }

    private PhysicalLocationInt getSamRecordLocation(final SAMRecord rec) {
        final PhysicalLocationInt pos = new PhysicalLocationInt();
        readNameParser.addLocationInformation(rec.getReadName(), pos);
        return pos;
    }

    /*
     * The reads are selected depending on whether they are in a periodically repeating circle whose representative
     * overlaps the boundary of the tile. The circle is chosen as to have an area of FRACTION and the also a overlap
     * of FRACTION with both the bottom and left edges of the unit square ([0,1] x [0,1]), which defines it uniquely.
     * Finally the repeating pattern is there to make sure that the mask on the flowcell also has minimal boundary.
     *
     * The position of the reads is mapped into the unit square using the min/max coordinates of the tile prior to
     * masking
     *
     * This choice of a mask is intended to accomplish several goasl:
     *  - pick out a fraction FRACTION of the reads
     *  - pick out a fraction FRACTION of the reads that are near the boundaries
     *  - keep nearby reads (in a tile) together by minimizing the boundary of the mask itself
     *  - keep nearby reads (in neighboring tiles) together (since they might be optical duplicates) by keeping that the
     *  mask is the same on all tiles, and by having the same mask on the left edge as on the right (same for top/bottom).
     */
    private class CircleSelector {

        private final double radiusSquared;
        private final double offset;
        private final boolean positiveSelection;

        CircleSelector(final double fraction) {

            final double p;
            if (fraction > 0.5) {
                p = 1 - fraction;
                positiveSelection = false;
            } else {
                p = fraction;
                positiveSelection = true;
            }
            radiusSquared = p / Math.PI; //thus the area is \pi r^2 = p, thus a fraction p of the unit square will be captured

            /* if offset is used as the center of the circle (both x and y), this makes the overlap
             region with each of the boundaries of the unit square have length p (and thus a fraction
             p of the boundaries of each tile will be removed) */

            if (p < 0) {
                // at this point a negative p will result in a square-root of a negative number in the next step.
                throw new PicardException("This shouldn't happen...");
            }
            offset = Math.sqrt(radiusSquared - p * p / 4);
        }

        private double roundedPart(final double x) {return x - Math.round(x);}

        // this function checks to see if the location of the read is within the masking circle
        private boolean select(final PhysicalLocationInt coord, final Coord tileCoord) {
            // r^2 = (x-x_0)^2 + (y-y_0)^2, where both x_0 and y_0 equal offset
            final double distanceSquared =
                            Math.pow(roundedPart(((coord.getX() - tileCoord.minX) / (double) (tileCoord.maxX - tileCoord.minX)) - offset), 2) +
                            Math.pow(roundedPart(((coord.getY() - tileCoord.minY) / (double) (tileCoord.maxY - tileCoord.minY)) - offset), 2);

            return (distanceSquared > radiusSquared) ^ positiveSelection;
        }
    }

    private class Coord {
        public int minX;
        public int minY;
        public int maxX;
        public int maxY;
        public long count;

        public Coord() {
            count = 0;
            minX = 0;
            minY = 0;
            maxX = 0;
            maxY = 0;
        }
    }
}
