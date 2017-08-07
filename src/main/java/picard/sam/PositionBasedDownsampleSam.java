/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;
import picard.sam.util.PhysicalLocationInt;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Class to downsample a BAM file while respecting that we should either get rid
 * of both ends of a pair or neither end of the pair. In addition, this program uses the read-name
 * and extracts the position within the tile whence the read came from. The downsampling is based on this position.
 * <p/>
 * Note 1: This is technology and read-name dependent. If your read-names do not have coordinate information, or if your
 * BAM contains reads from multiple technologies (flowcell versions, sequencing machines) this will not work properly.
 * This has been designed with Illumina MiSeq/HiSeq in mind.
 * <p/>
 * Note 2: The downsampling is _not_ random. It is deterministically dependent on the position of the read within its tile. Specifically,
 * it draws out an ellipse that covers a FRACTION fraction of the area and each of the edges and uses this to determine whether to keep the
 * record. Since reads with the same name have the same position (mates, secondary and supplemental alignments), the decision will be the same for all of them.
 * <p/>
 * Finally, the code has been designed to simulate sequencing less as accurately as possible, not for getting an exact downsample fraction.
 * In particular, since the reads may be distributed non-evenly within the lanes/tiles, the resulting downsampling percentage will not be accurately
 * determined by the input argument FRACTION. One should re-MarkDuplicates after downsampling in order to "expose" the duplicates whose representative has
 * been downsampled away.
 *
 * @author Yossi Farjoun
 */
@CommandLineProgramProperties(
        summary = "Class to downsample a BAM file while respecting that we should either get rid of both ends of a pair or neither \n" +
                "end of the pair. In addition, this program uses the read-name and extracts the position within the tile whence \n" +
                "the read came from. The downsampling is based on this position. Results with the exact same input will produce the \n" +
                "same results.\n" +
                "\n" +
                "Note 1: This is technology and read-name dependent. If your read-names do not have coordinate information, or if your\n" +
                "BAM contains reads from multiple technologies (flowcell versions, sequencing machines) this will not work properly. \n" +
                "This has been designed with Illumina MiSeq/HiSeq in mind.\n" +
                "Note 2: The downsampling is not random. It is deterministically dependent on the position of the read within its tile.\n" +
                "Note 3: Downsampling twice with this program is not supported.\n" +
                "Note 4: You should call MarkDuplicates after downsampling.\n" +
                "\n" +
                "Finally, the code has been designed to simulate sequencing less as accurately as possible, not for getting an exact downsample \n" +
                "fraction. In particular, since the reads may be distributed non-evenly within the lanes/tiles, the resulting downsampling \n" +
                "percentage will not be accurately determined by the input argument FRACTION.",
        oneLineSummary = "Downsample a SAM or BAM file to retain a subset of the reads based on the reads location in each tile in the flowcell.",
        programGroup = SamOrBam.class
)
public class PositionBasedDownsampleSam extends CommandLineProgram {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to downsample.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output, downsampled, SAM or BAM file to write.")
    public File OUTPUT;

    @Argument(shortName = "F", doc = "The (approximate) fraction of reads to be kept, between 0 and 1.", optional = false)
    public Double FRACTION = null;

    @Argument(doc = "Stop after processing N reads, mainly for debugging.", optional = true)
    public Long STOP_AFTER = null;

    @Argument(doc = "Allow Downsampling again despite this being a bad idea with possibly unexpected results.", optional = true)
    public boolean ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS = false;

    @Argument(doc = "Determines whether the duplicate tag should be reset since the downsampling requires re-marking duplicates.")
    public boolean REMOVE_DUPLICATE_INFORMATION = true;

    private final Log log = Log.getInstance(PositionBasedDownsampleSam.class);

    private OpticalDuplicateFinder opticalDuplicateFinder;
    private long total = 0;
    private long kept = 0;
    public static String PG_PROGRAM_NAME = "PositionBasedDownsampleSam";
    private final static double ACCEPTABLE_FUDGE_FACTOR = 0.2;

    /* max-position in tile as a function of tile. We might need to
       look per-readgroup, but at this point I'm making the assumptions that I need to downsample a
       sample where all the readgroups came from the same type of flowcell. */

    CollectionUtil.DefaultingMap.Factory<Coord, Short> defaultingMapFactory = new CollectionUtil.DefaultingMap.Factory<Coord, Short>() {
        @Override
        public Coord make(final Short aShort) {
            return new Coord();
        }
    };

    final private Map<Short, Coord> tileCoord = new CollectionUtil.DefaultingMap<Short, Coord>(defaultingMapFactory, true);


    final Map<Short, Histogram<Short>> xPositions = new HashMap<Short, Histogram<Short>>();
    final Map<Short, Histogram<Short>> yPositions = new HashMap<Short, Histogram<Short>>();

    public static void main(final String[] args) {
        new PositionBasedDownsampleSam().instanceMainWithExit(args);
    }

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errors = new ArrayList<String>();

        if (FRACTION < 0 || FRACTION > 1) {
            errors.add("FRACTION must be a value between 0 and 1, found: " + FRACTION);
        }

        if (errors.isEmpty()) {
            return null;
        } else {
            return errors.toArray(new String[errors.size()]);
        }
    }

    @Override
    protected int doWork() {

        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        log.info("Checking to see if input file has been downsampled with this program before.");
        checkProgramRecords();

        opticalDuplicateFinder = new OpticalDuplicateFinder();

        log.info("Starting first pass. Examining read distribution in tiles.");
        fillTileMinMaxCoord();
        log.info("First pass done.");

        log.info("Starting second pass. Outputting reads.");
        outputSamRecords();
        log.info("Second pass done.");

        final double finalP = kept / (double) total;
        if (Math.abs(finalP - FRACTION) / (Math.min(finalP, FRACTION) + 1e-10) > ACCEPTABLE_FUDGE_FACTOR) {
            log.warn(String.format("You've requested FRACTION=%g, the resulting downsampling resulted in a rate of %f.", FRACTION, finalP));
        }
        log.info(String.format("Finished! Kept %d out of %d reads (P=%g).", kept, total, finalP));

        return 0;
    }

    private void outputSamRecords() {

        final ProgressLogger progress = new ProgressLogger(log, (int) 1e7);
        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);

        final SAMFileHeader header = in.getFileHeader().clone();
        final SAMFileHeader.PgIdGenerator pgIdGenerator = new SAMFileHeader.PgIdGenerator(header);
        final SAMProgramRecord programRecord = new SAMProgramRecord(pgIdGenerator.getNonCollidingId(PG_PROGRAM_NAME));

        programRecord.setProgramName(PG_PROGRAM_NAME);
        programRecord.setCommandLine(getCommandLine());
        programRecord.setProgramVersion(getVersion());
        header.addProgramRecord(programRecord);

        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);

        final CircleSelector selector = new CircleSelector(FRACTION);

        for (final SAMRecord rec : in) {
            if (STOP_AFTER != null && total >= STOP_AFTER) break;

            total++;

            final PhysicalLocationInt pos = getSamRecordLocation(rec);

            if (!xPositions.containsKey(pos.getTile())) {
                xPositions.put(pos.getTile(), new Histogram<Short>(pos.getTile() + "-xpos", "count"));
            }
            if (!yPositions.containsKey(pos.getTile())) {
                yPositions.put(pos.getTile(), new Histogram<Short>(pos.getTile() + "-ypos", "count"));
            }

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

    private void checkProgramRecords() {
        final SamReader in = SamReaderFactory
                .makeDefault()
                .referenceSequence(REFERENCE_SEQUENCE)
                .open(INPUT);

        for (final SAMProgramRecord pg : in.getFileHeader().getProgramRecords()) {
            if (pg.getProgramName() != null && pg.getProgramName().equals(PG_PROGRAM_NAME)) {

                final String outText = "Found previous Program Record that indicates that this BAM has been downsampled already with this program. Operation not supported! Previous PG: " + pg.toString();

                if (ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS) {
                    log.warn(outText);
                } else {
                    log.error(outText);
                    throw new PicardException(outText);
                }
            }
        }
        CloserUtil.close(in);
    }

    // scan all the tiles and find the smallest and largest coordinate (x & y) in that tile.
    private void fillTileMinMaxCoord() {

        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);

        final ProgressLogger progress = new ProgressLogger(log, (int) 1e7, "Read");

        int total = 0;

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

    private PhysicalLocationInt getSamRecordLocation(final SAMRecord rec) {
        final PhysicalLocationInt pos = new PhysicalLocationInt();
        opticalDuplicateFinder.addLocationInformation(rec.getReadName(), pos);
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
        public int count;

        public Coord() {
            count = 0;
            minX = 0;
            minY = 0;
            maxX = 0;
            maxY = 0;
        }
    }
}
