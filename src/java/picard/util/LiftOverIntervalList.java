/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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
package picard.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Intervals;

import java.io.File;
import java.util.List;

/**
 * @author alecw@broadinstitute.org
 */
@CommandLineProgramProperties(
        usage = "Lifts over an interval list from one reference build to another. Based on UCSC liftOver." +
                " Uses a UCSC chain file to guide the liftOver.",
        usageShort = "Lifts over an interval list from one reference build to another",
        programGroup = Intervals.class
)
public class LiftOverIntervalList extends CommandLineProgram {

    private static final Log LOG = Log.getInstance(LiftOverIntervalList.class);

    @Option(doc = "Interval list to be lifted over.", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Option(doc = "Where to write lifted-over interval list.", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Option(doc = "Sequence dictionary to write into the output interval list.",
            shortName = StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME)
    public File SEQUENCE_DICTIONARY;

    @Option(doc = "Chain file that guides LiftOver.")
    public File CHAIN;

    @Option(doc = "Minimum percentage of bases in each input interval that must map to output interval.")
    public double MIN_LIFTOVER_PCT = LiftOver.DEFAULT_LIFTOVER_MINMATCH;

    public static void main(final String[] argv) {
        new LiftOverIntervalList().instanceMainWithExit(argv);
    }

    /**
     * Do the work after command line has been parsed. RuntimeException may be
     * thrown by this method, and are reported appropriately.
     *
     * @return program exit status.
     */
    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(SEQUENCE_DICTIONARY);
        IOUtil.assertFileIsReadable(CHAIN);
        IOUtil.assertFileIsWritable(OUTPUT);

        final LiftOver liftOver = new LiftOver(CHAIN);
        liftOver.setLiftOverMinMatch(MIN_LIFTOVER_PCT);

        final IntervalList fromIntervals = IntervalList.fromFile(INPUT);
        final SAMFileHeader toHeader = SamReaderFactory.makeDefault().open(SEQUENCE_DICTIONARY).getFileHeader();
        liftOver.validateToSequences(toHeader.getSequenceDictionary());
        final IntervalList toIntervals = new IntervalList(toHeader);
        boolean anyFailed = false;
        for (final Interval fromInterval : fromIntervals) {
            final Interval toInterval = liftOver.liftOver(fromInterval);
            if (toInterval != null) {
                toIntervals.add(toInterval);
            } else {
                anyFailed = true;
                LOG.warn("Liftover failed for ", fromInterval, "(len ", fromInterval.length(), ")");
                final List<LiftOver.PartialLiftover> partials = liftOver.diagnosticLiftover(fromInterval);
                for (final LiftOver.PartialLiftover partial : partials) {
                    LOG.info(partial);
                }
            }
        }

        toIntervals.sorted();
        toIntervals.write(OUTPUT);
        return anyFailed ? 1 : 0;
    }
}
