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
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.IntervalsManipulationProgramGroup;

import java.io.File;
import java.util.List;

/**
 * This tool adjusts the coordinates in an interval list on one reference to its homologous interval list on another
 * reference, based on a chain file that describes the correspondence between the two references. It is based on the
 * <a href="http://genome.ucsc.edu/cgi-bin/hgLiftOver">UCSC LiftOver tool</a> and uses a UCSC chain file to guide its operation.
 * It accepts both Picard interval_list files or VCF files as interval inputs.
 * <br />
 * <h3>Usage example:</h3>
 * <pre>
 * java -jar picard.jar LiftOverIntervalList \
 *       I=input.interval_list \
 *       O=output.interval_list \
 *       SD=reference_sequence.dict \
 *       CHAIN=build.chain
 * </pre>
 * <p>
 * <h3>Return codes</h3>
 * If all the intervals lifted over successfully, program will return 0. It will return 1 otherwise.
 * <p>
 * <h3>Caveats</h3>
 * An interval is "lifted" in its entirety, but it might intersect (a "hit") with multiple chain-blocks.
 * Instead of placing the interval in multiple hits, it is lifted over using the first hit that passes the
 * threshold of {@link #MIN_LIFTOVER_PCT}. For large enough {@link #MIN_LIFTOVER_PCT} this is non-ambiguous,
 * but if one uses small values of {@link #MIN_LIFTOVER_PCT} (perhaps in order to increase the rate of successful
 * hits...) the liftover could end up going to the smaller of two good hits. On the other hand, if none of the hits
 * pass the threshold a warning will be emitted and the interval will not be lifted.
 *
 * @author alecw@broadinstitute.org
 */
@CommandLineProgramProperties(
        summary = LiftOverIntervalList.USAGE_SUMMARY + LiftOverIntervalList.USAGE_DETAILS,
        oneLineSummary = LiftOverIntervalList.USAGE_SUMMARY,
        programGroup = IntervalsManipulationProgramGroup.class
)
@DocumentedFeature
public class LiftOverIntervalList extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Lifts over an interval list from one reference build to another. ";
    static final String USAGE_DETAILS = "This tool adjusts the coordinates in an interval list on one reference to its homologous " +
            "interval list on another " +
            "reference, based on a chain file that describes the correspondence between the two references. It is based on the " +
            "UCSC LiftOver tool (see: http://genome.ucsc.edu/cgi-bin/hgLiftOver) and uses a UCSC chain file to guide its operation. " +
            "It accepts both Picard interval_list files or VCF files as interval inputs.\n" +
            "\n" +
            "<h3>Usage example:</h3>" +
            "java -jar picard.jar LiftOverIntervalList \\\n" +
            "      I=input.interval_list \\\n" +
            "      O=output.interval_list \\\n" +
            "      SD=reference_sequence.dict \\\n" +
            "      CHAIN=build.chain" +
            "</pre>" +
            "\n" +
            "<h3>Return codes</h3>\n" +
            "If all the intervals lifted over successfully, program will return 0. It will return 1 otherwise.\n" +
            "\n" +
            "<h3>Caveats</h3>\n" +
            "An interval is \"lifted\" in its entirety, but it might intersect (a \"hit\") with multiple chain-blocks. " +
            "Instead of placing the interval in multiple hits, it is lifted over using the first hit that passes the " +
            "threshold of MIN_LIFTOVER_PCT. For large enough MIN_LIFTOVER_PCT this is non-ambiguous, but if one uses small values of MIN_LIFTOVER_PCT " +
            "(perhaps in order to increase the rate of successful hits...) the liftover could end up going to the smaller of two " +
            "good hits. On the other hand, if none of the hits pass the threshold a warning will be emitted and the interval will " +
            "not be lifted.";
    private static final Log LOG = Log.getInstance(LiftOverIntervalList.class);

    @Argument(doc = "The input interval list to be lifted over.", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "The output interval list file.", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Argument(doc = "Sequence dictionary to place in the output interval list. (This should be the dictionary of the target reference.)",
            shortName = StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME)
    public File SEQUENCE_DICTIONARY;

    @Argument(doc = "Chain file that guides the LiftOver process.")
    public File CHAIN;

    @Argument(doc = "Minimum percentage of bases in each input interval that must map to output interval for liftover of " +
            "that interval to occur. If the program fails to find a good target for an interval, a warning will be emitted " +
            "and the interval will be dropped from the output. ")
    public double MIN_LIFTOVER_PCT = LiftOver.DEFAULT_LIFTOVER_MINMATCH;

    @Argument(doc = "Interval List file for intervals that were rejected", optional = true)
    public File REJECT = null;

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
        if (REJECT != null) IOUtil.assertFileIsWritable(REJECT);

        final LiftOver liftOver = new LiftOver(CHAIN);
        liftOver.setLiftOverMinMatch(MIN_LIFTOVER_PCT);

        final IntervalList intervalList = IntervalList.fromFile(INPUT);
        final IntervalList rejects = new IntervalList(intervalList.getHeader());

        final long baseCount = intervalList.getBaseCount();
        LOG.info("Lifting over " + intervalList.getIntervals().size() + " intervals, encompassing " +
                baseCount + " bases.");

        final SAMFileHeader toHeader = SamReaderFactory.makeDefault().getFileHeader(SEQUENCE_DICTIONARY);
        liftOver.validateToSequences(toHeader.getSequenceDictionary());
        final IntervalList toIntervals = new IntervalList(toHeader);
        for (final Interval fromInterval : intervalList) {
            final Interval toInterval = liftOver.liftOver(fromInterval);
            if (toInterval != null) {
                toIntervals.add(toInterval);
            } else {
                rejects.add(fromInterval);
                LOG.warn("Liftover failed for ", fromInterval, " (len ", fromInterval.length(), ")");
                final List<LiftOver.PartialLiftover> partials = liftOver.diagnosticLiftover(fromInterval);
                for (final LiftOver.PartialLiftover partial : partials) {
                    LOG.info(partial);
                }
            }
        }

        toIntervals.sorted().write(OUTPUT);

        if (REJECT != null) {
            rejects.write(REJECT);
        }
        final long rejectBaseCount = rejects.getBaseCount();

        LOG.info(String.format("Liftover Complete. \n" +
                        "%d of %d intervals failed (%g%%) to liftover, encompassing %d of %d bases (%g%%).",
                rejects.getIntervals().size(), intervalList.getIntervals().size(),
                100 * rejects.getIntervals().size() / (double) intervalList.getIntervals().size(),
                rejectBaseCount, baseCount,
                100 * rejectBaseCount / (double) baseCount
        ));

        return rejects.getIntervals().isEmpty() ? 0 : 1;
    }
}
