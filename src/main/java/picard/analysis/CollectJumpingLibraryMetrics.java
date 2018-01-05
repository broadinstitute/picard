/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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

package picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import picard.sam.DuplicationMetrics;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Command-line program to compute metrics about outward-facing pairs, inward-facing
 * pairs, and chimeras in a jumping library.
 *
 * @author ktibbett@broadinstitute.org
 */
@CommandLineProgramProperties(
        summary = CollectJumpingLibraryMetrics.USAGE_SUMMARY + CollectJumpingLibraryMetrics.USAGE_DETAILS,
        oneLineSummary = CollectJumpingLibraryMetrics.USAGE_SUMMARY,
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@DocumentedFeature
public class CollectJumpingLibraryMetrics extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Collect jumping library metrics. ";
    static final String USAGE_DETAILS = "<p>This tool collects high-level metrics about the " +
"presence of outward-facing (jumping) and inward-facing (non-jumping) read pairs within a SAM or BAM file." +
"For a brief primer on jumping libraries, see the GATK "+
"<a href='https://www.broadinstitute.org/gatk/guide/article?id=6326'>Dictionary</a></p>." +

"<p>This program gets all data for computation from the first read in each pair in which the mapping quality (MQ) tag " +
"is set with the mate's mapping quality.  If the MQ tag is not set, then the program assumes that the mate's MQ is " +
"greater than or equal to MINIMUM_MAPPING_QUALITY (default value is 0).</p> "+

"<p>Note: Metrics labeled as percentages are actually expressed as fractions!</p>" +

"<h4>Usage example:</h4>" +
"<pre>" +
"java -jar picard.jar CollectJumpingLibraryMetrics \\<br />" +
"      I=input.bam  \\<br />" +
"      O=jumping_metrics.txt" +
"</pre>" +

"Please see the output metrics documentation on "+
"<a href='https://broadinstitute.github.io/picard/picard-metric-definitions.html#JumpingLibraryMetrics'>JumpingLibraryMetrics</a> "+
"for detailed explanations of the output metrics."+
"<hr />";

    // Usage and parameters

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "BAM file(s) of reads with duplicates marked")
    public List<File> INPUT = new ArrayList<File>();
    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "File to which metrics should be written")
    public File OUTPUT;
    @Argument(shortName = StandardOptionDefinitions.MINIMUM_MAPPING_QUALITY_SHORT_NAME, doc = "Mapping quality minimum cutoff")
    public Integer MINIMUM_MAPPING_QUALITY = 0;
    @Argument(shortName = "T", doc = "When calculating mean and stdev stop when the bins in the tail of the distribution " +
            "contain fewer than mode/TAIL_LIMIT items")
    public int TAIL_LIMIT = 10000;
    @Argument(doc = "Jumps greater than or equal to the greater of this value or 2 times the mode of the " +
            "outward-facing pairs are considered chimeras")
    public int CHIMERA_KB_MIN = 100000;

    private static final int SAMPLE_FOR_MODE = 50000; // How many outward-facing pairs to sample to determine the mode

    /** Stock main method. */
    public static void main(String[] args) {
        System.exit(new CollectJumpingLibraryMetrics().instanceMain(args));
    }

    /**
     * Calculates the detailed statistics about the jumping library and then generates the results.
     */
    protected int doWork() {

        for (File f : INPUT) {
            IOUtil.assertFileIsReadable(f);
        }
        IOUtil.assertFileIsWritable(OUTPUT);

        Histogram<Integer> innieHistogram = new Histogram<Integer>();
        Histogram<Integer> outieHistogram = new Histogram<Integer>();

        int fragments = 0;
        int innies = 0;
        int outies = 0;
        int innieDupes = 0;
        int outieDupes = 0;
        int crossChromPairs = 0;
        int superSized = 0;
        int tandemPairs = 0;
        double chimeraSizeMinimum = Math.max(getOutieMode(), (double) CHIMERA_KB_MIN);

        for (File f : INPUT) {
            SamReader reader = SamReaderFactory.makeDefault().open(f);

            if (reader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
                throw new PicardException("SAM file must " + f.getName() + " must be sorted in coordintate order");
            }

            for (SAMRecord sam : reader) {

                // We're getting all our info from the first of each pair.
                if (!sam.getFirstOfPairFlag()) {
                    continue;
                }

                // Ignore unmapped read pairs
                if (sam.getReadUnmappedFlag()) {
                    if (!sam.getMateUnmappedFlag()) {
                        fragments++;
                        continue;
                    }

                    // If both ends are unmapped and we've hit unaligned reads we're done
                    if (sam.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
                        break;
                    }
                    continue;
                }

                if (sam.getMateUnmappedFlag()) {
                    fragments++;
                    continue;
                }

                // Ignore low-quality reads.  If we don't have the mate mapping quality, assume it's OK
                if ((sam.getAttribute(SAMTag.MQ.name()) != null &&
                        sam.getIntegerAttribute(SAMTag.MQ.name()) < MINIMUM_MAPPING_QUALITY) ||
                        sam.getMappingQuality() < MINIMUM_MAPPING_QUALITY) {
                    continue;
                }

                final int absInsertSize = Math.abs(sam.getInferredInsertSize());
                if (absInsertSize > chimeraSizeMinimum) {
                    superSized++;
                } else if (sam.getMateNegativeStrandFlag() == sam.getReadNegativeStrandFlag()) {
                    tandemPairs++;
                } else if (!sam.getMateReferenceIndex().equals(sam.getReferenceIndex())) {
                    crossChromPairs++;
                } else {
                    final PairOrientation pairOrientation = SamPairUtil.getPairOrientation(sam);
                    if (pairOrientation == PairOrientation.RF) {
                        outieHistogram.increment(absInsertSize);
                        outies++;
                        if (sam.getDuplicateReadFlag()) {
                            outieDupes++;
                        }
                    } else if (pairOrientation == PairOrientation.FR) {
                        innieHistogram.increment(absInsertSize);
                        innies++;
                        if (sam.getDuplicateReadFlag()) {
                            innieDupes++;
                        }
                    } else {
                        throw new IllegalStateException("This should never happen");
                    }
                }
            }
            CloserUtil.close(reader);
        }

        MetricsFile<JumpingLibraryMetrics, Integer> metricsFile = getMetricsFile();
        JumpingLibraryMetrics metrics = new JumpingLibraryMetrics();
        metrics.JUMP_PAIRS = outies;
        metrics.JUMP_DUPLICATE_PAIRS = outieDupes;
        metrics.JUMP_DUPLICATE_PCT = outies != 0 ? outieDupes / (double) outies : 0;
        metrics.JUMP_LIBRARY_SIZE = (outies > 0 && outieDupes > 0) ? DuplicationMetrics.estimateLibrarySize(outies, outies - outieDupes) : 0;
        outieHistogram.trimByTailLimit(TAIL_LIMIT);
        metrics.JUMP_MEAN_INSERT_SIZE = outieHistogram.getMean();
        metrics.JUMP_STDEV_INSERT_SIZE = outieHistogram.getStandardDeviation();
        metrics.NONJUMP_PAIRS = innies;
        metrics.NONJUMP_DUPLICATE_PAIRS = innieDupes;
        metrics.NONJUMP_DUPLICATE_PCT = innies != 0 ? innieDupes / (double) innies : 0;
        metrics.NONJUMP_LIBRARY_SIZE = (innies > 0 && innieDupes > 0) ? DuplicationMetrics.estimateLibrarySize(innies, innies - innieDupes) : 0;
        innieHistogram.trimByTailLimit(TAIL_LIMIT);
        metrics.NONJUMP_MEAN_INSERT_SIZE = innieHistogram.getMean();
        metrics.NONJUMP_STDEV_INSERT_SIZE = innieHistogram.getStandardDeviation();
        metrics.CHIMERIC_PAIRS = crossChromPairs + superSized + tandemPairs;
        metrics.FRAGMENTS = fragments;
        double totalPairs = outies + innies + metrics.CHIMERIC_PAIRS;
        metrics.PCT_JUMPS = totalPairs != 0 ? outies / totalPairs : 0;
        metrics.PCT_NONJUMPS = totalPairs != 0 ? innies / totalPairs : 0;
        metrics.PCT_CHIMERAS = totalPairs != 0 ? metrics.CHIMERIC_PAIRS / totalPairs : 0;
        metricsFile.addMetric(metrics);
        metricsFile.write(OUTPUT);

        return 0;
    }

    /**
     * Calculates the mode for outward-facing pairs, using the first SAMPLE_FOR_MODE
     * outward-facing pairs found in INPUT
     */
    private double getOutieMode() {

        int samplePerFile = SAMPLE_FOR_MODE / INPUT.size();

        Histogram<Integer> histo = new Histogram<Integer>();

        for (File f : INPUT) {
            SamReader reader = SamReaderFactory.makeDefault().open(f);
            int sampled = 0;
            for (Iterator<SAMRecord> it = reader.iterator(); it.hasNext() && sampled < samplePerFile; ) {
                SAMRecord sam = it.next();
                if (!sam.getFirstOfPairFlag()) {
                    continue;
                }
                // If we get here we've hit the end of the aligned reads
                if (sam.getReadUnmappedFlag() && sam.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
                    break;
                } else if (sam.getReadUnmappedFlag() || sam.getMateUnmappedFlag()) {
                    continue;
                } else if ((sam.getAttribute(SAMTag.MQ.name()) == null ||
                        sam.getIntegerAttribute(SAMTag.MQ.name()) >= MINIMUM_MAPPING_QUALITY) &&
                        sam.getMappingQuality() >= MINIMUM_MAPPING_QUALITY &&
                        sam.getMateNegativeStrandFlag() != sam.getReadNegativeStrandFlag() &&
                        sam.getMateReferenceIndex().equals(sam.getReferenceIndex()) &&
                        SamPairUtil.getPairOrientation(sam) == PairOrientation.RF) {
                    histo.increment(Math.abs(sam.getInferredInsertSize()));
                    sampled++;
                }
            }
            CloserUtil.close(reader);
        }

        return histo.size() > 0 ? histo.getMode() : 0;
    }
}
