/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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

package picard.sam.markduplicates;

import htsjdk.samtools.DuplicateSet;
import htsjdk.samtools.DuplicateSetIterator;
import htsjdk.samtools.SAMRecordDuplicateComparator;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * <p>This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are
 * defined as originating from a single fragment of DNA. It is based on the {@link MarkDuplicatesWithMateCigar} tool, with added logic
 * to leverage Unique Molecular Identifier (UMI) information.</p>
 *
 * <p>It makes use of the fact that duplicate sets with UMIs can be broken up into subsets based on
 * information contained in the UMI.  In addition to assuming that all members of a duplicate set must have the same start and end position, it imposes that
 * they must also have sufficiently similar UMIs. In this context, 'sufficiently similar' is parameterized by the command line
 * argument MAX_EDIT_DISTANCE_TO_JOIN, which sets the edit distance between UMIs that will be considered to be part of the same
 * original molecule. This logic allows for sequencing errors in UMIs.</p>
 * <p> If UMIs contain dashes, the dashes will be ignored. If UMIs contain Ns, these UMIs will not contribute to UMI metrics
 * associated with each record. If the MAX_EDIT_DISTANCE_TO_JOIN allows, UMIs with Ns will be included in the duplicate set and
 * the UMI metrics associated with each duplicate set. Ns are counted as an edit distance from other bases {ATCG}, but are not
 * considered different from each other.</p>
 * <p>This tool is NOT intended to be used on data without UMIs; for marking duplicates in non-UMI data, see {@link MarkDuplicates} or
 * {@link MarkDuplicatesWithMateCigar}. Mixed data (where some reads have UMIs and others do not) is not supported.</p>
 *
 * <p>Note also that this tool will not work with alignments that have large gaps or deletions, such as those from RNA-seq data.
 * This is due to the need to buffer small genomic windows to ensure integrity of the duplicate marking, while large skips
 * (ex. skipping introns) in the alignment records would force making that window very large, thus exhausting memory. </p>
 *
 * <p>Note: Metrics labeled as percentages are actually expressed as fractions!</p>
 * <h4>Usage example:</h4>
 * <pre>
 *
 * java -jar picard.jar UmiAwareMarkDuplicatesWithMateCigar \\<br />
 *      I=input.bam \\<br />
 *      O=output.bam \\<br />
 *      M=output_duplicate_metrics.txt \\<br />
 *      UMI_METRICS=output_umi_metrics.txt
 * </pre>
 * <hr />
 *
 * @author fleharty
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = UmiAwareMarkDuplicatesWithMateCigar.USAGE_SUMMARY + UmiAwareMarkDuplicatesWithMateCigar.USAGE_DETAILS,
        oneLineSummary = UmiAwareMarkDuplicatesWithMateCigar.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class
)
public class UmiAwareMarkDuplicatesWithMateCigar extends SimpleMarkDuplicatesWithMateCigar {
    static final String USAGE_SUMMARY = "Identifies duplicate reads using information from read positions and UMIs. ";
    static final String USAGE_DETAILS = "<p>This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are" +
            "defined as originating from a single fragment of DNA. It is based on the {@link MarkDuplicatesWithMateCigar} tool, with added logic" +
            "to leverage Unique Molecular Identifier (UMI) information.</p>" +
            "<p>It makes use of the fact that duplicate sets with UMIs can be broken up into subsets based on" +
            "information contained in the UMI.  In addition to assuming that all members of a duplicate set must have the same start and end position, it imposes that" +
            "they must also have sufficiently similar UMIs. In this context, 'sufficiently similar' is parameterized by the command line" +
            "argument MAX_EDIT_DISTANCE_TO_JOIN, which sets the edit distance between UMIs that will be considered to be part of the same" +
            "original molecule. This logic allows for sequencing errors in UMIs.</p>" +
            "<p> If UMIs contain dashes, the dashes will be ignored. If UMIs contain Ns, these UMIs will not contribute to UMI metrics" +
            "associated with each record. If the MAX_EDIT_DISTANCE_TO_JOIN allows, UMIs with Ns will be included in the duplicate set and" +
            "the UMI metrics associated with each duplicate set. Ns are counted as an edit distance from other bases {ATCG}, but are not" +
            "considered different from each other.</p>" +
            "<p>This tool is NOT intended to be used on data without UMIs; for marking duplicates in non-UMI data, see {@link MarkDuplicates} or" +
            "{@link MarkDuplicatesWithMateCigar}. Mixed data (where some reads have UMIs and others do not) is not supported.</p>" +
            "<p>Note also that this tool will not work with alignments that have large gaps or deletions, such as those from RNA-seq data." +
            "This is due to the need to buffer small genomic windows to ensure integrity of the duplicate marking, while large skips" +
            "(ex. skipping introns) in the alignment records would force making that window very large, thus exhausting memory. </p>" +
            "<p>Note: Metrics labeled as percentages are actually expressed as fractions!</p>" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar UmiAwareMarkDuplicatesWithMateCigar <br />" +
            "      I=input.bam <br />" +
            "      O=output.bam <br />" +
            "      M=output_duplicate_metrics.txt <br />" +
            "      UMI_METRICS=output_umi_metrics.txt" +
            "</pre>" +
            "<hr />";

    @Argument(shortName = "MAX_EDIT_DISTANCE_TO_JOIN", doc = "Largest edit distance that UMIs must have in order to be considered as coming from distinct source molecules.", optional = true)
    public int MAX_EDIT_DISTANCE_TO_JOIN = 1;

    // The UMI_METRICS file provides various statistical measurements collected about the UMIs during deduplication.
    @Argument(shortName = "UMI_METRICS", doc = "UMI Metrics")
    public File UMI_METRICS_FILE;

    @Argument(shortName = "UMI_TAG_NAME", doc = "Tag name to use for UMI", optional = true)
    public String UMI_TAG_NAME = "RX";

    @Argument(shortName = "ASSIGNED_UMI_TAG", doc = "Tag name to use for assigned UMI", optional = true)
    public String ASSIGNED_UMI_TAG = "MI";

    // Since we inherit from SimpleMarkDuplicatesWithMateCigar, it is useful for us to also inherit the tests
    // which do not contain UMIs.  By default, we don't allow for missing UMIs, but for the inherited tests
    // we allow for missing UMIs.
    @Argument(doc = "FOR TESTING ONLY: allow for missing UMIs if data doesn't have UMIs. This option is intended to be used ONLY for testing the code. Use MarkDuplicatesWithMateCigar if data has no UMIs. Mixed data (where some reads have UMIs and others do not) is not supported.", optional = true)
    public boolean ALLOW_MISSING_UMIS = false;

    private final Log log = Log.getInstance(UmiAwareMarkDuplicatesWithMateCigar.class);
    private final Map<String, UmiMetrics> metrics = new HashMap<>();

    @Override
    protected int doWork() {
        // Before we do anything, make sure the UMI_METRICS_FILE can be written to.
        IOUtil.assertFileIsWritable(UMI_METRICS_FILE);

        // Perform Mark Duplicates work
        final int retval = super.doWork();

        // Add results in metrics to the metricsFile
        final MetricsFile<UmiMetrics, Double> metricsFile = getMetricsFile();
        for (final UmiMetrics metric : metrics.values()) {
            metricsFile.addMetric(metric);
        }

        metricsFile.write(UMI_METRICS_FILE);
        return retval;
    }

    @Override
    protected CloseableIterator<DuplicateSet> getDuplicateSetIterator(final SamHeaderAndIterator headerAndIterator, final SAMRecordDuplicateComparator comparator) {
        return new UmiAwareDuplicateSetIterator(
                    new DuplicateSetIterator(headerAndIterator.iterator,
                    headerAndIterator.header,
                    false,
                    comparator), MAX_EDIT_DISTANCE_TO_JOIN, UMI_TAG_NAME, ASSIGNED_UMI_TAG, ALLOW_MISSING_UMIS, metrics);
    }
}
