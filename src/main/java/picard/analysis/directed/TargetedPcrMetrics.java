/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
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

package picard.analysis.directed;


import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.util.help.HelpConstants;

/** Metrics class for the analysis of reads obtained from targeted pcr experiments e.g. the TruSeq Custom Amplicon
 * (TSCA) kit (Illumina).  */
@DocumentedFeature(groupName = HelpConstants.DOC_CAT_METRICS, summary = HelpConstants.DOC_CAT_METRICS_SUMMARY)
public class TargetedPcrMetrics extends TargetMetricsBase {

    /**  The name of the amplicon set used in this metrics collection run */
    public String CUSTOM_AMPLICON_SET;

    /** The number of unique bases covered by the intervals of all amplicons in the amplicon set */
    public long AMPLICON_TERRITORY;

    /** The number of PF_BASES_ALIGNED that mapped to an amplified region of the genome. */
    public long ON_AMPLICON_BASES;

    /** The number of PF_BASES_ALIGNED that mapped to within a fixed interval of an amplified region, but not on a
     * baited region. */
    public long NEAR_AMPLICON_BASES;

    /** The number of PF_BASES_ALIGNED that mapped neither on or near an amplicon. */
    public long OFF_AMPLICON_BASES;

    /** The fraction of PF_BASES_ALIGNED that mapped to or near an amplicon, (ON_AMPLICON_BASES +
     * NEAR_AMPLICON_BASES)/PF_BASES_ALIGNED. */
    public double PCT_AMPLIFIED_BASES;

    /** The fraction of PF_BASES_ALIGNED that mapped neither onto or near an amplicon,
     * OFF_AMPLICON_BASES/PF_BASES_ALIGNED */
    public double PCT_OFF_AMPLICON;

    /**
     * The fraction of bases mapping to regions on or near amplicons, which mapped directly to but not near
     * amplicons, ON_AMPLICON_BASES/(NEAR_AMPLICON_BASES + ON_AMPLICON_BASES)
     * */
    public double ON_AMPLICON_VS_SELECTED;

    /** The mean read coverage of all amplicon regions in the experiment. */
    public double MEAN_AMPLICON_COVERAGE;

    /** The fold by which the amplicon region has been amplified above genomic background. */
    public double FOLD_ENRICHMENT;
}
