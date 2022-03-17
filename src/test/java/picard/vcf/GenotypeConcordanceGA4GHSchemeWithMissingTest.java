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

package picard.vcf;

import htsjdk.samtools.util.FormatUtil;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.vcf.GenotypeConcordanceStates.TruthAndCallStates;
import java.io.File;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.List;

/**
 * Created by kbergin on 6/24/15.
 */
public class GenotypeConcordanceGA4GHSchemeWithMissingTest{
    private static final File TEST_DATA_PATH = new File("testdata/picard/vcf/");
    final File ceuTrioSNPSVCF = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");
    final File nistTruthVCF = new File(TEST_DATA_PATH, "NIST.selected.vcf");
    private final GenotypeConcordanceSchemeFactory schemeFactory = new GenotypeConcordanceSchemeFactory();
    final GenotypeConcordanceScheme scheme = schemeFactory.getScheme(true);
    final List<File> intervalList = Collections.singletonList(new File(TEST_DATA_PATH, "IntervalList1PerChrom.interval_list"));

    @Test
    public void testMissingHomRefSchemeWithIntervals() {
        final GenotypeConcordanceCounts concordanceCounts = GenotypeConcordanceTest.getGenotypeConcordanceCounts(nistTruthVCF, ceuTrioSNPSVCF, "NA12878", true, intervalList);
        concordanceCounts.validateCountsAgainstScheme(scheme);

        final Map<GenotypeConcordanceStates.TruthAndCallStates, Integer> nonZeroCounts = new HashMap<GenotypeConcordanceStates.TruthAndCallStates, Integer>();

        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.MISSING, GenotypeConcordanceStates.CallState.MISSING), 63001594);
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.MISSING, GenotypeConcordanceStates.CallState.HET_REF_VAR1), 1);
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.MISSING, GenotypeConcordanceStates.CallState.VC_FILTERED), 2);
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.HET_REF_VAR1, GenotypeConcordanceStates.CallState.MISSING), 1);
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.HET_REF_VAR1, GenotypeConcordanceStates.CallState.HET_REF_VAR1), 11);
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.HOM_VAR1, GenotypeConcordanceStates.CallState.MISSING), 4);
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.HOM_VAR1, GenotypeConcordanceStates.CallState.HOM_VAR1), 5);

        GenotypeConcordanceTest.assertNonZeroCountsAgree(concordanceCounts, nonZeroCounts);

        final FormatUtil fmt = new FormatUtil();

        Assert.assertEquals(fmt.format(concordanceCounts.getSensitivity(scheme, GenotypeConcordanceCounts.HET_TRUTH_STATES)), "0.916667");
        Assert.assertEquals(fmt.format(concordanceCounts.Ppv(scheme, GenotypeConcordanceCounts.HET_CALL_STATES)), "0.916667");
        Assert.assertEquals(fmt.format(concordanceCounts.getSensitivity(scheme, GenotypeConcordanceCounts.HOM_VAR_TRUTH_STATES)), "0.555556");
        Assert.assertEquals(fmt.format(concordanceCounts.Ppv(scheme, GenotypeConcordanceCounts.HOM_VAR_CALL_STATES)), "1");
        Assert.assertEquals(fmt.format(concordanceCounts.getSpecificity(scheme, GenotypeConcordanceCounts.HOM_VAR_TRUTH_STATES)), "?");
        Assert.assertEquals(fmt.format(concordanceCounts.getSensitivity(scheme, GenotypeConcordanceCounts.VAR_TRUTH_STATES)), "0.761905");
        Assert.assertEquals(fmt.format(concordanceCounts.Ppv(scheme, GenotypeConcordanceCounts.VAR_CALL_STATES)), "0.941176");
        Assert.assertEquals(fmt.format(concordanceCounts.getSpecificity(scheme, GenotypeConcordanceCounts.VAR_TRUTH_STATES)), "1");
        Assert.assertEquals(fmt.format(concordanceCounts.calculateGenotypeConcordance(scheme, true)), "1");
        Assert.assertEquals(fmt.format(concordanceCounts.calculateNonRefGenotypeConcordance(scheme, true)), "0.727273");
    }

}
