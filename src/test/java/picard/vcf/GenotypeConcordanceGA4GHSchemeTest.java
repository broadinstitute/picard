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
public class GenotypeConcordanceGA4GHSchemeTest{
    private static final File TEST_DATA_PATH = new File("testdata/picard/vcf/");
    final File ceuTrioSnpsVcf = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");
    private final GenotypeConcordanceSchemeFactory schemeFactory = new GenotypeConcordanceSchemeFactory();
    final GenotypeConcordanceScheme scheme = schemeFactory.getScheme(false);
    final List<File> intervalList = Collections.singletonList(new File(TEST_DATA_PATH, "IntervalListChr1Small.interval_list"));

    @Test
    public void testGA4GHScheme() throws Exception {

        final GenotypeConcordanceCounts concordanceCounts = GenotypeConcordanceTest.getGenotypeConcordanceCounts(ceuTrioSnpsVcf, ceuTrioSnpsVcf, "NA12878", false, null);
        concordanceCounts.validateCountsAgainstScheme(scheme);

        final Map<GenotypeConcordanceStates.TruthAndCallStates, Integer> nonZeroCounts = new HashMap<GenotypeConcordanceStates.TruthAndCallStates, Integer>();
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.HET_REF_VAR1, GenotypeConcordanceStates.CallState.HET_REF_VAR1), 104);
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.HOM_VAR1, GenotypeConcordanceStates.CallState.HOM_VAR1), 59);
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.VC_FILTERED, GenotypeConcordanceStates.CallState.VC_FILTERED), 40);

        GenotypeConcordanceTest.assertNonZeroCountsAgree(concordanceCounts, nonZeroCounts);

        final FormatUtil fmt = new FormatUtil();

        Assert.assertEquals(fmt.format(concordanceCounts.getSensitivity(scheme, GenotypeConcordanceCounts.HET_TRUTH_STATES)), "1");
        Assert.assertEquals(fmt.format(concordanceCounts.Ppv(scheme, GenotypeConcordanceCounts.HET_CALL_STATES)), "1");
        Assert.assertEquals(fmt.format(concordanceCounts.getSensitivity(scheme, GenotypeConcordanceCounts.HOM_VAR_TRUTH_STATES)), "1");
        Assert.assertEquals(fmt.format(concordanceCounts.Ppv(scheme, GenotypeConcordanceCounts.HOM_VAR_CALL_STATES)), "1");
        Assert.assertEquals(fmt.format(concordanceCounts.getSpecificity(scheme, GenotypeConcordanceCounts.HOM_VAR_TRUTH_STATES)), "?");
        Assert.assertEquals(fmt.format(concordanceCounts.getSensitivity(scheme, GenotypeConcordanceCounts.VAR_TRUTH_STATES)), "1");
        Assert.assertEquals(fmt.format(concordanceCounts.Ppv(scheme, GenotypeConcordanceCounts.VAR_CALL_STATES)), "1");
        Assert.assertEquals(fmt.format(concordanceCounts.getSpecificity(scheme, GenotypeConcordanceCounts.VAR_TRUTH_STATES)), "1");
        Assert.assertEquals(fmt.format(concordanceCounts.calculateGenotypeConcordance(scheme, true)), "1");
        Assert.assertEquals(fmt.format(concordanceCounts.calculateNonRefGenotypeConcordance(scheme, true)), "1");
    }

    @Test
    public void testGA4GHSchemeDiffSamples() throws Exception {

        final GenotypeConcordanceCounts concordanceCounts = GenotypeConcordanceTest.getGenotypeConcordanceCounts(ceuTrioSnpsVcf, ceuTrioSnpsVcf, "NA12891", false, null);
        concordanceCounts.validateCountsAgainstScheme(scheme);

        final Map<GenotypeConcordanceStates.TruthAndCallStates, Integer> nonZeroCounts = new HashMap<GenotypeConcordanceStates.TruthAndCallStates, Integer>();
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.HOM_REF, GenotypeConcordanceStates.CallState.HET_REF_VAR1), 31);
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.HET_REF_VAR1, GenotypeConcordanceStates.CallState.HOM_REF), 30);
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.HET_REF_VAR1, GenotypeConcordanceStates.CallState.HET_REF_VAR1), 50);
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.HET_REF_VAR1, GenotypeConcordanceStates.CallState.HOM_VAR1), 24);
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.HOM_VAR1, GenotypeConcordanceStates.CallState.HET_REF_VAR1), 18);
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.HOM_VAR1, GenotypeConcordanceStates.CallState.HOM_VAR1), 41);
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.VC_FILTERED, GenotypeConcordanceStates.CallState.VC_FILTERED), 49);

        GenotypeConcordanceTest.assertNonZeroCountsAgree(concordanceCounts, nonZeroCounts);

        final FormatUtil fmt = new FormatUtil();

        Assert.assertEquals(fmt.format(concordanceCounts.getSensitivity(scheme, GenotypeConcordanceCounts.HET_TRUTH_STATES)), "0.711538");
        Assert.assertEquals(fmt.format(concordanceCounts.Ppv(scheme, GenotypeConcordanceCounts.HET_CALL_STATES)), "0.686869");
        Assert.assertEquals(fmt.format(concordanceCounts.getSensitivity(scheme, GenotypeConcordanceCounts.HOM_VAR_TRUTH_STATES)), "0.766234");
        Assert.assertEquals(fmt.format(concordanceCounts.Ppv(scheme, GenotypeConcordanceCounts.HOM_VAR_CALL_STATES)), "0.730337");
        Assert.assertEquals(fmt.format(concordanceCounts.getSpecificity(scheme, GenotypeConcordanceCounts.HOM_VAR_TRUTH_STATES)), "?");
        Assert.assertEquals(fmt.format(concordanceCounts.getSensitivity(scheme, GenotypeConcordanceCounts.VAR_TRUTH_STATES)), "0.734807");
        Assert.assertEquals(fmt.format(concordanceCounts.Ppv(scheme, GenotypeConcordanceCounts.VAR_CALL_STATES)), "0.707447");
        Assert.assertEquals(fmt.format(concordanceCounts.getSpecificity(scheme, GenotypeConcordanceCounts.VAR_TRUTH_STATES)), "0.668675");
        Assert.assertEquals(fmt.format(concordanceCounts.calculateGenotypeConcordance(scheme, true)), "0.576132");
        Assert.assertEquals(fmt.format(concordanceCounts.calculateNonRefGenotypeConcordance(scheme, true)), "0.469072");
    }

    @Test
    public void testGA4GHSchemeWithIntervals() throws Exception {

        final GenotypeConcordanceCounts concordanceCounts = GenotypeConcordanceTest.getGenotypeConcordanceCounts(ceuTrioSnpsVcf, ceuTrioSnpsVcf, "NA12878", false, intervalList);
        concordanceCounts.validateCountsAgainstScheme(scheme);

        final Map<GenotypeConcordanceStates.TruthAndCallStates, Integer> nonZeroCounts = new HashMap<GenotypeConcordanceStates.TruthAndCallStates, Integer>();
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.HET_REF_VAR1, GenotypeConcordanceStates.CallState.HET_REF_VAR1), 1);
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.VC_FILTERED, GenotypeConcordanceStates.CallState.VC_FILTERED), 2);

        GenotypeConcordanceTest.assertNonZeroCountsAgree(concordanceCounts, nonZeroCounts);

        final FormatUtil fmt = new FormatUtil();

        Assert.assertEquals(fmt.format(concordanceCounts.getSensitivity(scheme, GenotypeConcordanceCounts.HET_TRUTH_STATES)), "1");
        Assert.assertEquals(fmt.format(concordanceCounts.Ppv(scheme, GenotypeConcordanceCounts.HET_CALL_STATES)), "1");
        Assert.assertEquals(fmt.format(concordanceCounts.getSensitivity(scheme, GenotypeConcordanceCounts.HOM_VAR_TRUTH_STATES)), "?");
        Assert.assertEquals(fmt.format(concordanceCounts.Ppv(scheme, GenotypeConcordanceCounts.HOM_VAR_CALL_STATES)), "?");
        Assert.assertEquals(fmt.format(concordanceCounts.getSpecificity(scheme, GenotypeConcordanceCounts.HOM_VAR_TRUTH_STATES)), "?");
        Assert.assertEquals(fmt.format(concordanceCounts.getSensitivity(scheme, GenotypeConcordanceCounts.VAR_TRUTH_STATES)), "1");
        Assert.assertEquals(fmt.format(concordanceCounts.Ppv(scheme, GenotypeConcordanceCounts.VAR_CALL_STATES)), "1");
        Assert.assertEquals(fmt.format(concordanceCounts.getSpecificity(scheme, GenotypeConcordanceCounts.VAR_TRUTH_STATES)), "1");
        Assert.assertEquals(fmt.format(concordanceCounts.calculateGenotypeConcordance(scheme, true)), "1");
        Assert.assertEquals(fmt.format(concordanceCounts.calculateNonRefGenotypeConcordance(scheme, true)), "1");
    }

    @Test
    public void testGA4GHSchemeDiffSamplesWithIntervals() throws Exception {

        final GenotypeConcordanceCounts concordanceCounts = GenotypeConcordanceTest.getGenotypeConcordanceCounts(ceuTrioSnpsVcf, ceuTrioSnpsVcf, "NA12891", false, intervalList);
        concordanceCounts.validateCountsAgainstScheme(scheme);

        final Map<GenotypeConcordanceStates.TruthAndCallStates, Integer> nonZeroCounts = new HashMap<GenotypeConcordanceStates.TruthAndCallStates, Integer>();
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.HOM_REF, GenotypeConcordanceStates.CallState.HET_REF_VAR1), 1);
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.HET_REF_VAR1, GenotypeConcordanceStates.CallState.HET_REF_VAR1), 1);
        nonZeroCounts.put(new TruthAndCallStates(GenotypeConcordanceStates.TruthState.VC_FILTERED, GenotypeConcordanceStates.CallState.VC_FILTERED), 2);

        GenotypeConcordanceTest.assertNonZeroCountsAgree(concordanceCounts, nonZeroCounts);

        final FormatUtil fmt = new FormatUtil();

        Assert.assertEquals(fmt.format(concordanceCounts.getSensitivity(scheme, GenotypeConcordanceCounts.HET_TRUTH_STATES)), "1");
        Assert.assertEquals(fmt.format(concordanceCounts.Ppv(scheme, GenotypeConcordanceCounts.HET_CALL_STATES)), "0.5");
        Assert.assertEquals(fmt.format(concordanceCounts.getSensitivity(scheme, GenotypeConcordanceCounts.HOM_VAR_TRUTH_STATES)), "?");
        Assert.assertEquals(fmt.format(concordanceCounts.Ppv(scheme, GenotypeConcordanceCounts.HOM_VAR_CALL_STATES)), "?");
        Assert.assertEquals(fmt.format(concordanceCounts.getSpecificity(scheme, GenotypeConcordanceCounts.HOM_VAR_TRUTH_STATES)), "?");
        Assert.assertEquals(fmt.format(concordanceCounts.getSensitivity(scheme, GenotypeConcordanceCounts.VAR_TRUTH_STATES)), "1");
        Assert.assertEquals(fmt.format(concordanceCounts.Ppv(scheme, GenotypeConcordanceCounts.VAR_CALL_STATES)), "0.5");
        Assert.assertEquals(fmt.format(concordanceCounts.getSpecificity(scheme, GenotypeConcordanceCounts.VAR_TRUTH_STATES)), "0.666667");
        Assert.assertEquals(fmt.format(concordanceCounts.calculateGenotypeConcordance(scheme, true)), "0.75");
        Assert.assertEquals(fmt.format(concordanceCounts.calculateNonRefGenotypeConcordance(scheme, true)), "0.5");
    }
}