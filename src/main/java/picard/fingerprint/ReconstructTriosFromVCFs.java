/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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
package picard.fingerprint;

import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.util.Collections;
import java.util.Map;
import java.util.concurrent.TimeUnit;

/**
 *
 * @author Tim Fennell
 * @author Jonathan Barlev
 * @author Yossi Farjoun
 */
@CommandLineProgramProperties(
        summary = "A program that infers trios in a cohort given by a VCF of variant calls. The program infers the sex of each sample " +
                "(by looking at AD fields in the VCF), find possible parents for each sample (by calculating the probability that " +
                "other sample is a parent), and then from the list of possible mothers and fathers tries to come up with the most likely trio.",
        oneLineSummary = "Infers trios in a cohort given by a VCF.",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@ExperimentalFeature
public class ReconstructTriosFromVCFs extends ReconstructTrios {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input VCF file.")
    public File INPUT;

    @Argument(shortName = "AUTO_V", doc = "The number of variants which will be sampled for coverage on each non-sex chromosome. " +
            "(Attempt to) sample coverage at AUTOSOMAL_VARIANTS evenly spaced variants on each non-Sex chromosome.")
    public int AUTOSOMAL_VARIANTS = 1000;

    @Argument(shortName = "ALLO_V", doc = "The number of variants which will be sampled for coverage on each sex chromosome. " +
            "(Attempt to) sample coverage at ALLOSOMAL_VARIANTS evenly spaced variants on each Sex chromosome. ")
    public int ALLOSOMAL_VARIANTS = 10000;

    @Argument(shortName = "PARIL", doc = "An Interval List contain the pseudoautosomal region. " +
            "Used for masking out that region since the regions on one chromosome will get mapped to the other and thus can change the apparent coverage.", optional = true)
    public File PSEUDOAUTOSOMAL_REGION_INTERVAL_LIST = null;

    @Argument(doc = "Some VCFs do not have variants on Male chromosomes. Consequently, when determining sex, one may wish to give less weight to their coverage than" +
            "the female chromosomes. We do this by dividing the male chromosome coverage values by a factor of Y_COVERAGE_SHRINK_FACTOR.")
    public double Y_COVERAGE_SHRINK_FACTOR = 2.0;

    @Argument(doc = "The number of threads to use when fingerprinting files.")
    public int NUM_THREADS = 1;

    @Argument(doc = "The maximum runtime in hours that should be used for fingerprinting files.")
    public int MAX_RUNTIME_HOURS = 5 * 24;

    @Override
    protected SexInferenceEngine getSexInferencer() {
        return new SexInferenceEngineFromVCF(MALE_CHROMOSOMES, FEMALE_CHROMOSOMES, INPUT, AUTOSOMAL_VARIANTS, ALLOSOMAL_VARIANTS, PSEUDOAUTOSOMAL_REGION_INTERVAL_LIST, Y_COVERAGE_SHRINK_FACTOR);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        return super.doWork();
    }

    @Override
    public Map<FingerprintIdDetails, Fingerprint> getFingerprints() {
        final FingerprintChecker checker = new FingerprintChecker(hapmap);
        return checker.fingerprintFiles(Collections.singleton(INPUT.toPath()), NUM_THREADS, MAX_RUNTIME_HOURS, TimeUnit.HOURS);
    }
}
