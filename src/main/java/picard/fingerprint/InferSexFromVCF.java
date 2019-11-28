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

package picard.fingerprint;

import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;

/**
 * Program to Infer Sex of a cohort of samples from a vcf. Needs a large cohort with both sexes to
 * work since it uses a clustering-based algorithm
 *
 * See InferSex for more details.
 *
 * This class looks at depth markers on spaced variants in the VCF and provide information from that to the base class.
 *
 * Created by jbarlev on 5/30/14.
 *
 * @author Jonathan Barlev
 * @author Yossi Farjoun
 */
@CommandLineProgramProperties(
        summary = "A program that can infer sample sex from a VCF file." +
                "It compares the coverage over the male and female chromosomes to that over the rest of the" +
                "genome, and performs clustering to find the answer.",
        oneLineSummary = "Infer sample sex from a VCF file",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@ExperimentalFeature
public class InferSexFromVCF extends InferSex {
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input VCF file.")
    public File INPUT;

    @Argument(shortName = "AUTO_V", doc = "Number determining how many variants are sampled for coverage on each non-sex chromosome. " +
            "Program will (attempt to) sample coverage at AUTOSOMAL_VARIANTS evenly spaced variants on each non-sex chromosome.")
    public int AUTOSOMAL_VARIANTS = 10;

    @Argument(shortName = "ALLO_V", doc = "Number determining how many variants are sampled for coverage on each sex chromosome. " +
            "Program will (attempt to) sample coverage at ALLOSOMAL_VARIANTS evenly spaced variants on each sex chromosome. ")
    public int ALLOSOMAL_VARIANTS = 100;

    @Argument(shortName = "PAR", doc = "An IntervalList containing the pseudoautosomal region. Used for masking out that region since the regions on one " +
            "chromosome may get mapped to the other and thus can change the apparent coverage.", optional = true)
    public File PSEUDOAUTOSOMAL_REGION = null;

    @Argument(doc="Some VCFs do not have variants (called) on Male chromosomes. Consequently, when determining sex, one may wish to give less weight to their Coverage than" +
            "the memale chromosomes. Do this by dividing the male chromosome coverage values by a factor of Y_COVERAGE_SHRINK_FACTOR.")
    public double Y_COVERAGE_SHRINK_FACTOR = 2.0;

    @Override
    SexInferenceEngine getSexInference() {
        return new SexInferenceEngineFromVCF(MALE_CHROMS, FEMALE_CHROMS, INPUT, AUTOSOMAL_VARIANTS, ALLOSOMAL_VARIANTS, PSEUDOAUTOSOMAL_REGION, Y_COVERAGE_SHRINK_FACTOR);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        return super.doWork();
    }
}
