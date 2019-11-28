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
import java.util.List;

/**
 * Program to Infer Sex of a cohort of samples from bams. Needs a large cohort with both sexes to
 * work since it uses a clustering-based algorithm
 *
 * See InferSex for more details.
 *
 * This class looks at the index information to get average coverage over the sex chromosomes.
 * it provides that information to the base class.
 *
 * Created by jbarlev on 5/30/14.
 *
 * @author Jonathan Barlev
 * @author Yossi Farjoun
 */
@CommandLineProgramProperties(
        summary = "A program that can infer sample sex from a collection of BAM files." +
                "It compares the coverage over the male and female chromosomes to that over the rest of the" +
                "genome, and performs clustering to find the answer. It uses the bam index to avoid having to " +
                "iterate over the whole BAM file.",
        oneLineSummary = "Infer sample sex from a collection of BAM files",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@ExperimentalFeature
public class InferSexFromBAM extends InferSex {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input BAM (not SAM) file. Must be indexed.")
    public List<File> INPUT;

    @Override
    SexInferenceEngine getSexInference() {
        return new SexInferenceEngineFromBAM(MALE_CHROMS, FEMALE_CHROMS, INPUT);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFilesAreReadable(INPUT);

        return super.doWork();
    }
}
