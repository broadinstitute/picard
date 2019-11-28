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

import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

/**
 * A Program that infers trios in a cohort from a collection of BAM files",
 */
@CommandLineProgramProperties(
        summary = "Attempts to use a set of genome data samples to determine trio relationships within the cohort provided.\n" +
                " " +
                "Roughly speaking it:\n" +
                " - Determines genders\n" +
                " - Performs pairwise tests between samples for 1st vs. 2nd degree relatedness\n" +
                " - Tests all plausible trios from 1st degree related samples of appropriate genders\n" +
                "Emits likely trios in PED file.",
        oneLineSummary = "Infers trios in a cohort from a collection of BAM files.",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@ExperimentalFeature
public class ReconstructTriosFromBAMs extends ReconstructTrios {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input BAM files, or files listing input file locations.")
    public List<File> INPUT;

    @Argument(doc = "The number of threads to use when fingerprinting files.")
    public int NUM_THREADS = 1;

    @Argument(doc = "The maximum runtime in hours that should be used for fingerprinting files.")
    public int MAX_RUNTIME_HOURS = 5 * 24;

    @Override
    protected SexInferenceEngine getSexInferencer() {
        return new SexInferenceEngineFromBAM(MALE_CHROMOSOMES, FEMALE_CHROMOSOMES, INPUT);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFilesAreReadable(INPUT);
        return super.doWork();
    }

    @Override
    public final Map<FingerprintIdDetails, Fingerprint> getFingerprints() {
        final List<Path> bams = IOUtil.unrollPaths(INPUT.stream().map(File::toPath).collect(Collectors.toList()), FileExtensions.BAM);
        final FingerprintChecker checker = new FingerprintChecker(hapmap);
        final Map<FingerprintIdDetails, Fingerprint> fpMap = checker.fingerprintFiles(bams, NUM_THREADS, MAX_RUNTIME_HOURS, TimeUnit.HOURS);
        return Fingerprint.mergeFingerprintsBy(fpMap, Fingerprint.getFingerprintIdDetailsStringFunction(CrosscheckMetric.DataType.SAMPLE));
    }
}
