/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
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

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

/**
 * Calculates vrious metrics on a sample fingerprint, indicating whether the fingerprint satisfies the assumptions we have.
 * For example, if too many sites are hetrozygous, that would get flagged.
 *
 * @author Yossi Farjoun
 */

@CommandLineProgramProperties(
        summary = CalculateFingerprintMetrics.USAGE_SUMMARY + CalculateFingerprintMetrics.USAGE_DETAILS,
        oneLineSummary = CalculateFingerprintMetrics.USAGE_SUMMARY,
        programGroup = DiagnosticsAndQCProgramGroup.class
)

@DocumentedFeature
public class CalculateFingerprintMetrics extends CommandLineProgram {

    static final String USAGE_DETAILS =
            "This tools collects various statistics that pertain to a single fingerprint (<b>not</b> the comparison, or " +
                    "\'fingerprinting\' of two distinct samples) and reports the results in a metrics file. " +
                    "<p>" +
                    "The statistics collected are p-values, where the null-hypothesis is that the fingerprint is collected from " +
                    "a non-contaminated, diploid human, whose genotypes are modelled by the probabilities given in the " +
                    "HAPLOTYPE_MAP file." +
                    "<p>" +
                    "<h3>Example</h3>\n" +
                    "<pre>\" +\n" +
                    "java -jar picard.jar CalculateFingerprintMetrics \\\n" +
                    "      INPUT=sample.bam \\\n" +
                    "      HAPLOTYPE_DATABASE=fingerprinting_haplotype_database.txt \\\n" +
                    "      OUTPUT=sample.fingerprint_metrics\n" +
                    " </pre>\n";

    static final String USAGE_SUMMARY="Calculate statistics on fingerprints, checking their viability";

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "One or more input files (SAM/BAM/CRAM or VCF).")
    public List<String> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output files to write.")
    public File OUTPUT;

    @Argument(shortName = "H", doc = "The file lists a set of SNPs, optionally arranged in high-LD blocks, to be used for fingerprinting. See " +
            "https://software.broadinstitute.org/gatk/documentation/article?id=9526 for details.")
    public File HAPLOTYPE_MAP;

    @Argument(doc = "Specificies which data-type should be used as the basic unit. Fingerprints from readgroups can " +
            "be \"rolled-up\" to the LIBRARY, SAMPLE, or FILE level before being used." +
            " Fingerprints from VCF can be be examined by SAMPLE or FILE.")
    public CrosscheckMetric.DataType CHECK_BY = CrosscheckMetric.DataType.READGROUP;

    @Override
    protected int doWork() {

        final List<Path> inputPaths = IOUtil.getPaths(INPUT);
        IOUtil.assertPathsAreReadable(inputPaths);
        IOUtil.assertFileIsReadable(HAPLOTYPE_MAP);
        IOUtil.assertFileIsWritable(OUTPUT);

        final FingerprintChecker checker = new FingerprintChecker(HAPLOTYPE_MAP);

        final MetricsFile<FingerprintMetrics, ?> metricsFile = getMetricsFile();
        final Map<FingerprintIdDetails, Fingerprint> fpMap = checker.fingerprintFiles(inputPaths, 1, 1, TimeUnit.DAYS);

        final Map<FingerprintIdDetails, Fingerprint> mergedFpMap = Fingerprint.mergeFingerprintsBy(fpMap,Fingerprint.getFingerprintIdDetailsStringFunction(CHECK_BY));

        metricsFile.addAllMetrics(mergedFpMap.values().stream().map(Fingerprint::getFingerprintMetrics).collect(Collectors.toList()));
        metricsFile.write(OUTPUT);

        return 0;
    }
}
