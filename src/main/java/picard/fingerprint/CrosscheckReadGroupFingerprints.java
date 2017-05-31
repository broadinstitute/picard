
/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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

import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Fingerprinting;

import java.io.File;
import java.util.List;

/**
 * Program to check that all read groups within the set of BAM files appear to come from the same
 * individual.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage = "DEPRECATED: USE CrosscheckFingerprints. Checks if all read groups within a set of BAM files appear to come from the same individual",
        usageShort = "DEPRECATED: USE CrosscheckFingerprints. Checks if all read groups appear to come from the same individual",
        programGroup = Fingerprinting.class
)
@Deprecated // USE CrosscheckFingerprints instead.
public class CrosscheckReadGroupFingerprints extends CommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "One or more input BAM files (or lists of BAM files) to compare fingerprints for.")
    public List<File> INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional = true,
            doc = "Optional output file to write metrics to. Default is to write to stdout.")
    public File OUTPUT;

    @Option(shortName = "H", doc = "The file lists a set of SNPs, optionally arranged in high-LD blocks, to be used for fingerprinting. See " +
            "https://software.broadinstitute.org/gatk/documentation/article?id=9526 for details.")
    public File HAPLOTYPE_MAP;

    @Option(shortName = "LOD",
            doc = "If any two read groups match with a LOD score lower than the threshold the program will exit " +
                    "with a non-zero code to indicate error. 0 means equal probability the read groups match vs. " +
                    "come from different individuals, negative numbers mean N logs more likely that the read groups " +
                    "are from different individuals and positive numbers mean N logs more likely that the read groups " +
                    "are from the sample individual.")
    public double LOD_THRESHOLD = 0;

    @Option(doc = "Instead of producing the normal comparison of read-groups, roll fingerprints up to the sample level " +
            "and print out a sample x sample matrix with LOD scores.")
    public boolean CROSSCHECK_SAMPLES = false;

    @Option(doc = "Instead of producing the normal comparison of read-groups, roll fingerprints up to the library level " +
            "and print out a library x library matrix with LOD scores.")
    public boolean CROSSCHECK_LIBRARIES = false;

    @Option(doc = "The number of threads to use to process BAM files and generate Fingerprints.")
    public int NUM_THREADS = 1;

    @Option(doc = "Allow the use of duplicate reads in performing the comparison. Can be useful when duplicate " +
            "marking has been overly aggressive and coverage is low.")
    public boolean ALLOW_DUPLICATE_READS = false;

    @Option(doc = "Assumed genotyping error rate that provides a floor on the probability that a genotype comes from" +
            " the expected sample.")
    public double GENOTYPING_ERROR_RATE = 0.01;

    @Option(doc = "If true then only read groups that do not relate to each other as expected will have their LODs reported.")
    public boolean OUTPUT_ERRORS_ONLY = false;

    @Option(doc = "The rate at which a het in a normal sample turns into a hom in the tumor.", optional = true)
    public double LOSS_OF_HET_RATE = 0.5;

    @Option(doc = "Expect all read groups' fingerprints to match, irrespective of their sample names.  By default (with this value set to " +
            "false), read groups with different sample names are expected to mismatch, and those with the same sample name are expected " +
            "to match.")
    public boolean EXPECT_ALL_READ_GROUPS_TO_MATCH = false;

    @Option(doc = "When one or more mismatches between read groups are detected, exit with this value instead of 0.")
    public int EXIT_CODE_WHEN_MISMATCH = 1;

    @Override
    protected int doWork() {

        CrosscheckFingerprints crosscheckFingerprints = new CrosscheckFingerprints();

        if (this.CROSSCHECK_LIBRARIES) crosscheckFingerprints.CROSSCHECK_BY = CrosscheckMetric.DataType.LIBRARY;
        else if (this.CROSSCHECK_SAMPLES) crosscheckFingerprints.CROSSCHECK_BY = CrosscheckMetric.DataType.SAMPLE;
        else crosscheckFingerprints.CROSSCHECK_BY = CrosscheckMetric.DataType.READGROUP;

        switch (crosscheckFingerprints.CROSSCHECK_BY) {
            case LIBRARY:
            case SAMPLE:
                crosscheckFingerprints.OUTPUT = this.OUTPUT;
                break;
            case READGROUP:
                crosscheckFingerprints.MATRIX_OUTPUT = this.OUTPUT;
            default:
                throw new PicardException("Unpossible!");
        }

        crosscheckFingerprints.ALLOW_DUPLICATE_READS = this.ALLOW_DUPLICATE_READS;
        crosscheckFingerprints.GENOTYPING_ERROR_RATE = this.GENOTYPING_ERROR_RATE;
        crosscheckFingerprints.EXIT_CODE_WHEN_MISMATCH = this.EXIT_CODE_WHEN_MISMATCH;
        crosscheckFingerprints.EXPECT_ALL_GROUPS_TO_MATCH = this.EXPECT_ALL_READ_GROUPS_TO_MATCH;
        crosscheckFingerprints.HAPLOTYPE_MAP = this.HAPLOTYPE_MAP;
        crosscheckFingerprints.INPUT = this.INPUT;
        crosscheckFingerprints.LOD_THRESHOLD = this.LOD_THRESHOLD;
        crosscheckFingerprints.LOSS_OF_HET_RATE = this.LOSS_OF_HET_RATE;
        crosscheckFingerprints.NUM_THREADS = this.NUM_THREADS;
        crosscheckFingerprints.OUTPUT_ERRORS_ONLY = this.OUTPUT_ERRORS_ONLY;
        crosscheckFingerprints.GA4GH_CLIENT_SECRETS = this.GA4GH_CLIENT_SECRETS;
        crosscheckFingerprints.REFERENCE_SEQUENCE = this.REFERENCE_SEQUENCE;
        crosscheckFingerprints.VERBOSITY = this.VERBOSITY;
        crosscheckFingerprints.QUIET = this.QUIET;
        crosscheckFingerprints.MAX_RECORDS_IN_RAM = this.MAX_RECORDS_IN_RAM;
        crosscheckFingerprints.TMP_DIR = this.TMP_DIR;
        crosscheckFingerprints.VALIDATION_STRINGENCY = this.VALIDATION_STRINGENCY;
        crosscheckFingerprints.CREATE_MD5_FILE = this.CREATE_MD5_FILE;
        crosscheckFingerprints.CREATE_INDEX = this.CREATE_INDEX;

        return crosscheckFingerprints.doWork();
    }
}
     