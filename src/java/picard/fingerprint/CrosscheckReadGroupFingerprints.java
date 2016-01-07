
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

import htsjdk.samtools.BamFileIoUtils;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.FormatUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.SAMReadGroupRecord;
import picard.cmdline.programgroups.Alpha;

import java.io.File;
import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.TimeUnit;

/**
 * Program to check that all read groups within the set of BAM files appear to come from the same
 * individual.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage = "Checks if all read groups within a set of BAM files appear to come from the same individual",
        usageShort = "Checks if all read groups appear to come from the same individual",
        programGroup = Alpha.class  // TODO -- when mature please move to a to-be-created Fingerprinting.class
)
public class CrosscheckReadGroupFingerprints extends CommandLineProgram {

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc="One or more input BAM files (or lists of BAM files) to compare fingerprints for.")
    public List<File> INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional=true,
            doc="Optional output file to write metrics to. Default is to write to stdout.")
    public File OUTPUT;

    @Option(shortName="H", doc="The file of haplotype data to use to pick SNPs to fingerprint")
    public File HAPLOTYPE_MAP;

    @Option(shortName="LOD",
            doc="If any two read groups match with a LOD score lower than the threshold the program will exit " +
                "with a non-zero code to indicate error. 0 means equal probability the read groups match vs. " +
                "come from different individuals, negative numbers mean N logs more likely that the read groups " +
                "are from different individuals and positive numbers mean N logs more likely that the read groups " +
                "are from the sample individual.")
    public double LOD_THRESHOLD = 0;

	@Option(doc="Instead of producing the normal comparison of read-groups, roll fingerprints up to the sample level " +
            "and print out a sample x sample matrix with LOD scores.")
	public boolean CROSSCHECK_SAMPLES = false;

    @Option(doc="Instead of producing the normal comparison of read-groups, roll fingerprints up to the library level " +
            "and print out a library x library matrix with LOD scores.")
    public boolean CROSSCHECK_LIBRARIES = false;

	@Option(doc="The number of threads to use to process BAM files and generate Fingerprints.")
	public int NUM_THREADS = 1;

    @Option(doc="Allow the use of duplicate reads in performing the comparison. Can be useful when duplicate " +
            "marking has been overly aggressive and coverage is low.")
    public boolean ALLOW_DUPLICATE_READS = false;

    @Option(doc="Assumed genotyping error rate that provides a floor on the probability that a genotype comes from" +
            " the expected sample.")
    public double GENOTYPING_ERROR_RATE = 0.01;

    @Option(doc="If true then only read groups that do not relate to each other as expected will have their LODs reported.")
    public boolean OUTPUT_ERRORS_ONLY = false;

    @Option(doc ="The rate at which a het in a normal sample turns into a hom in the tumor.", optional = true)
    public double LOSS_OF_HET_RATE = 0.5;

    @Option(doc="Expect all read groups' fingerprints to match, irrespective of their sample names.  By default (with this value set to " +
            "false), read groups with different sample names are expected to mismatch, and those with the same sample name are expected " +
            "to match.")
    public boolean EXPECT_ALL_READ_GROUPS_TO_MATCH = false;

    @Option(doc="When one or more mismatches between read groups are detected, exit with this value instead of 0.")
    public int EXIT_CODE_WHEN_MISMATCH = 1;

    private final Log log = Log.getInstance(CrosscheckReadGroupFingerprints.class);

    private final FormatUtil formatUtil = new FormatUtil();

    // These are public so that other programs can parse status from the crosscheck file
    public static final String EXPECTED_MATCH = "EXPECTED MATCH";
    public static final String EXPECTED_MISMATCH = "EXPECTED MISMATCH";
    public static final String UNEXPECTED_MATCH = "UNEXPECTED MATCH";
    public static final String UNEXPECTED_MISMATCH = "UNEXPECTED MISMATCH";

    /** Stock main method. */
    public static void main(final String[] args) {
        new CrosscheckReadGroupFingerprints().instanceMainWithExit(args);
    }

    @Override protected int doWork() {
        // Check inputs
        for (final File f : INPUT) IOUtil.assertFileIsReadable(f);
        IOUtil.assertFileIsReadable(HAPLOTYPE_MAP);
        if (OUTPUT != null) IOUtil.assertFileIsWritable(OUTPUT);

        final HaplotypeMap map = new HaplotypeMap(HAPLOTYPE_MAP);
        final FingerprintChecker checker = new FingerprintChecker(map);

        checker.setAllowDuplicateReads(ALLOW_DUPLICATE_READS);

		log.info("Done checking input files, moving onto fingerprinting files.");

        List<File> unrolledFiles = IOUtil.unrollFiles(INPUT, BamFileIoUtils.BAM_FILE_EXTENSION, IOUtil.SAM_FILE_EXTENSION);
        final Map<SAMReadGroupRecord, Fingerprint> fpMap = checker.fingerprintSamFiles(unrolledFiles, NUM_THREADS, 1, TimeUnit.DAYS);
        final List<Fingerprint> fingerprints = new ArrayList<>(fpMap.values());

		log.info("Finished generating fingerprints from BAM files, moving on to cross-checking.");

        // Setup the output
        final PrintStream out;
        if (OUTPUT != null) {
            out = new PrintStream(IOUtil.openFileForWriting(OUTPUT), true);
        }
        else {
            out = System.out;
        }

		if (this.CROSSCHECK_SAMPLES) {
			crossCheckSamples(fingerprints, out);
			return 0;
		}
        else if (this.CROSSCHECK_LIBRARIES) {
            crossCheckLibraries(fpMap, out);
            return 0;
        }
		else {
			return crossCheckReadGroups(fpMap, out);
		}
    }

	/**
	 * Method that combines the fingerprint evidence across all the read groups for the same sample
	 * and then produces a matrix of LOD scores for comparing every sample with every other sample.
	 */
	private void crossCheckSamples(final List<Fingerprint> fingerprints, final PrintStream out) {
		final SortedMap<String,Fingerprint> sampleFps = FingerprintChecker.mergeFingerprintsBySample(fingerprints);
		final SortedSet<String> samples = (SortedSet<String>) sampleFps.keySet();

		// Print header row
		out.print("\t");
		for (final String sample : samples) { out.print(sample); out.print("\t"); }
		out.println();

		// Print results rows
		for (final String sample : samples) {
			out.print(sample);
			final Fingerprint fp = sampleFps.get(sample);

			for (final String otherSample : samples) {
				final MatchResults results = FingerprintChecker.calculateMatchResults(fp, sampleFps.get(otherSample), GENOTYPING_ERROR_RATE, LOSS_OF_HET_RATE);
				out.print("\t");
				out.print(formatUtil.format(results.getLOD()));
			}

			out.println();
		}
	}

    /**
     * Method that combines the fingerprint evidence across all the read groups for the same library
     * and then produces a matrix of LOD scores for comparing every library with every other library.
     */
    private void crossCheckLibraries(final Map<SAMReadGroupRecord,Fingerprint> fingerprints, final PrintStream out) {
        final List<Fingerprint> fixedFps = new ArrayList<>();
        for (final SAMReadGroupRecord rg : fingerprints.keySet()) {
            final Fingerprint old = fingerprints.get(rg);
            final String name = rg.getSample() + "::" + rg.getLibrary();
            final Fingerprint newFp = new Fingerprint(name, old.getSource(), old.getInfo());
            newFp.putAll(old);

            fixedFps.add(newFp);
        }

        crossCheckSamples(fixedFps, out);
    }

	/**
	 * Method that pairwise checks every pair of read groups and reports a LOD score for the two read groups
	 * coming from the same sample.
	 */
	private int crossCheckReadGroups(final Map<SAMReadGroupRecord,Fingerprint> fingerprints, final PrintStream out) {
		int mismatches = 0;
		int unexpectedMatches = 0;

		final List<SAMReadGroupRecord> readGroupRecords = new ArrayList<>(fingerprints.keySet());
		final List<String> output = new ArrayList<>();

		for (int i = 0; i < readGroupRecords.size(); i++) {
			final SAMReadGroupRecord lhsRg = readGroupRecords.get(i);
			for (int j= i+1; j < readGroupRecords.size(); j++) {
				final SAMReadGroupRecord rhsRg = readGroupRecords.get(j);
				final boolean expectedToMatch = EXPECT_ALL_READ_GROUPS_TO_MATCH || lhsRg.getSample().equals(rhsRg.getSample());

				final MatchResults results = FingerprintChecker.calculateMatchResults(fingerprints.get(lhsRg), fingerprints.get(rhsRg), GENOTYPING_ERROR_RATE, LOSS_OF_HET_RATE);
                if (expectedToMatch) {
                    if (results.getLOD() < LOD_THRESHOLD) {
                        mismatches++;
                        output.add(getMatchDetails(UNEXPECTED_MISMATCH, results, lhsRg, rhsRg));
                    } else {
                        if (!OUTPUT_ERRORS_ONLY) {
                            output.add(getMatchDetails(EXPECTED_MATCH, results, lhsRg, rhsRg));
                        }
                    }
                } else {
                    if (results.getLOD() > -LOD_THRESHOLD) {
                        unexpectedMatches++;
                        output.add(getMatchDetails(UNEXPECTED_MATCH, results, lhsRg, rhsRg));
                    } else {
                        if (!OUTPUT_ERRORS_ONLY) {
                            output.add(getMatchDetails(EXPECTED_MISMATCH, results, lhsRg, rhsRg));
                        }
                    }
                }
			}
		}

		if (!output.isEmpty()) {
			out.println("RESULT\tLOD_SCORE\tLOD_SCORE_TUMOR_NORMAL\tLOD_SCORE_NORMAL_TUMOR\tLEFT_RUN_BARCODE\tLEFT_LANE\tLEFT_MOLECULAR_BARCODE_SEQUENCE\tLEFT_LIBRARY\tLEFT_SAMPLE\t" +
				"RIGHT_RUN_BARCODE\tRIGHT_LANE\tRIGHT_MOLECULAR_BARCODE_SEQUENCE\tRIGHT_LIBRARY\tRIGHT_SAMPLE");
			out.println(String.join("\n", output));
		}

		if (mismatches + unexpectedMatches > 0) {
			log.info("WARNING: At least two read groups did not relate as expected.");
			return EXIT_CODE_WHEN_MISMATCH;
		}
		else {
			log.info("All read groups related as expected.");
			return 0;
		}
	}

    /**
     * Generates tab delimited string containing details about a possible match between fingerprints on two different SAMReadGroupRecords
     * @param matchResult String describing the match type.
     * @param results MatchResults object
     * @param left left hand side SAMReadGroupRecord
     * @param right right hand side SAMReadGroupRecord
     * @return tab delimited string containing details about a possible match
     */
    private String getMatchDetails(final String matchResult, final MatchResults results, final SAMReadGroupRecord left, final SAMReadGroupRecord right) {
        final List<String> elements = new ArrayList<>(4);
        elements.add(matchResult);
        elements.add(formatUtil.format(results.getLOD()));
        elements.add(formatUtil.format(results.getLodTN()));
        elements.add(formatUtil.format(results.getLodNT()));
        elements.add(getReadGroupDetails(left));
        elements.add(getReadGroupDetails(right));
        return String.join("\t", elements);
    }

    /**
     * Generates tab delimited string containing details about the passed SAMReadGroupRecord
     * @param readGroupRecord record
     * @return tab delimited string containing details about the SAMReadGroupRecord
     */
    private String getReadGroupDetails(final SAMReadGroupRecord readGroupRecord) {
        final List<String> elements = new ArrayList<>(5);

        final String tmp[] = readGroupRecord.getPlatformUnit().split("\\.");    // Expect to look like: D047KACXX110901.1.ACCAACTG
        String runBarcode = "?";
        String lane = "?";
        String molBarcode = "?";
        if ((tmp.length == 3) || (tmp.length == 2)) {
            runBarcode = tmp[0];
            lane = tmp[1];
            molBarcode = (tmp.length == 3) ? tmp[2] : "";       // In older BAMS there may be no molecular barcode sequence
        } else {
            log.error("Unexpected format " + readGroupRecord.getPlatformUnit() + " for PU attribute");
        }
        elements.add(runBarcode);
        elements.add(lane);
        elements.add(molBarcode);
        elements.add(readGroupRecord.getLibrary());
        elements.add(readGroupRecord.getSample());
        return String.join("\t", elements);
    }
}
