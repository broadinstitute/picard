/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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

package picard.sam;

import htsjdk.samtools.SAMValidationError;
import htsjdk.samtools.SamFileValidator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.BamIndexValidator.IndexValidationStringency;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import htsjdk.samtools.util.Log;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * Command line program wrapping SamFileValidator.
 *
 * @author Doug Voet
 */
@CommandLineProgramProperties(
        summary = ValidateSamFile.USAGE_SUMMARY + ValidateSamFile.USAGE_DETAILS,
        oneLineSummary = ValidateSamFile.USAGE_SUMMARY,
        programGroup = SamOrBam.class
)
public class ValidateSamFile extends CommandLineProgram {
    private static final Log log = Log.getInstance(ValidateSamFile.class);
    static final String USAGE_SUMMARY = "Validates a SAM or BAM file.  ";
    static final String USAGE_DETAILS = "<p>This tool reports on the validity of a SAM or BAM file relative to the SAM format " +
            "specification.  This is useful for troubleshooting errors encountered with other tools that may be caused by improper " +
            "formatting, faulty alignments, incorrect flag values, etc. </p> " +

            "<p>By default, the tool runs in VERBOSE mode and will exit after finding 100 errors and output them to the console (stdout). " +
            "Therefore, it is often more practical to run this tool initially using the MODE=SUMMARY option.  This mode outputs a summary " +
            "table listing the numbers of all 'errors' and 'warnings'.</p> "+

            "<p>When fixing errors in your file, it is often useful to prioritize the severe validation errors and ignore the " +
            "errors/warnings of lesser concern.  This can be done using the IGNORE and/or IGNORE_WARNINGS arguments.  For helpful " +
            "suggestions on error prioritization, please follow this link to obtain additional documentation on <a href='https://www.broadinstitute.org/gatk/guide/article?id=7571'>ValidateSamFile</a>.</p>" +

            "<p>After identifying and fixing your 'warnings/errors', we recommend that you rerun this tool to validate your SAM/BAM " +
            "file prior to proceeding with your downstream analysis.  This will verify that all problems in your file have been addressed.</p>" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar ValidateSamFile \\<br />" +
            "      I=input.bam \\<br />" +
            "      MODE=SUMMARY" +
            "</pre>" +
            "<p>To obtain a complete list with descriptions of both 'ERROR' and 'WARNING' messages, please see our additional " +
            " <a href='https://www.broadinstitute.org/gatk/guide/article?id=7571'>documentation</a> for this tool.</p>" +
            ""+
            "<hr />";
    public enum Mode {VERBOSE, SUMMARY}

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "Input SAM/BAM file")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file or standard out if missing",
            optional = true)
    public File OUTPUT;

    @Argument(shortName = "M",
            doc = "Mode of output")
    public Mode MODE = Mode.VERBOSE;

    @Argument(doc = "List of validation error types to ignore.", optional = true)
    public List<SAMValidationError.Type> IGNORE = new ArrayList<SAMValidationError.Type>();

    @Argument(shortName = "MO",
            doc = "The maximum number of lines output in verbose mode")
    public Integer MAX_OUTPUT = 100;

    @Argument(doc = "If true, only report errors and ignore warnings.")
    public boolean IGNORE_WARNINGS = false;

    @Argument(doc = "DEPRECATED.  Use INDEX_VALIDATION_STRINGENCY instead.  If true and input is " +
            "a BAM file with an index file, also validates the index.  Until this parameter is retired " +
            "VALIDATE INDEX and INDEX_VALIDATION_STRINGENCY must agree on whether to validate the index.")
    public boolean VALIDATE_INDEX = true;

    @Argument(doc = "If set to anything other than IndexValidationStringency.NONE and input is " +
            "a BAM file with an index file, also validates the index at the specified stringency. " +
            "Until VALIDATE_INDEX is retired, VALIDATE INDEX and INDEX_VALIDATION_STRINGENCY " +
            "must agree on whether to validate the index.")
    public IndexValidationStringency INDEX_VALIDATION_STRINGENCY = IndexValidationStringency.EXHAUSTIVE;

    @Argument(shortName = "BISULFITE",
            doc = "Whether the SAM or BAM file consists of bisulfite sequenced reads. " +
                    "If so, C->T is not counted as an error in computing the value of the NM tag.")
    public boolean IS_BISULFITE_SEQUENCED = false;

    @Argument(doc = "Relevant for a coordinate-sorted file containing read pairs only. " +
            "Maximum number of file handles to keep open when spilling mate info to disk. " +
            "Set this number a little lower than the per-process maximum number of file that may be open. " +
            "This number can be found by executing the 'ulimit -n' command on a Unix system.")
    public int MAX_OPEN_TEMP_FILES = 8000;

    public static void main(final String[] args) {
        System.exit(new ValidateSamFile().instanceMain(args));
    }

    /**
     * Return types for doWork()
     */
    enum ReturnTypes {
        FAILED(-1),     // failed to complete execution
        SUCCESSFUL(0),  // ran successfully
        WARNINGS(1),    // warnings but no errors
        ERRORS_WARNINGS(2),      // errors and warnings
        ERRORS(3);      // errors but no warnings

        private final int value;

        ReturnTypes(final int value) {
            this.value = value;
        }

        int value() {
            return value;
        }
    }

    /**
     * Do the work after command line has been parsed.
     *
     * @return -1 - failed to complete execution, 0 - ran successfully without warnings or errors, 1 - warnings but no errors,
     * 2 - errors and warnings, 3 - errors but no warnings
     *
     */
    @Override
    protected int doWork() {
        try {
            IOUtil.assertFileIsReadable(INPUT);
            ReferenceSequenceFile reference = null;
            if (REFERENCE_SEQUENCE != null) {
                IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
                reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);

            }
            final PrintWriter out;
            if (OUTPUT != null) {
                IOUtil.assertFileIsWritable(OUTPUT);
                try {
                    out = new PrintWriter(OUTPUT);
                } catch (FileNotFoundException e) {
                    // we already asserted this so we should not get here
                    throw new PicardException("Unexpected exception", e);
                }
            } else {
                out = new PrintWriter(System.out);
            }

            boolean result;

            final SamReaderFactory factory = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE)
                    .validationStringency(ValidationStringency.SILENT)
                    .enable(SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS);
            final SamReader samReader = factory.open(INPUT);

            if (samReader.type() != SamReader.Type.BAM_TYPE) VALIDATE_INDEX = false;

            factory.setOption(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES, VALIDATE_INDEX);
            factory.reapplyOptions(samReader);

            final SamFileValidator validator = new SamFileValidator(out, MAX_OPEN_TEMP_FILES);
            validator.setErrorsToIgnore(IGNORE);

            if (IGNORE_WARNINGS) {
                validator.setIgnoreWarnings(IGNORE_WARNINGS);
            }
            if (MODE == Mode.SUMMARY) {
                validator.setVerbose(false, 0);
            } else {
                validator.setVerbose(true, MAX_OUTPUT);
            }
            if (IS_BISULFITE_SEQUENCED) {
                validator.setBisulfiteSequenced(IS_BISULFITE_SEQUENCED);
            }
            if (VALIDATE_INDEX) {
                validator.setIndexValidationStringency(VALIDATE_INDEX ? IndexValidationStringency.EXHAUSTIVE : IndexValidationStringency.NONE);
            }
            if (IOUtil.isRegularPath(INPUT)) {
                // Do not check termination if reading from a stream
                validator.validateBamFileTermination(INPUT);
            }

            result = false;

            switch (MODE) {
                case SUMMARY:
                    result = validator.validateSamFileSummary(samReader, reference);
                    break;
                case VERBOSE:
                    result = validator.validateSamFileVerbose(samReader, reference);
                    break;
            }
            out.flush();

            if (result) {
                return ReturnTypes.SUCCESSFUL.value();  // ran successfully with no warnings or errors
            } else {
                if (validator.getNumErrors() == 0) {
                    if (validator.getNumWarnings() > 0) {
                        return ReturnTypes.WARNINGS.value();   // warnings but no errors
                    } else {
                        log.error("SAM file validation fails without warnings or errors.");
                        return ReturnTypes.FAILED.value();
                    }
                } else {
                    if (validator.getNumWarnings() > 0) {
                        return ReturnTypes.ERRORS_WARNINGS.value();   // errors and warnings
                    } else {
                        return ReturnTypes.ERRORS.value();  // errors but no warnings
                    }
                }
            }
        } catch (Exception e) {
            log.error(e.getMessage());
            return ReturnTypes.FAILED.value();  // failed to complete execution
        }
    }

    @Override
    protected String[] customCommandLineValidation() {
        if ((!VALIDATE_INDEX && INDEX_VALIDATION_STRINGENCY != IndexValidationStringency.NONE) ||
            (VALIDATE_INDEX && INDEX_VALIDATION_STRINGENCY == IndexValidationStringency.NONE)) {
            return new String[]{"VALIDATE_INDEX and INDEX_VALIDATION_STRINGENCY must be consistent: " +
                    "VALIDATE_INDEX is " + VALIDATE_INDEX + " and INDEX_VALIDATION_STRINGENCY is " +
                    INDEX_VALIDATION_STRINGENCY};
        }

        return super.customCommandLineValidation();
    }
}
