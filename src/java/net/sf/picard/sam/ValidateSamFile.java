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

package net.sf.picard.sam;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.List;
import java.util.ArrayList;

import net.sf.picard.cmdline.Usage;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMValidationError;
import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.io.IoUtil;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;

/**
 * Commandline program wrapping SamFileValidator.
 *
 * @author Doug Voet
 */
public class ValidateSamFile extends CommandLineProgram {
    @Usage
    public final String usage = getStandardUsagePreamble() + "Read a SAM or BAM file and report on its validity.";
    
    public enum Mode { VERBOSE, SUMMARY }
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input SAM/BAM file")
    public File INPUT;
    
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file or standard out if missing", optional=true)
    public File OUTPUT;
    
    @Option(shortName="M", doc="Mode of output")
    public Mode MODE = Mode.VERBOSE;

    @Option(doc="List of validation error types to ignore.")
    public List<SAMValidationError.Type> IGNORE = new ArrayList<SAMValidationError.Type>();
    
    @Option(shortName="MO", doc="The maximum number of lines output in verbose mode")
    public Integer MAX_OUTPUT = 100;
    
    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, 
            doc="Reference sequence file, the NM tag check will be skipped if this is missing", 
            optional=true)
    public File REFERENCE_SEQUENCE;

    @Option(doc="If true, only report errors, and ignore warnings.")
    public boolean IGNORE_WARNINGS = false;

    @Option(doc="If true and input is a BAM file with an index file, also validates the index.")
    public boolean VALIDATE_INDEX = true;

    @Option (shortName="BISULFITE", doc="Whether the SAM or BAM file consists of bisulfite sequenced reads. " +
            "If so, C->T is not counted as an error in computing the value of the NM tag.")
    public boolean IS_BISULFITE_SEQUENCED = false;

    
    @Option(doc="Relevant for a coordinate-sorted file containing read pairs only.  " +
            "Maximum number of file handles to keep open when spilling mate info to disk.  " + "" +
            "Set this number a little lower than the per-process maximum number of file that may be open.  " +
            "This number can be found by executing the 'ulimit -n' command on a Unix system.")
    public int MAX_OPEN_TEMP_FILES = 8000;

    public static void main(final String[] args) {
        System.exit(new ValidateSamFile().instanceMain(args));
    }

    @Override
    protected int doWork() {
        IoUtil.assertFileIsReadable(INPUT);
        ReferenceSequenceFile reference = null;
        if (REFERENCE_SEQUENCE != null) {
            IoUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
            reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);

        }
        final PrintWriter out;
        if (OUTPUT != null) {
            IoUtil.assertFileIsWritable(OUTPUT);
            try {
                out = new PrintWriter(OUTPUT);
            }
            catch (FileNotFoundException e) {
                // we already asserted this so we should not get here
                throw new PicardException("Unexpected exception", e);
            }
        }
        else {
            out = new PrintWriter(System.out);
        }

        final SAMFileReader.ValidationStringency originalStringency = SAMFileReader.getDefaultValidationStringency();
        SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        boolean result;
        try {
            final SAMFileReader samReader = new SAMFileReader(INPUT);
            if (!samReader.isBinary()) VALIDATE_INDEX=false;
            if (VALIDATE_INDEX){
                samReader.enableIndexCaching(true);
            }
            samReader.enableCrcChecking(true);

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
            if (VALIDATE_INDEX){
                validator.setValidateIndex(VALIDATE_INDEX);
            }
            validator.validateBamFileTermination(INPUT);

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
        } finally {
            SAMFileReader.setDefaultValidationStringency(originalStringency);
        }

        return result ? 0 : 1;
    }
}
