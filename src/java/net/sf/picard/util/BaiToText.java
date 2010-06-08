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

package net.sf.picard.util;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.samtools.BAMFileIndexWriter;
import net.sf.samtools.BAMIndexTextWriter;

import java.io.File;

/**
 * Command line program to generate textual form of a BAM index (.bai) file
 *
 * @author Martha Borkan
 */
public class BaiToText extends CommandLineProgram {
    @Usage(programVersion="1.1")
    public String USAGE = getStandardUsagePreamble() + "Generates a textual form of a BAM index (.bai) file.";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc="A BAM Index file to process.")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc="The textual BAM index file output", optional=true)
    public File OUTPUT;

    @Option(doc = "Whether to overwrite an existing bai.txt file", shortName = "OVERWRITE")
    public Boolean OVERWRITE_EXISTING_BAI_FILE = false;

    @Option(doc = "Whether to sort the bins in sequence order")
    public Boolean SORT = true;

    @Option(doc = "Whether to write textual or binary representation", shortName = "TEXT")
    public Boolean TEXTUAL = true;  // binary output is used to test. Should generate the same file as input

    /** Stock main method for a command line program. */
    public static void main(final String[] argv) {
        System.exit(new BaiToText().instanceMain(argv));
    }

    /**
     * Main method for the program.  Checks that all input files are present and
     * readable and that the output file can be written to.  Then iterates through
     * all the records generating a BAM Index, then writes the bai file.
     */
    protected int doWork() {
        final Log log = Log.getInstance(getClass());

        // Some quick parameter checking
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);

        log.info("Reading input file and building index.");

        try {

            final BAMIndexTextWriter instance = new BAMIndexTextWriter(INPUT, OUTPUT);

            if (TEXTUAL){
                instance.writeText(SORT);
            } else {
                instance.writeBinary(SORT, 0);
            }

        } catch (Exception e) {
            log.error(e.getMessage() + " exception when writing output file " + OUTPUT + e);
            e.printStackTrace();
            return 1;
        }
        log.info("Successfully wrote bam index file " + OUTPUT);

        return 0;
    }

    @Override
    protected String[] customCommandLineValidation() {
        // set default output file - input-file.bai
        if (OUTPUT == null){
            OUTPUT = new File (INPUT.getName() + ".txt");
        }
        // check OVERWRITE
        if (!OVERWRITE_EXISTING_BAI_FILE && OUTPUT.exists()){
            return new String[] {"Output file already exists.  Use option OVERWRITE=true to replace " + OUTPUT.toString()};
        }
        return null;
    }
}