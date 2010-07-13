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
import net.sf.samtools.*;
import net.sf.samtools.util.Log;

import java.io.File;

/**
 * Command line program to generate a BAM index (.bai) file from an existing BAM (.bam) file
 *
 * @author Martha Borkan
 */
public class BuildBamIndex extends CommandLineProgram {
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Generates a BAM index (.bai) file. " +
            "Input BAM file must be sorted in coordinate order. " +
            "Output file defaults to INPUT.bai (or if TEXTUAL, INPUT.bai.txt";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc="A BAM file to process.")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc="The BAM index file", optional=true)
    public File OUTPUT;

    @Option(doc = "Whether to overwrite an existing bai file", shortName = "OVERWRITE")
    public Boolean OVERWRITE_EXISTING_BAI_FILE = false;

    @Option(doc = "Whether to sort the bins in sequence order")
    public Boolean SORT = true;

    @Option(doc = "Whether to write textual or binary representation", shortName = "TEXT")
    public Boolean TEXTUAL = false;       // todo is there a way to hide this option

    /** Stock main method for a command line program. */
    public static void main(final String[] argv) {
        System.exit(new BuildBamIndex().instanceMain(argv));
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
            final SAMFileReader bam = new SAMFileReader(INPUT);
            bam.enableFileSource(true);

            // check sort order in the header - must be coordinate
            if (!bam.getFileHeader().getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
                throw new PicardException("Input BAM file must be sorted by coordinates");
            }

            final BAMIndexer instance = new BAMIndexer(INPUT, OUTPUT,
                    bam.getFileHeader().getSequenceDictionary().size(), SORT, TEXTUAL);

            instance.createIndex();

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
            String baseFileName = INPUT.getAbsolutePath();
            if (baseFileName.endsWith(".bam")){
                final int lastDot = baseFileName.lastIndexOf('.');
                if (lastDot != -1){
                   baseFileName = baseFileName.substring(0,lastDot);
                }
            }
            if (TEXTUAL){
                OUTPUT = new File (baseFileName + BAMIndex.BAMIndexSuffix + ".txt");
            } else {
                OUTPUT = new File (baseFileName + BAMIndex.BAMIndexSuffix);
            }
        }

        // make sure input is BAM file not SAM
        final SAMFileReader bam = new SAMFileReader(INPUT);
        if (!bam.isBinary()){
            return new String[] {"Input file must be bam file, not sam file"};
        }

        // check OVERWRITE
        if (!OVERWRITE_EXISTING_BAI_FILE && OUTPUT.exists()){
            return new String[] {"Output file already exists.  Use option OVERWRITE=true to replace " + OUTPUT.toString()};
        } else if (OVERWRITE_EXISTING_BAI_FILE && OUTPUT.exists()){
            OUTPUT.delete();
        }
        return null;
    }
}