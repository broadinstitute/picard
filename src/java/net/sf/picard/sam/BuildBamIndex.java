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

package net.sf.picard.sam;

import net.sf.picard.cmdline.*;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.*;

import java.io.File;
import java.net.URL;

/**
 * Command line program to generate a BAM index (.bai) file from a BAM (.bam) file
 *
 * @author Martha Borkan
 */
public class BuildBamIndex extends CommandLineProgram {

    private static final Log log = Log.getInstance(BAMIndexer.class);

    @Usage
    public String USAGE = getStandardUsagePreamble() + "Generates a BAM index (.bai) file. " +
            "Input BAM file must be sorted in coordinate order. " +
            "Output file defaults to x.bai if INPUT is x.bam, otherwise INPUT.bai)";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc="A BAM file or URL to process.")
    public String INPUT;

    URL inputUrl = null;   // INPUT as URL
    File inputFile = null; // INPUT as File, if it can't be interpreted as a valid URL

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc="The BAM index file.", optional=true)
    public File OUTPUT;

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

        try {
            inputUrl = new URL(INPUT);
        } catch (java.net.MalformedURLException e) {
            inputFile = new File(INPUT);
        }

        // set default output file - input-file.bai
        if (OUTPUT == null) {

            final String baseFileName;
            if (inputUrl != null) {
                String path = inputUrl.getPath();
                int lastSlash = path.lastIndexOf("/");
                baseFileName = path.substring(lastSlash + 1, path.length());
            } else {
                baseFileName = inputFile.getAbsolutePath();
            }

            if (baseFileName.endsWith(".bam")) {

                final int index = baseFileName.lastIndexOf(".");
                OUTPUT = new File(baseFileName.substring(0, index) + BAMIndex.BAMIndexSuffix);

            } else {
                OUTPUT = new File(baseFileName + BAMIndex.BAMIndexSuffix);
            }
        }

        IoUtil.assertFileIsWritable(OUTPUT);
        final SAMFileReader bam;

        if (inputUrl != null) {
            // remote input
            bam = new SAMFileReader(inputUrl, null, false);
        } else {
            // input from a normal file
            IoUtil.assertFileIsReadable(inputFile);
            bam = new SAMFileReader(inputFile);
        }

        if (!bam.isBinary()) {
            throw new SAMException("Input file must be bam file, not sam file.");
        }

        if (!bam.getFileHeader().getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
            throw new SAMException("Input bam file must be sorted by coordinates");
        }

        bam.enableFileSource(true);

        createIndex(bam, OUTPUT);

        log.info("Successfully wrote bam index file " + OUTPUT);
        return 0;
    }

    /**
     * Generates a BAM index file from an input BAM file
     *
     * @param reader SAMFileReader for input BAM file
     * @para output  File for output index file
     */
    public static void createIndex(SAMFileReader reader, File output) {

        BAMIndexer indexer = new BAMIndexer(output, reader.getFileHeader());

        reader.enableFileSource(true);
        int totalRecords = 0;

        // create and write the content
        for (SAMRecord rec : reader) {
            if (++totalRecords % 1000000 == 0) {
                log.info(totalRecords + " reads processed ...");
            }
            indexer.processAlignment(rec);
        }
        indexer.finish();
    }
}