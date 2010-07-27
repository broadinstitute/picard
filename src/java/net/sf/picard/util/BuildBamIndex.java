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

import net.sf.picard.cmdline.*;
import net.sf.picard.io.IoUtil;
import net.sf.samtools.*;

import java.io.File;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.List;

/**
 * Command line program to generate a BAM index (.bai) file from an existing BAM (.bam) file
 *
 * @author Martha Borkan
 */
public class BuildBamIndex extends CommandLineProgram {
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Generates a BAM index (.bai) file. " +
            "Input BAM file must be sorted in coordinate order. " +
            "Output file defaults to INPUT.bai (or x.bai if INPUT is x.bam)\n";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc="A BAM file to process.", optional=true, mutex="INPUT_URL")
    public File INPUT;

    @Option(shortName= "URL",
            doc="A BAM file to process.", optional=true, mutex="INPUT")
    public String INPUT_URL;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc="The BAM index file.", optional=true)
    public File OUTPUT;

    @Option(doc = "Whether to sort the bins in bin number order.")
    public Boolean SORT = false;

    @PositionalArguments(minElements = 0, maxElements = 2)
    public List<File> IN_OUT;

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

        IoUtil.assertFileIsWritable(OUTPUT);
        final SAMFileReader bam;

        if (INPUT_URL != null) {
            // remote input
            final URL bamURL;
            try {
                bamURL = new URL(INPUT_URL);
            } catch (MalformedURLException e) {
                throw new SAMException(e);
            }
            bam = new SAMFileReader(bamURL, null, false);

        } else {
            // input from a normal file
            IoUtil.assertFileIsReadable(INPUT);
            bam = new SAMFileReader(INPUT);
        }

        if (!bam.isBinary()){
            throw new SAMException ("Input file must be bam file, not sam file.");
        }

        if (!bam.getFileHeader().getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
            throw new SAMException ("Input bam file must be sorted by coordinates");
        }

        bam.enableFileSource(true);

        BAMIndexer instance = new BAMIndexer(OUTPUT,
                    bam.getFileHeader().getSequenceDictionary().size(), SORT);
        instance.createIndex(bam);

        log.info("Successfully wrote bam index file " + OUTPUT);
        return 0;
    }

    @Override
    protected String[] customCommandLineValidation() {
        // allow positional as alternative to i=input o=output
        if (IN_OUT.size() > 0) {
            if (INPUT != null) {
                return new String[]{"Can't specify both named and positional INPUT parameter"};
            }
            INPUT = IN_OUT.get(0);
            if (IN_OUT.size() > 1) {
                if (OUTPUT != null) {
                    return new String[]{"Can't specify both named and positional OUTPUT parameter"};
                }
                OUTPUT = IN_OUT.get(1);
            }
        }
        if (INPUT == null && INPUT_URL == null)
            return new String[]{"INPUT BAM file unspecified"};

        // Todo allow output to a stream, e.g. BuildBamIndex > BaiToText
        if (OUTPUT == null) {
            OUTPUT = new File(IoUtil.basename(INPUT) + BAMIndex.BAMIndexSuffix);
        }
        if (OUTPUT.exists()) {
            OUTPUT.delete();
        }
        return null;
    }
}