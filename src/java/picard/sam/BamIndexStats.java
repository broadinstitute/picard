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

package picard.sam;

import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;

/**
 * Command line program to print statistics from BAM index (.bai) file
 * Statistics include count of aligned and unaligned reads for each reference sequence
 * and a count of all records with no start coordinate.
 * Similar to the 'samtools idxstats' command.
 *
 * @author Martha Borkan
 */
@CommandLineProgramProperties(
        usage = "Generates BAM index statistics, including the number of aligned and unaligned SAMRecords for each reference sequence, " +
                "and the number of SAMRecords with no coordinate." +
                "Input BAM file must have a corresponding index file.\n",
        usageShort = "Generates index statistics from a BAM file",
        programGroup = SamOrBam.class
)
public class BamIndexStats extends CommandLineProgram {

    private static final Log log = Log.getInstance(BamIndexStats.class);

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc="A BAM file to process.")
    public File INPUT;

    /** Stock main method for a command line program. */
    public static void main(final String[] argv) {
        System.exit(new BamIndexStats().instanceMain(argv));
    }

    /**
     * Main method for the program.  Checks that input file is present and
     * readable, then iterates through the index printing meta data to stdout.
     */
    protected int doWork() {

        if (INPUT.getName().endsWith(BAMIndex.BAMIndexSuffix))
               log.warn("INPUT should be BAM file not index file");
        IOUtil.assertFileIsReadable(INPUT);
        BAMIndexMetaData.printIndexStats(INPUT);

        return 0;
    }
}
