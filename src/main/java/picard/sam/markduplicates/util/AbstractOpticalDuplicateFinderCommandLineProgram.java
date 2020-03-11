/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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

package picard.sam.markduplicates.util;

import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.CommandLineProgram;

/**
 * Abstract class that holds parameters and methods common to classes that optical duplicate detection.  We put them here so that
 * the explanation about how read names are parsed is in once place
 *
 * @author Tim Fennell
 */
public abstract class AbstractOpticalDuplicateFinderCommandLineProgram extends CommandLineProgram {
    protected static Log LOG = Log.getInstance(AbstractOpticalDuplicateFinderCommandLineProgram.class);

    @Argument(doc = "MarkDuplicates can use the tile and cluster positions to estimate the rate of optical duplication " +
            "in addition to the dominant source of duplication, PCR, to provide a more accurate estimation of library size. " +
            "By default (with no READ_NAME_REGEX specified), MarkDuplicates will attempt to extract coordinates " +
            "using a split on ':' (see Note below).  " +
            "Set READ_NAME_REGEX to 'null' to disable optical duplicate detection. " +
            "Note that without optical duplicate counts, library size estimation will be less accurate. " +
            "If the read name does not follow a standard Illumina colon-separation convention, but does contain tile and x,y coordinates, " + 
            "a regular expression can be specified to extract three variables: tile/region, x coordinate and y coordinate from a read name. " +
            "The regular expression must contain three capture groups for the three variables, in order. " +
            "It must match the entire read name. " +
            "  e.g. if field names were separated by semi-colon (';') this example regex could be specified " +
            "     (?:.*;)?([0-9]+)[^;]*;([0-9]+)[^;]*;([0-9]+)[^;]*$ " +
            "Note that if no READ_NAME_REGEX is specified, the read name is split on ':'. " +
            "  For 5 element names, the 3rd, 4th and 5th elements are assumed to be tile, x and y values. " +
            "  For 7 element names (CASAVA 1.8), the 5th, 6th, and 7th elements are assumed to be tile, x and y values.",
            optional = true)
    public String READ_NAME_REGEX = OpticalDuplicateFinder.DEFAULT_READ_NAME_REGEX;

    @Argument(doc = "The maximum offset between two duplicate clusters in order to consider them optical duplicates. The default " +
            "is appropriate for unpatterned versions of the Illumina platform. For the patterned flowcell models, 2500 is more" +
            "appropriate. For other platforms and models, users should experiment to find what works best.")
    public int OPTICAL_DUPLICATE_PIXEL_DISTANCE = OpticalDuplicateFinder.DEFAULT_OPTICAL_DUPLICATE_DISTANCE;

    @Argument(doc = "This number is the maximum size of a set of duplicate reads for which we will attempt to determine " +
            "which are optical duplicates.  Please be aware that if you raise this value too high and do encounter a very " +
            "large set of duplicate reads, it will severely affect the runtime of this tool.  To completely disable this check, " +
            "set the value to -1.")
    public long MAX_OPTICAL_DUPLICATE_SET_SIZE = OpticalDuplicateFinder.DEFAULT_MAX_DUPLICATE_SET_SIZE;

    // The tool with which to find optical duplicates
    protected OpticalDuplicateFinder opticalDuplicateFinder = null;

    // Needed for testing
    public void setupOpticalDuplicateFinder() {
        this.opticalDuplicateFinder = new OpticalDuplicateFinder(READ_NAME_REGEX, OPTICAL_DUPLICATE_PIXEL_DISTANCE, MAX_OPTICAL_DUPLICATE_SET_SIZE, LOG);
    }

    @Override
    protected String[] customCommandLineValidation() {
        setupOpticalDuplicateFinder();
        return super.customCommandLineValidation();
    }
}
