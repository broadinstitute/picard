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

import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.*;
import htsjdk.samtools.util.Log;

/**
 * Abstract class that holds parameters and methods common to classes that optical duplicate detection.  We put them here so that
 * the explanation about how read names are parsed is in once place
 *
 * @author Tim Fennell
 */
public abstract class AbstractOpticalDuplicateFinderCommandLineProgram extends CommandLineProgram {
    protected static Log LOG = Log.getInstance(AbstractOpticalDuplicateFinderCommandLineProgram.class);


    @Argument(doc = "Regular expression that can be used to parse read names in the incoming SAM file. Read names are " +
            "parsed to extract three variables: tile/region, x coordinate and y coordinate. These values are used " +
            "to estimate the rate of optical duplication in order to give a more accurate estimated library size. " +
            "Set this option to null to disable optical duplicate detection, e.g. for RNA-seq " +
            "or other data where duplicate sets are extremely large and estimating library complexity is not an aim. " +
            "Note that without optical duplicate counts, library size estimation will be inaccurate. " +
            "The regular expression should contain three capture groups for the three variables, in order. " +
            "It must match the entire read name. " +
            "Note that if the default regex is specified, a regex match is not actually done, but instead the read name " +
            " is split on colon character. " +
            "For 5 element names, the 3rd, 4th and 5th elements are assumed to be tile, x and y values. " +
            "For 7 element names (CASAVA 1.8), the 5th, 6th, and 7th elements are assumed to be tile, x and y values.",
            optional = true)
    public String READ_NAME_REGEX = OpticalDuplicateFinder.DEFAULT_READ_NAME_REGEX;

    @Argument(doc = "The maximum offset between two duplicate clusters in order to consider them optical duplicates. The default " +
            "is appropriate for unpatterned versions of the Illumina platform. For the patterned flowcell models, 2500 is more" +
            "appropriate. For other platforms and models, users should experiment to find what works best.")
    public int OPTICAL_DUPLICATE_PIXEL_DISTANCE = OpticalDuplicateFinder.DEFAULT_OPTICAL_DUPLICATE_DISTANCE;

    // The tool with which to find optical duplicates
    protected OpticalDuplicateFinder opticalDuplicateFinder = null;

    // Needed for testing
    public void setupOpticalDuplicateFinder() {
        this.opticalDuplicateFinder = new OpticalDuplicateFinder(READ_NAME_REGEX, OPTICAL_DUPLICATE_PIXEL_DISTANCE, LOG);
    }

    @Override
    protected String[] customCommandLineValidation() {
        setupOpticalDuplicateFinder();
        return super.customCommandLineValidation();
    }
}
