/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.programgroups.SamOrBam;

/**
 * @author Yossi Farjoun
 */
@Deprecated
@CommandLineProgramProperties(
        summary = SetNmAndUqTags.USAGE_SUMMARY + SetNmMdAndUqTags.USAGE_DETAILS,
        oneLineSummary = SetNmAndUqTags.USAGE_SUMMARY,
        programGroup = SamOrBam.class
)
public class SetNmAndUqTags extends SetNmMdAndUqTags {
    static final String USAGE_SUMMARY = "DEPRECATED: Use SetNmMdAndUqTags instead.";
    static final String USAGE_DETAILS = "DEPRECATED: Use SetNmMdAndUqTags instead.  This tool takes in a SAM or BAM file (sorted by coordinate) and calculates the NM, MD, and UQ tags by comparing with the reference."+
            "<br />" +
            "This may be needed when MergeBamAlignment was run with SORT_ORDER different from 'coordinate' and thus could not fix\n"+
            "these tags then.<br />"+
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar SetNmAndUqTags \\<br />" +
            "      I=sorted.bam \\<br />" +
            "      O=fixed.bam \\<br />"+
            "</pre>" +
            "<hr />";
}

