/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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

package picard.sam.markduplicates;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import picard.PicardException;

/**
 *
 * A collection of functions for use in processing UMIs
 *
 * @author mduran
 */

class UmiUtil {

    static final String DUPLEX_UMI_DELIMITER = "-";

    /**
     * Creates a top-strand normalized duplex UMI.
     * Duplex UMIs that come from a top strand read are by definition, top-strand normalized.
     * The UMI from a bottom strand can be normalized to be identical to the UMI
     * read from its corresponding top strand by swapping the content of the
     * UMI around the "-" found in duplex UMIs.  For example, a bottom strand
     * UMI reading ATC-CGG when top-strand normalized will read CGG-ATC.
     *
     * @param record SAM record to retrieve UMI from.
     * @param umiTag The tag used in the bam file that designates the UMI.
     * @return Normalized Duplex UMI.  If the UMI isn't duplex, it returns the UMI unaltered.
     */
    static String getTopStrandNormalizedDuplexUMI(final SAMRecord record, final String umiTag, final boolean duplexUmi) {
        final String umi = record.getStringAttribute(umiTag);

        if (umi == null) return null;

        if (duplexUmi) {
            final String[] split = umi.split(DUPLEX_UMI_DELIMITER);
            if (split.length != 2) {
                throw new PicardException("Duplex UMIs must be of the form X-Y where X and Y are equal length UMIs, for example AT-GA.  Found UMI, " + umi);
            }

            if (isTopStrand(record)) {
                return split[0] + DUPLEX_UMI_DELIMITER + split[1];
            } else {
                return split[1] + DUPLEX_UMI_DELIMITER + split[0];
            }
        }
        else {
            return umi;
        }
    }

    /**
     * Determines if the read represented by a SAM record belongs to the top or bottom strand.
     * Top strand is defined as having a Read 1 unclipped 5' coordinate
     * less than the Read 2 unclipped 5' coordinate.  If a read is unmapped
     * it is considered to have an unclipped 5' coordinate of 0.
     * @param rec Record to determine top or bottom strand
     * @return Top or bottom strand, true (top), false (bottom).
     */
    static boolean isTopStrand(final SAMRecord rec) {

        final int read5PrimeStart = (rec.getReadNegativeStrandFlag()) ? rec.getUnclippedEnd() : rec.getUnclippedStart();
        final int mate5PrimeStart = (rec.getMateNegativeStrandFlag()) ? SAMUtils.getMateUnclippedEnd(rec) : SAMUtils.getMateUnclippedStart(rec);
        return rec.getFirstOfPairFlag() == (read5PrimeStart < mate5PrimeStart);
    }
}
