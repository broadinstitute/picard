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
import org.apache.commons.lang3.StringUtils;
import picard.PicardException;

/**
 *
 * A collection of functions for use in processing UMIs
 *
 * @author mduran
 */

class UmiUtil {

    public static final String DUPLEX_UMI_DELIMITER = "-";
    public static final String TOP_STRAND_DUPLEX = "/A";
    public static final String BOTTOM_STRAND_DUPLEX = "/B";
    public static final String CONTIG_SEPARATOR = ":";
    public static final String UMI_NAME_SEPARATOR = "/";
    static final String INFERRED_UMI_TAG = "inferredUmi";

    /**
     * Creates a top-strand normalized duplex UMI.
     * Single stranded UMIs are by definition already top-strand normalized, they require no transformation.
     * Duplex UMIs that come from a top strand read are also by definition, top-strand normalized.
     * A duplex UMI from a bottom strand can be normalized to be identical to the UMI
     * read from its corresponding top strand by swapping the content of the
     * UMI around the "-" found in duplex UMIs.  For example, a bottom strand
     * duplex UMI reading ATC-CGG when top-strand normalized will read CGG-ATC.
     *
     * @param record SAM record to retrieve UMI from.
     * @param umiTag The tag used in the bam file that designates the UMI.
     * @return Normalized Duplex UMI.  If the UMI isn't duplex, it returns the UMI unaltered.
     */
    static String getTopStrandNormalizedUmi(final SAMRecord record, final String umiTag, final boolean duplexUmi) {
        if (umiTag == null) return null;

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
        } else {
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

    /**
     * Creates a string that uniquely identifies a molecule using its UMI and alignment start position
     * @param rec SAMRecord to create molecular identifier of
     * @param umiTag Tag that contains the UMI record
     * @return String that uniquely identifies a fragment
     */
    static String molecularIdentifierString(final SAMRecord rec, final String umiTag) {
        return rec.getContig() + ":" + rec.getAlignmentStart() + getTopStrandNormalizedUmi(rec, umiTag, true);
    }

    /**
     * Set molecular index tag of record using UMI and alignment start position along with top and bottom strand information
     * @param rec
     * @param assignedUmi
     * @param molecularIdentifierTag
     * @param duplexUmis
     */
    static void setMolecularIndex(final SAMRecord rec, final String assignedUmi, final String molecularIdentifierTag, final boolean duplexUmis) {

        final String fragmentStartPosition;
        if (rec.getReadNegativeStrandFlag()) {
            fragmentStartPosition = rec.getContig() + CONTIG_SEPARATOR + rec.getAlignmentStart();
        } else {
            fragmentStartPosition = rec.getContig() + CONTIG_SEPARATOR + rec.getMateAlignmentStart();
        }

        final String strandPosition;
        if (duplexUmis) {
            if (isTopStrand(rec)) {
                strandPosition = TOP_STRAND_DUPLEX;
            } else {
                strandPosition = BOTTOM_STRAND_DUPLEX;
            }
            rec.setAttribute(molecularIdentifierTag, fragmentStartPosition + UMI_NAME_SEPARATOR + assignedUmi + strandPosition);
        } else {
            rec.setAttribute(molecularIdentifierTag, fragmentStartPosition + UMI_NAME_SEPARATOR + assignedUmi);
        }
    }

    /**
     * Get the length of UMI ignoring the effect of hyphens
     * @param umi String that represents a unique molecular index
     * @return Length of the UMI ignoring hyphen characters
     */
    static int getUmiLength(final String umi) {
        return umi.length() - StringUtils.countMatches(umi, UmiUtil.DUPLEX_UMI_DELIMITER);
    }
}
