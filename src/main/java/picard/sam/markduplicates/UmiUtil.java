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
import picard.PicardException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import org.apache.commons.lang3.StringUtils;
import picard.sam.markduplicates.util.MarkDuplicatesUtil;

import java.util.regex.Pattern;

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
    static final String INFERRED_UMI_TRANSIENT_TAG = "inferredUmi";
    static final Pattern ALLOWED_UMI = Pattern.compile("^[ATCGNatcgn-]*$");
    /**
     * Creates a top-strand normalized duplex UMI.
     * Single stranded UMIs are by definition already top-strand normalized, they require no transformation.
     * Duplex UMIs that come from a top strand read are also by definition, top-strand normalized. A duplex
     * UMI from a bottom strand can be normalized to be identical to the read from its corresponding top strand
     * by swapping the content of the UMI around the "-" found in duplex UMIs.  For example, a bottom strand
     * duplex UMI reading ATC-CGG when top-strand normalized will read CGG-ATC.
     *
     * @param record SAM record to retrieve UMI from.
     * @param umiTag The tag used in the bam file that designates the UMI, null returns null
     * @return Normalized Duplex UMI.  If the UMI isn't duplex, it returns the UMI unaltered.
     */
    static String getTopStrandNormalizedUmi(final SAMRecord record, final String umiTag, final boolean duplexUmi) {
        if (umiTag == null) {
            return null;
        }

        final String umi = record.getStringAttribute(umiTag);

        if (umi == null) {
            return null;
        }

        if (!ALLOWED_UMI.matcher(umi).matches()) {
            throw new PicardException("UMI found with illegal characters.  UMIs must match the regular expression ^[ATCGNatcgn-]*$.");
        }

        if (!duplexUmi) {
            return umi;
        }

        final String[] split = umi.split(DUPLEX_UMI_DELIMITER);
        if (split.length != 2) {
            throw new PicardException("Duplex UMIs must be of the form X-Y where X and Y are equal length UMIs, for example AT-GA.  Found UMI, " + umi);
        }

        switch (getStrand(record)) {
            case BOTTOM:
                return split[1] + DUPLEX_UMI_DELIMITER + split[0];
                default:
                    return umi;
        }
    }

    /**
     * Determines if the read represented by a SAM record belongs to the top or bottom strand
     * or if it cannot determine strand position due to one of the reads being unmapped.
     * Top strand is defined as having a read 1 unclipped 5' coordinate
     * less than the read 2 unclipped 5' coordinate.  If a read is unmapped
     * we do not attempt to determine the strand to which the read or its mate belongs.
     * If the mate belongs to a different contig from the read, then the reference
     * index for the read and its mate is used in leu of the unclipped 5' coordinate.
     * @param rec Record to determine top or bottom strand
     * @return Top or bottom strand, unknown if it cannot be determined.
     */
    static ReadStrand getStrand(final SAMRecord rec) {
        if (!MarkDuplicatesUtil.pairedForMarkDuplicates(rec)) {
            return ReadStrand.UNKNOWN;
        }

        // If the read pair are aligned to different contigs we use
        // the reference index to determine relative 5' coordinate ordering.
        // Both the read and its mate should not have their unmapped flag set to true.
        if (!rec.getReferenceIndex().equals(rec.getMateReferenceIndex())) {
            if (rec.getFirstOfPairFlag() == (rec.getReferenceIndex() < rec.getMateReferenceIndex())) {
                return ReadStrand.TOP;
            } else {
                return ReadStrand.BOTTOM;
            }
        }

        final int read5PrimeStart = (rec.getReadNegativeStrandFlag()) ? rec.getUnclippedEnd() : rec.getUnclippedStart();
        final int mate5PrimeStart = (rec.getMateNegativeStrandFlag()) ? SAMUtils.getMateUnclippedEnd(rec) : SAMUtils.getMateUnclippedStart(rec);

        if (rec.getFirstOfPairFlag() == (read5PrimeStart <= mate5PrimeStart)) {
            return ReadStrand.TOP;
        } else {
            return ReadStrand.BOTTOM;
        }
    }

    /**
     * Set molecular identifier tag of record using UMI and alignment start position along with top and bottom strand information
     * @param rec SAMRecord to set molecular index of
     * @param assignedUmi Assigned or inferred UMI to use in the molecular index tag
     * @param molecularIdentifierTag SAM tag to use as molecular identifier, if null no modification to the record will be made
     * @param duplexUmis Treat UMI as duplex, if true /A and /B will be added to denote top and bottom strands respectively
     */
    static void setMolecularIdentifier(final SAMRecord rec, final String assignedUmi, final String molecularIdentifierTag, final boolean duplexUmis) {

        if (molecularIdentifierTag == null) {
            return;
        }

        final StringBuilder molecularIdentifier = new StringBuilder();
        molecularIdentifier.append(rec.getContig());
        molecularIdentifier.append(CONTIG_SEPARATOR);
        molecularIdentifier.append(rec.getReadNegativeStrandFlag() ? rec.getAlignmentStart() : rec.getMateAlignmentStart());
        molecularIdentifier.append(UMI_NAME_SEPARATOR);
        molecularIdentifier.append(assignedUmi);

        if (duplexUmis) {
            // Reads whose strand position can be determined will have their
            // molecularIdentifier set to include an identifier appended that
            // indicates top or bottom strand.
            ReadStrand strand = getStrand(rec);
            switch (strand) {
                case TOP:
                    molecularIdentifier.append(TOP_STRAND_DUPLEX);
                    break;
                case BOTTOM:
                    molecularIdentifier.append(BOTTOM_STRAND_DUPLEX);
                    break;
                default:
                    // If we can't determine strand position nothing
                    // is appended to the molecularIdentifier.
                    break;
            }
        }
        rec.setAttribute(molecularIdentifierTag, molecularIdentifier.toString());
    }

    /**
     * Get the length of UMI ignoring the effect of hyphens
     * @param umi String that represents a unique molecular index
     * @return Length of the UMI ignoring hyphen characters
     */
    static int getUmiLength(final String umi) {
        return umi.length() - StringUtils.countMatches(umi, UmiUtil.DUPLEX_UMI_DELIMITER);
    }

    /**
     * An enum to hold the strand position (TOP or BOTTOM) of a read.
     */
    public enum ReadStrand {
        TOP,
        BOTTOM,
        UNKNOWN
    }
}
