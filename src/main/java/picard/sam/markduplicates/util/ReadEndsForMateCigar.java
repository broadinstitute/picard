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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.SamRecordWithOrdinal;
import picard.PicardException;

import java.util.List;
import java.util.Set;

/**
 * A class to store individual records for MarkDuplicatesWithMateCigar.  This aids in comparing records to determine which need to
 * be compared when we mark duplicates.  We also store the original SAMRecord and its ordinal in the input file (in SamRecordWithOrdinal) to
 * access optional tags (mate cigar) and other information.
 */
public class ReadEndsForMateCigar extends ReadEnds {
    // to see if either end is unmapped
    byte hasUnmapped = 0;

    // we need this reference so we can access the mate cigar among other things
    public SamRecordWithOrdinal samRecordWithOrdinal = null;

    /**
     * Physical locations used for optical duplicate tracking.  This is only stored for paired end reads where both ends are mapped,
     * and when we see the first mate.
     */
    private PhysicalLocationForMateCigarSet locationSet = null;

    /** Builds a read ends object that represents a single read. */
    public ReadEndsForMateCigar(final SAMFileHeader header, final SamRecordWithOrdinal samRecordWithOrdinal,
                                final OpticalDuplicateFinder opticalDuplicateFinder, final short libraryId) {

        this.readGroup = -1;
        this.tile = -1;
        this.x = this.y = -1;
        this.read2ReferenceIndex = this.read2Coordinate = -1;
        this.hasUnmapped = 0;

        this.samRecordWithOrdinal = samRecordWithOrdinal;

        final SAMRecord record = this.samRecordWithOrdinal.getRecord();

        this.read1ReferenceIndex = record.getReferenceIndex();
        this.read1Coordinate = record.getReadNegativeStrandFlag() ? record.getUnclippedEnd() : record.getUnclippedStart();
        if (record.getReadUnmappedFlag()) {
            throw new PicardException("Found an unexpected unmapped read");
        }

        if (record.getReadPairedFlag() && !record.getReadUnmappedFlag() && !record.getMateUnmappedFlag()) {
            this.read2ReferenceIndex = record.getMateReferenceIndex();
            this.read2Coordinate = record.getMateNegativeStrandFlag() ? SAMUtils.getMateUnclippedEnd(record) : SAMUtils.getMateUnclippedStart(record);

            // set orientation
            this.orientation = ReadEnds.getOrientationByte(record.getReadNegativeStrandFlag(), record.getMateNegativeStrandFlag());

            // Set orientationForOpticalDuplicates, which always goes by the first then the second end for the strands.  NB: must do this
            // before updating the orientation later.
            if (record.getReadPairedFlag()) {
                if (record.getFirstOfPairFlag()) {
                    this.orientationForOpticalDuplicates = ReadEnds.getOrientationByte(record.getReadNegativeStrandFlag(), record.getMateNegativeStrandFlag());
                } else {
                    this.orientationForOpticalDuplicates = ReadEnds.getOrientationByte(record.getMateNegativeStrandFlag(), record.getReadNegativeStrandFlag());
                }
            }
        } else {
            this.orientation = record.getReadNegativeStrandFlag() ? ReadEndsForMateCigar.R : ReadEndsForMateCigar.F;
        }

        // Fill in the library ID
        this.libraryId = libraryId;

        // Is this unmapped or its mate?
        if (record.getReadUnmappedFlag() || (record.getReadPairedFlag() && record.getMateUnmappedFlag())) {
            this.hasUnmapped = 1;
        }

        // Fill in the location information for optical duplicates
        if (opticalDuplicateFinder.addLocationInformation(record.getReadName(), this)) {
            // calculate the RG number (nth in list)
            // NB: could this be faster if we used a hash?
            this.readGroup = 0;
            final String rg = (String) record.getAttribute("RG");
            final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
            if (rg != null && readGroups != null) {
                for (final SAMReadGroupRecord readGroup : readGroups) {
                    if (readGroup.getReadGroupId().equals(rg)) break;
                    else this.readGroup++;
                }
            }
        }
    }

    /** Creates a shallow copy from the "other" */
    public ReadEndsForMateCigar(final ReadEndsForMateCigar other, final SamRecordWithOrdinal samRecordWithOrdinal) {
        this.readGroup = other.readGroup;
        this.tile = other.tile;
        this.x = other.x;
        this.y = other.y;
        this.read1ReferenceIndex = other.read1ReferenceIndex;
        this.read1Coordinate = other.read1Coordinate;
        this.read2ReferenceIndex = other.read2ReferenceIndex;
        this.read2Coordinate = other.read2Coordinate;
        this.hasUnmapped = other.hasUnmapped;
        this.samRecordWithOrdinal = samRecordWithOrdinal;
        this.orientation = other.orientation;
        this.libraryId = other.libraryId;
    }

    /** A number of convenience functions */
    public SamRecordWithOrdinal getSamRecordIndex() { return this.samRecordWithOrdinal; }

    public SAMRecord getRecord() { return this.samRecordWithOrdinal.getRecord(); }

    public String getRecordReadName() { return this.samRecordWithOrdinal.getRecord().getReadName(); }

    @Override
    public boolean isPaired() { return this.getRecord().getReadPairedFlag(); }

    /** Gets the read ends for optical duplicate tracking */
    public Set<ReadEnds> getReadEndSetForOpticalDuplicates() {
        if (null == this.locationSet) throw new PicardException("Already called getReadEndSetForOpticalDuplicates");
        final Set<ReadEnds> locationSet = this.locationSet.getReadEnds();
        this.locationSet = null;
        return locationSet;
    }

    public PhysicalLocationForMateCigarSet getLocationSet() {
        return this.locationSet;
    }

    public PhysicalLocationForMateCigarSet removeLocationSet() {
        final PhysicalLocationForMateCigarSet locationSet = this.locationSet;
        this.locationSet = null;
        return locationSet;
    }

    public void setLocationSet(final PhysicalLocationForMateCigarSet locationSet) {
        this.locationSet = locationSet;
    }

}
