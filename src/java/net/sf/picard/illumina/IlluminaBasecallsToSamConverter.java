/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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
package net.sf.picard.illumina;

import net.sf.picard.illumina.parser.IlluminaFileUtil;
import net.sf.picard.filter.AggregateFilter;
import net.sf.picard.filter.SamRecordFilter;
import net.sf.picard.filter.SolexaNoiseFilter;
import net.sf.picard.sam.ReservedTagConstants;
import net.sf.picard.util.IlluminaUtil;
import net.sf.samtools.*;

import java.util.Arrays;

import net.sf.picard.illumina.parser.IlluminaReadData;
import net.sf.picard.illumina.parser.IlluminaEndData;
import net.sf.samtools.util.StringUtil;

/**
 * Convert IlluminaReadData into SAMRecord.
 * 
 * @author alecw@broadinstitute.org
 */
public class IlluminaBasecallsToSamConverter {

    private final String runBarcode;
    private final String readGroupId;
    private AggregateFilter filters;

    /**
     * Constructor
     *
     * @param runBarcode        Used to construct read names.
     * @param readGroupId       If non-null, set RG attribute on SAMRecord to this.
     */
    public IlluminaBasecallsToSamConverter(final String runBarcode,
                                           final String readGroupId) {
        this.runBarcode = runBarcode;
        this.readGroupId = readGroupId;
        initializeFilters();
    }

    private void initializeFilters() {
        filters = new AggregateFilter(Arrays.asList(
            (SamRecordFilter)new SolexaNoiseFilter()
        ));
    }
    private String createReadName(final IlluminaReadData ird) {
        return IlluminaUtil.makeReadName(runBarcode, ird.getLane(), ird.getTile(), ird.getX(), ird.getY());
    }

    /**
     * Creates a SAMRecord from Illumina Basecall data
     *
     * @param ird           The IlluminaReadData to use in populating the SAMRecord
     * @param isFirstRead   whether this is the first read of a pair
     * @param readName      The read name to use, or null if it should be constructed here.  For paired-end runs
     *                      the same read name can be used for both ends.
     * @return SAMRecord    fully populated SAMRecord
     */
    public SAMRecord createSamRecord(final IlluminaReadData ird, final boolean isFirstRead, final SAMFileHeader header, 
                                     final String readName) {
        final SAMRecord sam = new SAMRecord(header);
        sam.setReadName(readName != null? readName: createReadName(ird));
        final IlluminaEndData end = (isFirstRead? ird.getFirstEnd(): ird.getSecondEnd());
        sam.setReadBases(end.getBases());
        sam.setBaseQualities(end.getQualities());

        // Flag values
        sam.setReadPairedFlag(ird.isPairedEnd());
        sam.setReadUnmappedFlag(true);
        sam.setReadFailsVendorQualityCheckFlag(!ird.isPf());
        if (ird.isPairedEnd()) {
            sam.setMateUnmappedFlag(true);
            sam.setFirstOfPairFlag(isFirstRead);
            sam.setSecondOfPairFlag(!isFirstRead);
        }

        if (filters.filterOut(sam)) {
            sam.setAttribute(ReservedTagConstants.XN, 1);
        }

        if (this.readGroupId != null) {
            sam.setAttribute("RG", readGroupId);
        }

        // If it's a barcoded run and the read isn't assigned to a barcode, then add the barcode read as an optional tag
        if (ird.getMatchedBarcode() == null && ird.getBarcodeRead() != null) {
            sam.setAttribute("BC", StringUtil.bytesToString(ird.getBarcodeRead().getBases()).replace('.', 'N'));
        }

        return sam;
    }
}
