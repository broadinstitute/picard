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

import net.sf.picard.PicardException;
import net.sf.picard.filter.AggregateFilter;
import net.sf.picard.filter.SamRecordFilter;
import net.sf.picard.filter.SolexaNoiseFilter;
import net.sf.picard.illumina.parser.ReadStructure;
import net.sf.picard.sam.ReservedTagConstants;
import net.sf.picard.util.ClippingUtility;
import net.sf.picard.util.IlluminaUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.*;

import java.util.Arrays;

import net.sf.picard.illumina.parser.ClusterData;
import net.sf.picard.illumina.parser.ReadData;
import net.sf.samtools.util.StringUtil;

/**
 * Convert ClusterData into SAMRecord.
 * 
 * @author jburke@broadinstitute.org
 */
public class IlluminaBasecallsToSamConverter {

    private final String runBarcode;
    private final String readGroupId;
    private final ReadStructure readStructure;
    private AggregateFilter filters;
    private int barcodeIndex;
    private final boolean isPairedEnd;
    private final boolean isBarcoded;
    private static final Log log = Log.getInstance(IlluminaBasecallsToSamConverter.class);

    /**
     * Constructor
     *
     * @param runBarcode        Used to construct read names.
     * @param readGroupId       If non-null, set RG attribute on SAMRecord to this.
     */
    public IlluminaBasecallsToSamConverter(final String runBarcode,
                                           final String readGroupId,
                                           final ReadStructure readStructure) {
        this.runBarcode  = runBarcode;
        this.readGroupId = readGroupId;
        this.readStructure = readStructure;

        if(readStructure.numTemplates > 3) {
            throw new PicardException("IlluminaBasecallsToSamConverter does not support more than 2 template reads.  Number of template reads found in configuration: " + readStructure.numTemplates);
        }

        if(readStructure.numBarcodes > 1) {
            throw new PicardException("IlluminaBasecallsToSamConverter does not support more than 1 barcode read.  Number of template reads found in configuration: " + readStructure.numBarcodes);
        }

        this.isPairedEnd = readStructure.numTemplates == 2;
        this.isBarcoded  = readStructure.numBarcodes > 0;

        this.barcodeIndex = -1;
        if(readStructure.barcodeIndices.length > 0) { //Should only be one for now
            this.barcodeIndex = readStructure.barcodeIndices[0];
        }

        initializeFilters();
    }

    private void initializeFilters() {
        filters = new AggregateFilter(Arrays.asList(
            (SamRecordFilter)new SolexaNoiseFilter()
        ));
    }
    private String createReadName(final ClusterData cluster) {
        return IlluminaUtil.makeReadName(runBarcode, cluster.getLane(), cluster.getTile(), cluster.getX(), cluster.getY());
    }

    public int getNumRecordsPerCluster() {
        return readStructure.numTemplates;
    }

    private SAMRecord createSamRecord(final ReadData readData, final SAMFileHeader header, final String readName, final boolean isPf, final boolean firstOfPair, final ReadData unmatchedBarcodeRead) {
        final SAMRecord sam = new SAMRecord(header);
        sam.setReadName(readName);
        sam.setReadBases(readData.getBases());
        sam.setBaseQualities(readData.getQualities());

        // Flag values
        sam.setReadPairedFlag(isPairedEnd);
        sam.setReadUnmappedFlag(true);
        sam.setReadFailsVendorQualityCheckFlag(!isPf);
        if (isPairedEnd) {
            sam.setMateUnmappedFlag(true);
            sam.setFirstOfPairFlag(firstOfPair);
            sam.setSecondOfPairFlag(!firstOfPair);
        }

        if (filters.filterOut(sam)) {
            sam.setAttribute(ReservedTagConstants.XN, 1);
        }

        if (this.readGroupId != null) {
            sam.setAttribute("RG", readGroupId);
        }

        // If it's a barcoded run and the read isn't assigned to a barcode, then add the barcode read as an optional tag
        if (unmatchedBarcodeRead != null) {
            sam.setAttribute("BC", StringUtil.bytesToString(unmatchedBarcodeRead.getBases()).replace('.', 'N'));
        }

        return sam;
    }

    public void createSamRecords(final ClusterData cluster, final SAMFileHeader header, boolean markAdapter, final SAMRecord [] recordsOut) {
        final String readName = createReadName(cluster);

        SAMRecord firstOfPair  = createSamRecord(cluster.getRead(readStructure.templateIndices[0]), header, readName, cluster.isPf(), true,  (cluster.getMatchedBarcode() == null && isBarcoded) ? cluster.getRead(barcodeIndex) : null);
        recordsOut[0] = firstOfPair;

        SAMRecord secondOfPair = null;

        if(isPairedEnd) {
            secondOfPair  = createSamRecord(cluster.getRead(readStructure.templateIndices[1]), header, readName, cluster.isPf(), false, (cluster.getMatchedBarcode() == null && isBarcoded) ? cluster.getRead(barcodeIndex) : null);
            recordsOut[1] = secondOfPair;
        }

        if (markAdapter) {
            if (isPairedEnd){
                assert (firstOfPair.getFirstOfPairFlag() && secondOfPair.getSecondOfPairFlag());
                String warnString = ClippingUtility.adapterTrimIlluminaPairedReads(firstOfPair, secondOfPair,
                        isBarcoded ? IlluminaUtil.IlluminaAdapterPair.INDEXED.adapterPair
                                   : IlluminaUtil.IlluminaAdapterPair.PAIRED_END.adapterPair);
               if (warnString != null){
                    log.debug("Adapter trimming " + warnString);
                }
            } else {
                ClippingUtility.adapterTrimIlluminaSingleRead(firstOfPair,
                    isBarcoded   ? IlluminaUtil.IlluminaAdapterPair.INDEXED.adapterPair
                                 : IlluminaUtil.IlluminaAdapterPair.PAIRED_END.adapterPair);
                // note if not barcoded, it could instead be SINGLE_END
                // we're assuming one read of paired_end is more common
                // Note, the single_end adapters have first 13-18 bases in common with paired_end
                // We could, alternatively, try both and use the one that matched more adapter
            }
        }
    }


}
