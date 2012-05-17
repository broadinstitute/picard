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
import net.sf.picard.util.AdapterPair;
import net.sf.picard.util.ClippingUtility;
import net.sf.picard.util.IlluminaUtil;
import net.sf.samtools.*;

import java.util.*;

import net.sf.picard.illumina.parser.ClusterData;
import net.sf.picard.illumina.parser.ReadData;

/**
 * Takes ClusterData provided by an IlluminaDataProvider into one or two SAMRecords,
 * as appropriate, and optionally marking adapter sequence.  There is one converter per
 * IlluminaBasecallsToSam run, and all the TileProcessors use the same converter.
 * 
 * @author jburke@broadinstitute.org
 */
public class IlluminaBasecallsToSamConverter {


    private final String runBarcode;
    private final String readGroupId;
    private final ReadStructure readStructure;
    private final SamRecordFilter filters = new SolexaNoiseFilter();
    private final boolean isPairedEnd;
    private final boolean isBarcoded;
    private final int [] templateIndices;
    private final int [] barcodeIndices;
    private final AdapterPair[] adaptersToCheck;

    /**
     * Constructor
     *
     * @param runBarcode        Used to construct read names.
     * @param readGroupId       If non-null, set RG attribute on SAMRecord to this.
     * @param readStructure     The expected structure (number of reads and indexes,
     *                          and their length) in the read.
     * @param adapters          The list of adapters to check for in the read
     */
    public IlluminaBasecallsToSamConverter(final String runBarcode,
                                           final String readGroupId,
                                           final ReadStructure readStructure,
                                           final List<IlluminaUtil.IlluminaAdapterPair> adapters) {
        this.runBarcode  = runBarcode;
        this.readGroupId = readGroupId;
        this.readStructure = readStructure;

        this.isPairedEnd = readStructure.templates.length() == 2;
        this.isBarcoded  = !readStructure.barcodes.isEmpty();

        this.adaptersToCheck = new AdapterPair[adapters.size()];
        for (int i = 0; i < adapters.size(); i++) adaptersToCheck[i] = adapters.get(i);

        this.templateIndices = readStructure.templates.getIndices();
        this.barcodeIndices = readStructure.barcodes.getIndices();
    }

    /**
     * Gets the number of non-index reads per cluster
     */
    public int getNumRecordsPerCluster() {
        return readStructure.templates.length();
    }

    /**
     * Creates a new SAM record from the basecall data
     */
    private SAMRecord createSamRecord(final ReadData readData, final SAMFileHeader header, final String readName, final boolean isPf, final boolean firstOfPair, final String unmatchedBarcode) {
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

        // If it's a barcoded run and the read isn't assigned to a barcode, then add the barcode
        // that was read as an optional tag
        if (unmatchedBarcode != null) {
            sam.setAttribute("BC", unmatchedBarcode);
        }

        return sam;
    }

    /**
     * Creates the SAMRecord for each read in the cluster
     */
    public void createSamRecords(final ClusterData cluster, final SAMFileHeader header, final SAMRecord [] recordsOut) {

        final String readName = IlluminaUtil.makeReadName(runBarcode, cluster.getLane(), cluster.getTile(), cluster.getX(), cluster.getY());

        // Get and transform the unmatched barcode, if any, to store with the reads
        String unmatchedBarcode = null;
        if (isBarcoded && cluster.getMatchedBarcode() == null) {
            final byte barcode[][] = new byte[barcodeIndices.length][];
            for (int i = 0; i < barcodeIndices.length; i++) {
                barcode[i] = cluster.getRead(barcodeIndices[i]).getBases();
            }
            unmatchedBarcode = IlluminaUtil.barcodeSeqsToString(barcode).replace('.', 'N'); //TODO: This has a separator, where as in other places we do not use a separator
        }

        final SAMRecord firstOfPair  = createSamRecord(
            cluster.getRead(templateIndices[0]), header, readName, cluster.isPf(), true,unmatchedBarcode);
        recordsOut[0] = firstOfPair;

        SAMRecord secondOfPair = null;

        if(isPairedEnd) {
            secondOfPair  = createSamRecord(
                cluster.getRead(templateIndices[1]), header, readName, cluster.isPf(), false, unmatchedBarcode);
            recordsOut[1] = secondOfPair;
        }

        if (adaptersToCheck.length > 0) {
            // Clip the read
            if (isPairedEnd) {
                ClippingUtility.adapterTrimIlluminaPairedReads(firstOfPair, secondOfPair, adaptersToCheck);
            }
            else {
                ClippingUtility.adapterTrimIlluminaSingleRead(firstOfPair, adaptersToCheck);
            }
        }
    }


}
