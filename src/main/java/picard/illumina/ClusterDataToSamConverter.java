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
package picard.illumina;

import htsjdk.samtools.ReservedTagConstants;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SolexaNoiseFilter;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import picard.PicardException;
import picard.fastq.IlluminaReadNameEncoder;
import picard.fastq.ReadNameEncoder;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.ReadData;
import picard.illumina.parser.ReadStructure;
import picard.util.AdapterMarker;
import picard.util.AdapterPair;
import picard.util.IlluminaUtil;

import java.util.*;

/**
 * Takes ClusterData provided by an IlluminaDataProvider into one or two SAMRecords,
 * as appropriate, and optionally marking adapter sequence.  There is one converter per
 * IlluminaBasecallsToSam run, and all the TileProcessors use the same converter.
 *
 * @author jburke@broadinstitute.org
 */
public class ClusterDataToSamConverter implements
        BasecallsConverter.ClusterDataConverter<IlluminaBasecallsToSam.SAMRecordsForCluster> {

    private final String readGroupId;
    private final SamRecordFilter filters = new SolexaNoiseFilter();
    private final boolean isPairedEnd;
    private final boolean hasSampleBarcode;
    private final boolean hasMolecularBarcode;
    private final int[] templateIndices;
    private final int[] sampleBarcodeIndices;
    private final int[] molecularBarcodeIndices;

    private final AdapterMarker adapterMarker;
    private final int outputRecordsPerCluster;
    private final ReadNameEncoder readNameEncoder;

    // TODO: add RX and QX to the list of SAMTags and change this. initial discussion
    // TODO: here:
    // TODO: - https://github.com/broadinstitute/picard/issues/287
    // TODO: - HTS-spec issue: https://github.com/samtools/hts-specs/issues/109
    // TODO: - https://github.com/samtools/hts-specs/pull/119
    // TODO: also add ~ and - to some htsjdk file and reference here.
    private String MOLECULAR_INDEX_TAG = "RX";
    private String MOLECULAR_INDEX_QUALITY_TAG = "QX";
    private static final String MOLECULAR_INDEX_ = "-";
    private static final String MOLECULAR_INDEX_QUALITY_DELIMITER = "~";
    private static final Character MISSING_BARCODE = '.';
    private static final Character MISSING_BARCODE_BASE = 'N';

    private List<String> tagPerMolecularIndex = Collections.emptyList();

    private PopulateBarcode barcodePopulationStrategy;
    private boolean includeQualitiesWithBarcode;
    // used also in IlluminaBasecallsToSam
    enum PopulateBarcode implements CommandLineParser.ClpEnum {
        ORPHANS_ONLY("Put barcodes only into the records that were not assigned to any declared barcode."),
        INEXACT_MATCH("Put barcodes into records for which an exact match with a declared barcode was not found."),
        ALWAYS("Put barcodes into all the records.");

        private final String description;

        PopulateBarcode(final String description) {
            this.description = description;
        }

        @Override
        public String getHelpDoc() {
            return description;
        }
    }

    /**
     * Constructor
     *
     * @param runBarcode                Used to construct read names.
     * @param readGroupId               If non-null, set RG attribute on SAMRecord to this.
     * @param readStructure             The expected structure (number of reads and indexes,
     *                                  and their length) in the read.
     * @param adapters                  The list of adapters to check for in the read
     * @param barcodePopulationStrategy When to populate BC tag?
     */
    public ClusterDataToSamConverter(final String runBarcode,
                                     final String readGroupId,
                                     final ReadStructure readStructure,
                                     final List<AdapterPair> adapters,
                                     final PopulateBarcode barcodePopulationStrategy,
                                     final boolean includeQualitiesWithBarcode) {
        this.barcodePopulationStrategy = barcodePopulationStrategy;
        this.includeQualitiesWithBarcode = includeQualitiesWithBarcode;
        this.readGroupId = readGroupId;

        this.readNameEncoder = new IlluminaReadNameEncoder(runBarcode);

        this.isPairedEnd = readStructure.templates.length() == 2;
        this.hasSampleBarcode = readStructure.hasSampleBarcode();
        this.hasMolecularBarcode = !readStructure.molecularBarcode.isEmpty();

        if (adapters.isEmpty()) {
            this.adapterMarker = null;
        } else {
            this.adapterMarker = new AdapterMarker(adapters.toArray(new AdapterPair[adapters.size()]));
        }

        this.templateIndices = readStructure.templates.getIndices();
        this.sampleBarcodeIndices = readStructure.sampleBarcodes.getIndices();
        this.molecularBarcodeIndices = readStructure.molecularBarcode.getIndices();

        this.outputRecordsPerCluster = readStructure.templates.length();
    }

    /**
     * Sets the SAM tag to use to store the molecular index bases.  If multiple molecular indexes exist, it will concatenate them
     * and store them in this tag.
     */
    public ClusterDataToSamConverter withMolecularIndexTag(final String molecularIndexTag) {
        if (molecularIndexTag == null) throw new IllegalArgumentException("Molecular index tag was null");
        this.MOLECULAR_INDEX_TAG = molecularIndexTag;
        return this;
    }

    /**
     * Sets the SAM tag to use to store the molecular index base qualities.  If multiple molecular indexes exist, it will concatenate them
     * and store them in this tag.
     */
    public ClusterDataToSamConverter withMolecularIndexQualityTag(final String molecularIndexQualityTag) {
        if (molecularIndexQualityTag == null) {
            throw new IllegalArgumentException("Molecular index quality tag was null");
        }
        this.MOLECULAR_INDEX_QUALITY_TAG = molecularIndexQualityTag;
        return this;
    }

    /**
     * Sets the SAM tags to use to store the bases each molecular index.  This will only be used if there are more than one molecular
     * index. If fewer tags are given than molecular indexes found, then the remaining molecular indexes will be concatenated and stored
     * in the last tag.  If more tags are provided than molecular indexes found, the additional tags will not be used.
     */
    public ClusterDataToSamConverter withTagPerMolecularIndex(final List<String> tagPerMolecularIndex) {
        if (tagPerMolecularIndex == null) {
            throw new IllegalArgumentException("Null given for tagPerMolecularIndex");
        }
        this.tagPerMolecularIndex = tagPerMolecularIndex;
        return this;
    }

    /**
     * Creates a new SAM record from the basecall data
     */
    private SAMRecord createSamRecord(final ReadData readData, final String readName, final boolean isPf, final boolean firstOfPair,
                                      final String unmatchedBarcode, final String barcodeQuality,
                                      final List<String> molecularIndexes, final List<String> molecularIndexQualities) {
        final SAMRecord sam = new SAMRecord(null);
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
            sam.setAttribute(SAMTag.RG.name(), readGroupId);
        }

        // If it's a barcoded run and it has been decided that the original BC value should be added to the record, do it
        if (unmatchedBarcode != null) {
            sam.setAttribute(SAMTag.BC.name(), unmatchedBarcode);
            if (barcodeQuality != null ) {
                sam.setAttribute(SAMTag.QT.name(), barcodeQuality);
            }
        }

        if (!molecularIndexes.isEmpty()) {
            if (!this.MOLECULAR_INDEX_TAG.isEmpty()) {
                sam.setAttribute(this.MOLECULAR_INDEX_TAG, String.join(MOLECULAR_INDEX_, molecularIndexes));
            }
            if (!this.MOLECULAR_INDEX_QUALITY_TAG.isEmpty()) {
                sam.setAttribute(this.MOLECULAR_INDEX_QUALITY_TAG, String.join(MOLECULAR_INDEX_, molecularIndexQualities));
            }
            if (!this.tagPerMolecularIndex.isEmpty()) {
                if (tagPerMolecularIndex.size() != molecularIndexes.size()) {
                    throw new PicardException("Found " + molecularIndexes.size() + " molecular indexes but only " + tagPerMolecularIndex.size() + " SAM tags given.");
                }
                for (int i = 0; i < this.tagPerMolecularIndex.size(); i++) {
                    sam.setAttribute(this.tagPerMolecularIndex.get(i), molecularIndexes.get(i));
                }
            }
        }

        return sam;
    }

    /**
     * Creates the SAMRecord for each read in the cluster
     */
    public IlluminaBasecallsToSam.SAMRecordsForCluster convertClusterToOutputRecord(final ClusterData cluster) {

        final IlluminaBasecallsToSam.SAMRecordsForCluster ret = new IlluminaBasecallsToSam.SAMRecordsForCluster(outputRecordsPerCluster);
        final String readName = readNameEncoder.generateReadName(cluster, null); // Use null here to prevent /1 or /2 suffixes on read name.

        // Get and transform the unmatched barcode, if any, to store with the reads
        String unmatchedBarcode = null;
        if (hasSampleBarcode && (this.barcodePopulationStrategy == PopulateBarcode.ALWAYS ||
                this.barcodePopulationStrategy == PopulateBarcode.ORPHANS_ONLY &&
                        cluster.getMatchedBarcode() == null ||
                this.barcodePopulationStrategy == PopulateBarcode.INEXACT_MATCH &&
                        !IlluminaUtil.byteArrayToString(getBarcodeSeqs(cluster),"").equals(cluster.getMatchedBarcode()))) {
            unmatchedBarcode = getUnmatchedBarcode(cluster);
        }

        String barcodeQuality = null;
        if (unmatchedBarcode != null && includeQualitiesWithBarcode) {
            barcodeQuality = getBarcodeQuality(cluster);
        }


        final List<String> molecularIndexes = new ArrayList<>();
        final List<String> molecularIndexQualities = new ArrayList<>();
        if (hasMolecularBarcode) {
            for (int molecularBarcodeIndice : molecularBarcodeIndices) {
                molecularIndexes.add(convertMissingToNoCall(new String(cluster.getRead(molecularBarcodeIndice).getBases())));
                molecularIndexQualities.add(SAMUtils.phredToFastq(cluster.getRead(molecularBarcodeIndice).getQualities()));
            }
        }

        final SAMRecord firstOfPair = createSamRecord(
                cluster.getRead(templateIndices[0]),
                readName,
                cluster.isPf(),
                true,
                unmatchedBarcode,
                barcodeQuality,
                molecularIndexes,
                molecularIndexQualities);
        ret.records[0] = firstOfPair;


        SAMRecord secondOfPair = null;

        if (isPairedEnd) {
            secondOfPair = createSamRecord(
                    cluster.getRead(templateIndices[1]),
                    readName,
                    cluster.isPf(),
                    false,
                    unmatchedBarcode,
                    barcodeQuality,
                    molecularIndexes,
                    molecularIndexQualities);
            ret.records[1] = secondOfPair;
        }

        if (adapterMarker != null) {
            // Clip the read
            if (isPairedEnd) {
                adapterMarker.adapterTrimIlluminaPairedReads(firstOfPair, secondOfPair);
            } else {
                adapterMarker.adapterTrimIlluminaSingleRead(firstOfPair);
            }
        }
        return ret;
    }

    private String getBarcodeQuality(ClusterData cluster) {
        final StringJoiner barcodeQ = new StringJoiner(MOLECULAR_INDEX_QUALITY_DELIMITER);

        for (int sampleBarcodeIndex : sampleBarcodeIndices) {
            barcodeQ.add(SAMUtils.phredToFastq(cluster.getRead(sampleBarcodeIndex).getQualities()));
        }
        return barcodeQ.toString();
    }

    private String getUnmatchedBarcode(ClusterData cluster) {
        return convertMissingToNoCall(IlluminaUtil.barcodeSeqsToString(getBarcodeSeqs(cluster)));
    }

    private byte[][] getBarcodeSeqs(ClusterData cluster) {
        final byte[][] barcode = new byte[sampleBarcodeIndices.length][];
        for (int i = 0; i < sampleBarcodeIndices.length; i++) {
            barcode[i] = cluster.getRead(sampleBarcodeIndices[i]).getBases();
        }
        return barcode;
    }

    private static String convertMissingToNoCall(final String barcode){
        return barcode.replace(MISSING_BARCODE, MISSING_BARCODE_BASE);
    }
}
