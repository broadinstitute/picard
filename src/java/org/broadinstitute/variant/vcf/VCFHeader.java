/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.variant.vcf;

import org.broad.tribble.TribbleException;
import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.variant.utils.GeneralUtils;

import java.util.*;


/**
 * NOTE: This class allows duplicate entries in the metadata & stores header lines in
 * lots of places. The original author noted that this should be cleaned up at some point
 * in the future (jgentry - 5/2013)
 *
 * @author aaron
 *         <p/>
 *         Class VCFHeader
 *         <p/>
 *         A class representing the VCF header
 */
public class VCFHeader {

    // the mandatory header fields
    public enum HEADER_FIELDS {
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
    }

    // the associated meta data
    private final Set<VCFHeaderLine> mMetaData = new LinkedHashSet<VCFHeaderLine>();
    private final Map<String, VCFInfoHeaderLine> mInfoMetaData = new LinkedHashMap<String, VCFInfoHeaderLine>();
    private final Map<String, VCFFormatHeaderLine> mFormatMetaData = new LinkedHashMap<String, VCFFormatHeaderLine>();
    private final Map<String, VCFFilterHeaderLine> mFilterMetaData = new LinkedHashMap<String, VCFFilterHeaderLine>();
    private final Map<String, VCFHeaderLine> mOtherMetaData = new LinkedHashMap<String, VCFHeaderLine>();
    private final List<VCFContigHeaderLine> contigMetaData = new ArrayList<VCFContigHeaderLine>();

    // the list of auxillary tags
    private final List<String> mGenotypeSampleNames = new ArrayList<String>();

    // the character string that indicates meta data
    public static final String METADATA_INDICATOR = "##";

    // the header string indicator
    public static final String HEADER_INDICATOR = "#";

    public static final String SOURCE_KEY = "source";
    public static final String REFERENCE_KEY = "reference";
    public static final String CONTIG_KEY = "contig";
    public static final String INTERVALS_KEY = "intervals";
    public static final String EXCLUDE_INTERVALS_KEY = "excludeIntervals";
    public static final String INTERVAL_MERGING_KEY = "interval_merging";
    public static final String INTERVAL_SET_RULE_KEY = "interval_set_rule";
    public static final String INTERVAL_PADDING_KEY = "interval_padding";

    // were the input samples sorted originally (or are we sorting them)?
    private boolean samplesWereAlreadySorted = true;

    // cache for efficient conversion of VCF -> VariantContext
    private ArrayList<String> sampleNamesInOrder = null;
    private HashMap<String, Integer> sampleNameToOffset = null;

    private boolean writeEngineHeaders = true;
    private boolean writeCommandLine = true;

    /**
     * Create an empty VCF header with no header lines and no samples
     */
    public VCFHeader() {
        this(Collections.<VCFHeaderLine>emptySet(), Collections.<String>emptySet());
    }

    /**
     * create a VCF header, given a list of meta data and auxillary tags
     *
     * @param metaData     the meta data associated with this header
     */
    public VCFHeader(final Set<VCFHeaderLine> metaData) {
        mMetaData.addAll(metaData);
        loadVCFVersion();
        loadMetaDataMaps();
    }

    /**
     * Creates a shallow copy of the meta data in VCF header toCopy
     *
     * @param toCopy
     */
    public VCFHeader(final VCFHeader toCopy) {
        this(toCopy.mMetaData);
    }

    /**
     * create a VCF header, given a list of meta data and auxillary tags
     *
     * @param metaData            the meta data associated with this header
     * @param genotypeSampleNames the sample names
     */
    public VCFHeader(final Set<VCFHeaderLine> metaData, final Set<String> genotypeSampleNames) {
        this(metaData, new ArrayList<String>(genotypeSampleNames));
    }

    public VCFHeader(final Set<VCFHeaderLine> metaData, final List<String> genotypeSampleNames) {
        this(metaData);

        if ( genotypeSampleNames.size() != new HashSet<String>(genotypeSampleNames).size() )
            throw new TribbleException.InvalidHeader("BUG: VCF header has duplicate sample names");

        mGenotypeSampleNames.addAll(genotypeSampleNames);
        samplesWereAlreadySorted = ParsingUtils.isSorted(genotypeSampleNames);
        buildVCFReaderMaps(genotypeSampleNames);
    }

    /**
     * Tell this VCF header to use pre-calculated sample name ordering and the
     * sample name -> offset map.  This assumes that all VariantContext created
     * using this header (i.e., read by the VCFCodec) will have genotypes
     * occurring in the same order
     *
     * @param genotypeSampleNamesInAppearenceOrder genotype sample names, must iterator in order of appearance
     */
    private void buildVCFReaderMaps(final Collection<String> genotypeSampleNamesInAppearenceOrder) {
        sampleNamesInOrder = new ArrayList<String>(genotypeSampleNamesInAppearenceOrder.size());
        sampleNameToOffset = new HashMap<String, Integer>(genotypeSampleNamesInAppearenceOrder.size());

        int i = 0;
        for (final String name : genotypeSampleNamesInAppearenceOrder) {
            sampleNamesInOrder.add(name);
            sampleNameToOffset.put(name, i++);
        }
        Collections.sort(sampleNamesInOrder);
    }


    /**
     * Adds a header line to the header metadata.
     *
     * @param headerLine Line to add to the existing metadata component.
     */
    public void addMetaDataLine(final VCFHeaderLine headerLine) {
        mMetaData.add(headerLine);
        loadMetaDataMaps();
    }

    /**
     * @return all of the VCF header lines of the ##contig form in order, or an empty list if none were present
     */
    public List<VCFContigHeaderLine> getContigLines() {
        return Collections.unmodifiableList(contigMetaData);
    }


    /**
     * @return all of the VCF FILTER lines in their original file order, or an empty list if none were present
     */
    public List<VCFFilterHeaderLine> getFilterLines() {
        final List<VCFFilterHeaderLine> filters = new ArrayList<VCFFilterHeaderLine>();
        for (final VCFHeaderLine line : mMetaData) {
            if ( line instanceof VCFFilterHeaderLine )  {
                filters.add((VCFFilterHeaderLine)line);
            }
        }
        return filters;
    }

    /**
     * @return all of the VCF FILTER lines in their original file order, or an empty list if none were present
     */
    public List<VCFIDHeaderLine> getIDHeaderLines() {
        final List<VCFIDHeaderLine> filters = new ArrayList<VCFIDHeaderLine>();
        for (final VCFHeaderLine line : mMetaData) {
            if (line instanceof VCFIDHeaderLine)  {
                filters.add((VCFIDHeaderLine)line);
            }
        }
        return filters;
    }

    /**
     * check our metadata for a VCF version tag, and throw an exception if the version is out of date
     * or the version is not present
     */
    public void loadVCFVersion() {
        final List<VCFHeaderLine> toRemove = new ArrayList<VCFHeaderLine>();
        for (final VCFHeaderLine line : mMetaData)
            if (VCFHeaderVersion.isFormatString(line.getKey())) {
                toRemove.add(line);
            }
        // remove old header lines for now,
        mMetaData.removeAll(toRemove);

    }

    /**
     * load the format/info meta data maps (these are used for quick lookup by key name)
     */
    private void loadMetaDataMaps() {
        for (final VCFHeaderLine line : mMetaData) {
            if ( line instanceof VCFInfoHeaderLine )  {
                final VCFInfoHeaderLine infoLine = (VCFInfoHeaderLine)line;
                addMetaDataMapBinding(mInfoMetaData, infoLine);
            } else if ( line instanceof VCFFormatHeaderLine ) {
                final VCFFormatHeaderLine formatLine = (VCFFormatHeaderLine)line;
                addMetaDataMapBinding(mFormatMetaData, formatLine);
            } else if ( line instanceof VCFFilterHeaderLine ) {
                final VCFFilterHeaderLine filterLine = (VCFFilterHeaderLine)line;
                mFilterMetaData.put(filterLine.getID(), filterLine);
            } else if ( line instanceof VCFContigHeaderLine ) {
                contigMetaData.add((VCFContigHeaderLine)line);
            } else {
                mOtherMetaData.put(line.getKey(), line);
            }
        }

        if ( hasFormatLine(VCFConstants.GENOTYPE_LIKELIHOODS_KEY) && ! hasFormatLine(VCFConstants.GENOTYPE_PL_KEY) ) {
            if ( GeneralUtils.DEBUG_MODE_ENABLED ) {
                System.err.println("Found " + VCFConstants.GENOTYPE_LIKELIHOODS_KEY + " format, but no "
                                   + VCFConstants.GENOTYPE_PL_KEY + " field.  We now only manage PL fields internally"
                                   + " automatically adding a corresponding PL field to your VCF header");
            }
            addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_PL_KEY, VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification"));
        }
    }

    /**
     * Add line to map, issuing warnings about duplicates
     *
     * @param map
     * @param line
     * @param <T>
     */
    private final <T extends VCFCompoundHeaderLine> void addMetaDataMapBinding(final Map<String, T> map, T line) {
        final String key = line.getID();
        if ( map.containsKey(key) ) {
            if ( GeneralUtils.DEBUG_MODE_ENABLED ) {
                System.err.println("Found duplicate VCF header lines for " + key + "; keeping the first only" );
            }
        }
        else {
            map.put(key, line);
        }
    }

    /**
     * get the header fields in order they're presented in the input file (which is now required to be
     * the order presented in the spec).
     *
     * @return a set of the header fields, in order
     */
    public Set<HEADER_FIELDS> getHeaderFields() {
        return new LinkedHashSet<HEADER_FIELDS>(Arrays.asList(HEADER_FIELDS.values()));
    }

    /**
     * get the meta data, associated with this header, in sorted order
     *
     * @return a set of the meta data
     */
    public Set<VCFHeaderLine> getMetaDataInInputOrder() {
        return makeGetMetaDataSet(mMetaData);
    }

    public Set<VCFHeaderLine> getMetaDataInSortedOrder() {
        return makeGetMetaDataSet(new TreeSet<VCFHeaderLine>(mMetaData));
    }

    private static Set<VCFHeaderLine> makeGetMetaDataSet(final Set<VCFHeaderLine> headerLinesInSomeOrder) {
        final Set<VCFHeaderLine> lines = new LinkedHashSet<VCFHeaderLine>();
        lines.add(new VCFHeaderLine(VCFHeaderVersion.VCF4_1.getFormatString(), VCFHeaderVersion.VCF4_1.getVersionString()));
        lines.addAll(headerLinesInSomeOrder);
        return Collections.unmodifiableSet(lines);
    }

    /**
     * Get the VCFHeaderLine whose key equals key.  Returns null if no such line exists
     * @param key
     * @return
     */
    public VCFHeaderLine getMetaDataLine(final String key) {
        for (final VCFHeaderLine line: mMetaData) {
            if ( line.getKey().equals(key) )
                return line;
        }

        return null;
    }

    /**
     * get the genotyping sample names
     *
     * @return a list of the genotype column names, which may be empty if hasGenotypingData() returns false
     */
    public List<String> getGenotypeSamples() {
        return mGenotypeSampleNames;
    }

    public int getNGenotypeSamples() {
        return mGenotypeSampleNames.size();
    }

    /**
     * do we have genotyping data?
     *
     * @return true if we have genotyping columns, false otherwise
     */
    public boolean hasGenotypingData() {
        return getNGenotypeSamples() > 0;
    }

    /**
     * were the input samples sorted originally?
     *
     * @return true if the input samples were sorted originally, false otherwise
     */
    public boolean samplesWereAlreadySorted() {
        return samplesWereAlreadySorted;
    }

    /** @return the column count */
    public int getColumnCount() {
        return HEADER_FIELDS.values().length + (hasGenotypingData() ? mGenotypeSampleNames.size() + 1 : 0);
    }

    /**
     * Returns the INFO HeaderLines in their original ordering
     */
    public Collection<VCFInfoHeaderLine> getInfoHeaderLines() {
        return mInfoMetaData.values();
    }

    /**
     * Returns the FORMAT HeaderLines in their original ordering
     */
    public Collection<VCFFormatHeaderLine> getFormatHeaderLines() {
        return mFormatMetaData.values();
    }

    /**
     * @param id the header key name
     * @return the meta data line, or null if there is none
     */
    public VCFInfoHeaderLine getInfoHeaderLine(final String id) {
        return mInfoMetaData.get(id);
    }

    /**
     * @param id    the header key name
     * @return the meta data line, or null if there is none
     */
    public VCFFormatHeaderLine getFormatHeaderLine(final String id) {
        return mFormatMetaData.get(id);
    }

    /**
     * @param id    the header key name
     * @return the meta data line, or null if there is none
     */
    public VCFFilterHeaderLine getFilterHeaderLine(final String id) {
        return mFilterMetaData.get(id);
    }

    public boolean hasInfoLine(final String id) {
        return getInfoHeaderLine(id) != null;
    }

    public boolean hasFormatLine(final String id) {
        return getFormatHeaderLine(id) != null;
    }

    public boolean hasFilterLine(final String id) {
        return getFilterHeaderLine(id) != null;
    }

    /**
     * @param key    the header key name
     * @return the meta data line, or null if there is none
     */
    public VCFHeaderLine getOtherHeaderLine(final String key) {
        return mOtherMetaData.get(key);
    }

    /**
     * If true additional engine headers will be written to the VCF, otherwise only the walker headers will be output.
     * @return true if additional engine headers will be written to the VCF
     */
    public boolean isWriteEngineHeaders() {
        return writeEngineHeaders;
    }

    /**
     * If true additional engine headers will be written to the VCF, otherwise only the walker headers will be output.
     * @param writeEngineHeaders true if additional engine headers will be written to the VCF
     */
    public void setWriteEngineHeaders(final boolean writeEngineHeaders) {
        this.writeEngineHeaders = writeEngineHeaders;
    }

    /**
     * If true, and isWriteEngineHeaders also returns true, the command line will be written to the VCF.
     * @return true if the command line will be written to the VCF
     */
    public boolean isWriteCommandLine() {
        return writeCommandLine;
    }

    /**
     * If true, and isWriteEngineHeaders also returns true, the command line will be written to the VCF.
     * @param writeCommandLine true if the command line will be written to the VCF
     */
    public void setWriteCommandLine(final boolean writeCommandLine) {
        this.writeCommandLine = writeCommandLine;
    }

    public ArrayList<String> getSampleNamesInOrder() {
        return sampleNamesInOrder;
    }

    public HashMap<String, Integer> getSampleNameToOffset() {
        return sampleNameToOffset;
    }

    @Override
    public String toString() {
        final StringBuilder b = new StringBuilder();
        b.append("[VCFHeader:");
        for ( final VCFHeaderLine line : mMetaData )
            b.append("\n\t").append(line);
        return b.append("\n]").toString();
    }
}
