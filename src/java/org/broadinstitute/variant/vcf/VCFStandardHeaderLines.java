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

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broad.tribble.TribbleException;
import org.broadinstitute.variant.utils.GeneralUtils;

import java.util.*;

/**
 * Manages header lines for standard VCF INFO and FORMAT fields
 *
 * Provides simple mechanisms for registering standard lines,
 * looking them up, and adding them to headers
 *
 * @author Mark DePristo
 * @since 6/12
 */
public class VCFStandardHeaderLines {
    /**
     * Enabling this causes us to repair header lines even if only their descriptions differ
     */
    private final static boolean REPAIR_BAD_DESCRIPTIONS = false;
    private static Standards<VCFFormatHeaderLine> formatStandards = new Standards<VCFFormatHeaderLine>();
    private static Standards<VCFInfoHeaderLine> infoStandards = new Standards<VCFInfoHeaderLine>();

    /**
     * Walks over the VCF header and repairs the standard VCF header lines in it, returning a freshly
     * allocated VCFHeader with standard VCF header lines repaired as necessary
     *
     * @param header
     * @return
     */
    @Requires("header != null")
    @Ensures("result != null")
    public static VCFHeader repairStandardHeaderLines(final VCFHeader header) {
        final Set<VCFHeaderLine> newLines = new LinkedHashSet<VCFHeaderLine>(header.getMetaDataInInputOrder().size());
        for ( VCFHeaderLine line : header.getMetaDataInInputOrder() ) {
            if ( line instanceof VCFFormatHeaderLine ) {
                line = formatStandards.repair((VCFFormatHeaderLine) line);
            } else if ( line instanceof VCFInfoHeaderLine) {
                line = infoStandards.repair((VCFInfoHeaderLine) line);
            }

            newLines.add(line);
        }

        return new VCFHeader(newLines, header.getGenotypeSamples());
    }

    /**
     * Adds header lines for each of the format fields in IDs to header, returning the set of
     * IDs without standard descriptions, unless throwErrorForMissing is true, in which
     * case this situation results in a TribbleException
     *
     * @param IDs
     * @return
     */
    public static Set<String> addStandardFormatLines(final Set<VCFHeaderLine> headerLines, final boolean throwErrorForMissing, final Collection<String> IDs) {
        return formatStandards.addToHeader(headerLines, IDs, throwErrorForMissing);
    }

    /**
     * @see #addStandardFormatLines(java.util.Set, boolean, java.util.Collection)
     *
     * @param headerLines
     * @param throwErrorForMissing
     * @param IDs
     * @return
     */
    public static Set<String> addStandardFormatLines(final Set<VCFHeaderLine> headerLines, final boolean throwErrorForMissing, final String ... IDs) {
        return addStandardFormatLines(headerLines, throwErrorForMissing, Arrays.asList(IDs));
    }

    /**
     * Returns the standard format line for ID.  If none exists, return null or throw an exception, depending
     * on throwErrorForMissing
     *
     * @param ID
     * @param throwErrorForMissing
     * @return
     */
    public static VCFFormatHeaderLine getFormatLine(final String ID, final boolean throwErrorForMissing) {
        return formatStandards.get(ID, throwErrorForMissing);
    }

    /**
     * Returns the standard format line for ID.  If none exists throw an exception
     *
     * @param ID
     * @return
     */
    public static VCFFormatHeaderLine getFormatLine(final String ID) {
        return formatStandards.get(ID, true);
    }

    private static void registerStandard(final VCFFormatHeaderLine line) {
        formatStandards.add(line);
    }

    /**
     * Adds header lines for each of the info fields in IDs to header, returning the set of
     * IDs without standard descriptions, unless throwErrorForMissing is true, in which
     * case this situation results in a TribbleException
     *
     * @param IDs
     * @return
     */
    public static Set<String> addStandardInfoLines(final Set<VCFHeaderLine> headerLines, final boolean throwErrorForMissing, final Collection<String> IDs) {
        return infoStandards.addToHeader(headerLines, IDs, throwErrorForMissing);
    }

    /**
     * @see #addStandardFormatLines(java.util.Set, boolean, java.util.Collection)
     *
     * @param IDs
     * @return
     */
    public static Set<String> addStandardInfoLines(final Set<VCFHeaderLine> headerLines, final boolean throwErrorForMissing, final String ... IDs) {
        return addStandardInfoLines(headerLines, throwErrorForMissing, Arrays.asList(IDs));
    }

    /**
     * Returns the standard info line for ID.  If none exists, return null or throw an exception, depending
     * on throwErrorForMissing
     *
     * @param ID
     * @param throwErrorForMissing
     * @return
     */
    public static VCFInfoHeaderLine getInfoLine(final String ID, final boolean throwErrorForMissing) {
        return infoStandards.get(ID, throwErrorForMissing);
    }

    /**
     * Returns the standard info line for ID.  If none exists throw an exception
     *
     * @param ID
     * @return
     */
    public static VCFInfoHeaderLine getInfoLine(final String ID) {
        return getInfoLine(ID, true);
    }

    private static void registerStandard(final VCFInfoHeaderLine line) {
        infoStandards.add(line);
    }


    //
    // VCF header line constants
    //
    static {
        // FORMAT lines
        registerStandard(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1, VCFHeaderLineType.String, "Genotype"));
        registerStandard(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_QUALITY_KEY, 1, VCFHeaderLineType.Integer, "Genotype Quality"));
        registerStandard(new VCFFormatHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Approximate read depth (reads with MQ=255 or with bad mates are filtered)"));
        registerStandard(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_PL_KEY, VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification"));
        registerStandard(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Allelic depths for the ref and alt alleles in the order listed"));
        registerStandard(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_FILTER_KEY, 1, VCFHeaderLineType.String, "Genotype-level filter"));

        // INFO lines
        registerStandard(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1, VCFHeaderLineType.Integer, "Stop position of the interval"));
        registerStandard(new VCFInfoHeaderLine(VCFConstants.MLE_ALLELE_COUNT_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed"));
        registerStandard(new VCFInfoHeaderLine(VCFConstants.MLE_ALLELE_FREQUENCY_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed"));
        registerStandard(new VCFInfoHeaderLine(VCFConstants.DOWNSAMPLED_KEY, 0, VCFHeaderLineType.Flag, "Were any of the samples downsampled?"));
        registerStandard(new VCFInfoHeaderLine(VCFConstants.DBSNP_KEY, 0, VCFHeaderLineType.Flag, "dbSNP Membership"));
        registerStandard(new VCFInfoHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Approximate read depth; some reads may have been filtered"));
        registerStandard(new VCFInfoHeaderLine(VCFConstants.STRAND_BIAS_KEY, 1, VCFHeaderLineType.Float, "Strand Bias"));
        registerStandard(new VCFInfoHeaderLine(VCFConstants.ALLELE_FREQUENCY_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Allele Frequency, for each ALT allele, in the same order as listed"));
        registerStandard(new VCFInfoHeaderLine(VCFConstants.ALLELE_COUNT_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Allele count in genotypes, for each ALT allele, in the same order as listed"));
        registerStandard(new VCFInfoHeaderLine(VCFConstants.ALLELE_NUMBER_KEY, 1, VCFHeaderLineType.Integer, "Total number of alleles in called genotypes"));
        registerStandard(new VCFInfoHeaderLine(VCFConstants.MAPPING_QUALITY_ZERO_KEY, 1, VCFHeaderLineType.Integer, "Total Mapping Quality Zero Reads"));
        registerStandard(new VCFInfoHeaderLine(VCFConstants.RMS_MAPPING_QUALITY_KEY, 1, VCFHeaderLineType.Float, "RMS Mapping Quality"));
        registerStandard(new VCFInfoHeaderLine(VCFConstants.SOMATIC_KEY, 0, VCFHeaderLineType.Flag, "Somatic event"));
    }

    private static class Standards<T extends VCFCompoundHeaderLine> {
        private final Map<String, T> standards = new HashMap<String, T>();

        @Requires("line != null")
        @Ensures({"result != null", "result.getID().equals(line.getID())"})
        public T repair(final T line) {
            final T standard = get(line.getID(), false);
            if ( standard != null ) {
                final boolean badCountType = line.getCountType() != standard.getCountType();
                final boolean badCount = line.isFixedCount() && ! badCountType && line.getCount() != standard.getCount();
                final boolean badType = line.getType() != standard.getType();
                final boolean badDesc = ! line.getDescription().equals(standard.getDescription());
                final boolean needsRepair = badCountType || badCount || badType || (REPAIR_BAD_DESCRIPTIONS && badDesc);

                if ( needsRepair ) {
                    if ( GeneralUtils.DEBUG_MODE_ENABLED ) {
                        System.err.println("Repairing standard header line for field " + line.getID() + " because"
                                           + (badCountType ? " -- count types disagree; header has " + line.getCountType() + " but standard is " + standard.getCountType() : "")
                                           + (badType ? " -- type disagree; header has " + line.getType() + " but standard is " + standard.getType() : "")
                                           + (badCount ? " -- counts disagree; header has " + line.getCount() + " but standard is " + standard.getCount() : "")
                                           + (badDesc ? " -- descriptions disagree; header has '" + line.getDescription() + "' but standard is '" + standard.getDescription() + "'": ""));
                    }
                    return standard;
                } else
                    return line;
            } else
                return line;
        }

        @Requires("headerLines != null")
        @Ensures({"result != null", "result.isEmpty() || ! throwErrorForMissing", "IDs.containsAll(result)"})
        public Set<String> addToHeader(final Set<VCFHeaderLine> headerLines, final Collection<String> IDs, final boolean throwErrorForMissing) {
            final Set<String> missing = new HashSet<String>();
            for ( final String ID : IDs ) {
                final T line = get(ID, throwErrorForMissing);
                if ( line == null )
                    missing.add(ID);
                else
                    headerLines.add(line);
            }

            return missing;
        }

        @Requires("line != null")
        @Ensures({"standards.containsKey(line.getID())"})
        public void add(final T line) {
            if ( standards.containsKey(line.getID()) )
                throw new TribbleException("Attempting to add multiple standard header lines for ID " + line.getID());
            standards.put(line.getID(), line);
        }

        @Requires("ID != null")
        @Ensures({"result != null || ! throwErrorForMissing"})
        public T get(final String ID, final boolean throwErrorForMissing) {
            final T x = standards.get(ID);
            if ( throwErrorForMissing && x == null )
                throw new TribbleException("Couldn't find a standard VCF header line for field " + ID);
            return x;
        }
    }
}
