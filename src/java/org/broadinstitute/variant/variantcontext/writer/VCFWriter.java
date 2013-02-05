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

package org.broadinstitute.variant.variantcontext.writer;

import net.sf.samtools.SAMSequenceDictionary;
import org.broad.tribble.TribbleException;
import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.variant.vcf.*;
import org.broadinstitute.variant.variantcontext.*;

import java.io.*;
import java.lang.reflect.Array;
import java.nio.charset.Charset;
import java.util.*;

/**
 * this class writes VCF files
 */
class VCFWriter extends IndexingVariantContextWriter {
    private final static String VERSION_LINE = VCFHeader.METADATA_INDICATOR + VCFHeaderVersion.VCF4_1.getFormatString() + "=" + VCFHeaderVersion.VCF4_1.getVersionString();

    // should we write genotypes or just sites?
    final protected boolean doNotWriteGenotypes;

    // the VCF header we're storing
    protected VCFHeader mHeader = null;

    final private boolean allowMissingFieldsInHeader;

    /**
     * The VCF writer uses an internal Writer, based by the ByteArrayOutputStream lineBuffer,
     * to temp. buffer the header and per-site output before flushing the per line output
     * in one go to the super.getOutputStream.  This results in high-performance, proper encoding,
     * and allows us to avoid flushing explicitly the output stream getOutputStream, which
     * allows us to properly compress vcfs in gz format without breaking indexing on the fly
     * for uncompressed streams.
     */
    private static final int INITIAL_BUFFER_SIZE = 1024 * 16;
    private final ByteArrayOutputStream lineBuffer = new ByteArrayOutputStream(INITIAL_BUFFER_SIZE);
    private final Writer writer;

    /**
     * The encoding used for VCF files.  ISO-8859-1
     */
    final private Charset charset;

    private IntGenotypeFieldAccessors intGenotypeFieldAccessors = new IntGenotypeFieldAccessors();

    public VCFWriter(final File location, final OutputStream output, final SAMSequenceDictionary refDict,
                     final boolean enableOnTheFlyIndexing, boolean doNotWriteGenotypes,
                     final boolean allowMissingFieldsInHeader ) {
        super(writerName(location, output), location, output, refDict, enableOnTheFlyIndexing);
        this.doNotWriteGenotypes = doNotWriteGenotypes;
        this.allowMissingFieldsInHeader = allowMissingFieldsInHeader;
        this.charset = Charset.forName("ISO-8859-1");
        this.writer = new OutputStreamWriter(lineBuffer, charset);
    }

    // --------------------------------------------------------------------------------
    //
    // VCFWriter interface functions
    //
    // --------------------------------------------------------------------------------

    /**
     * Write String s to the internal buffered writer.
     *
     * flushBuffer() must be called to actually write the data to the true output stream.
     *
     * @param s the string to write
     * @throws IOException
     */
    private void write(final String s) throws IOException {
        writer.write(s);
    }

    /**
     * Actually write the line buffer contents to the destination output stream.
     *
     * After calling this function the line buffer is reset, so the contents of the buffer can be reused
     *
     * @throws IOException
     */
    private void flushBuffer() throws IOException {
        writer.flush();
        getOutputStream().write(lineBuffer.toByteArray());
        lineBuffer.reset();
    }

    @Override
    public void writeHeader(VCFHeader header) {
        // note we need to update the mHeader object after this call because they header
        // may have genotypes trimmed out of it, if doNotWriteGenotypes is true
        try {
            mHeader = writeHeader(header, writer, doNotWriteGenotypes, getVersionLine(), getStreamName());
            flushBuffer();
        } catch ( IOException e ) {
            throw new RuntimeException("Couldn't write file " + getStreamName(), e);
        }
    }

    public static String getVersionLine() {
        return VERSION_LINE;
    }

    public static VCFHeader writeHeader(VCFHeader header,
                                        final Writer writer,
                                        final boolean doNotWriteGenotypes,
                                        final String versionLine,
                                        final String streamNameForError) {
        header = doNotWriteGenotypes ? new VCFHeader(header.getMetaDataInSortedOrder()) : header;
        
        try {
            // the file format field needs to be written first
            writer.write(versionLine + "\n");

            for ( VCFHeaderLine line : header.getMetaDataInSortedOrder() ) {
                if ( VCFHeaderVersion.isFormatString(line.getKey()) )
                    continue;

                writer.write(VCFHeader.METADATA_INDICATOR);
                writer.write(line.toString());
                writer.write("\n");
            }

            // write out the column line
            writer.write(VCFHeader.HEADER_INDICATOR);
            boolean isFirst = true;
            for ( VCFHeader.HEADER_FIELDS field : header.getHeaderFields() ) {
                if ( isFirst )
                    isFirst = false; // don't write out a field separator
                else
                    writer.write(VCFConstants.FIELD_SEPARATOR);
                writer.write(field.toString());
            }

            if ( header.hasGenotypingData() ) {
                writer.write(VCFConstants.FIELD_SEPARATOR);
                writer.write("FORMAT");
                for ( String sample : header.getGenotypeSamples() ) {
                    writer.write(VCFConstants.FIELD_SEPARATOR);
                    writer.write(sample);
                }
            }

            writer.write("\n");
            writer.flush();  // necessary so that writing to an output stream will work
        }
        catch (IOException e) {
            throw new RuntimeException("IOException writing the VCF header to " + streamNameForError, e);
        }

        return header;
    }

    /**
     * attempt to close the VCF file
     */
    @Override
    public void close() {
        // try to close the vcf stream
        try {
            // TODO -- would it be useful to null out the line buffer so we don't have it around unnecessarily?
            writer.close();
        } catch (IOException e) {
            throw new RuntimeException("Unable to close " + getStreamName(), e);
        }

        super.close();
    }

    /**
     * add a record to the file
     *
     * @param vc      the Variant Context object
     */
    @Override
    public void add(VariantContext vc) {
        if ( mHeader == null )
            throw new IllegalStateException("The VCF Header must be written before records can be added: " + getStreamName());

        if ( doNotWriteGenotypes )
            vc = new VariantContextBuilder(vc).noGenotypes().make();

        try {
            super.add(vc);

            Map<Allele, String> alleleMap = buildAlleleMap(vc);

            // CHROM
            write(vc.getChr());
            write(VCFConstants.FIELD_SEPARATOR);

            // POS
            write(String.valueOf(vc.getStart()));
            write(VCFConstants.FIELD_SEPARATOR);

            // ID
            String ID = vc.getID();
            write(ID);
            write(VCFConstants.FIELD_SEPARATOR);

            // REF
            String refString = vc.getReference().getDisplayString();
            write(refString);
            write(VCFConstants.FIELD_SEPARATOR);

            // ALT
            if ( vc.isVariant() ) {
                Allele altAllele = vc.getAlternateAllele(0);
                String alt = altAllele.getDisplayString();
                write(alt);

                for (int i = 1; i < vc.getAlternateAlleles().size(); i++) {
                    altAllele = vc.getAlternateAllele(i);
                    alt = altAllele.getDisplayString();
                    write(",");
                    write(alt);
                }
            } else {
                write(VCFConstants.EMPTY_ALTERNATE_ALLELE_FIELD);
            }
            write(VCFConstants.FIELD_SEPARATOR);

            // QUAL
            if ( !vc.hasLog10PError() )
                write(VCFConstants.MISSING_VALUE_v4);
            else
                write(formatQualValue(vc.getPhredScaledQual()));
            write(VCFConstants.FIELD_SEPARATOR);

            // FILTER
            String filters = getFilterString(vc);
            write(filters);
            write(VCFConstants.FIELD_SEPARATOR);

            // INFO
            Map<String, String> infoFields = new TreeMap<String, String>();
            for ( Map.Entry<String, Object> field : vc.getAttributes().entrySet() ) {
                String key = field.getKey();

                if ( ! mHeader.hasInfoLine(key) )
                    fieldIsMissingFromHeaderError(vc, key, "INFO");

                String outputValue = formatVCFField(field.getValue());
                if ( outputValue != null )
                    infoFields.put(key, outputValue);
            }
            writeInfoString(infoFields);

            // FORMAT
            final GenotypesContext gc = vc.getGenotypes();
            if ( gc.isLazyWithData() && ((LazyGenotypesContext)gc).getUnparsedGenotypeData() instanceof String ) {
                write(VCFConstants.FIELD_SEPARATOR);
                write(((LazyGenotypesContext) gc).getUnparsedGenotypeData().toString());
            } else {
                List<String> genotypeAttributeKeys = calcVCFGenotypeKeys(vc, mHeader);
                if ( ! genotypeAttributeKeys.isEmpty() ) {
                    for ( final String format : genotypeAttributeKeys )
                        if ( ! mHeader.hasFormatLine(format) )
                            fieldIsMissingFromHeaderError(vc, format, "FORMAT");

                    final String genotypeFormatString = ParsingUtils.join(VCFConstants.GENOTYPE_FIELD_SEPARATOR, genotypeAttributeKeys);

                    write(VCFConstants.FIELD_SEPARATOR);
                    write(genotypeFormatString);

                    addGenotypeData(vc, alleleMap, genotypeAttributeKeys);
                }
            }
            
            write("\n");
            // note that we cannot call flush here if we want block gzipping to work properly
            // calling flush results in all gzipped blocks for each variant
            flushBuffer();
        } catch (IOException e) {
            throw new RuntimeException("Unable to write the VCF object to " + getStreamName(), e);
        }
    }

    private static Map<Allele, String> buildAlleleMap(final VariantContext vc) {
        final Map<Allele, String> alleleMap = new HashMap<Allele, String>(vc.getAlleles().size()+1);
        alleleMap.put(Allele.NO_CALL, VCFConstants.EMPTY_ALLELE); // convenience for lookup

        final List<Allele> alleles = vc.getAlleles();
        for ( int i = 0; i < alleles.size(); i++ ) {
            alleleMap.put(alleles.get(i), String.valueOf(i));
        }

        return alleleMap;
    }

    // --------------------------------------------------------------------------------
    //
    // implementation functions
    //
    // --------------------------------------------------------------------------------

    private final String getFilterString(final VariantContext vc) {
        if ( vc.isFiltered() ) {
            for ( final String filter : vc.getFilters() )
                if ( ! mHeader.hasFilterLine(filter) )
                    fieldIsMissingFromHeaderError(vc, filter, "FILTER");

            return ParsingUtils.join(";", ParsingUtils.sortList(vc.getFilters()));
        }
        else if ( vc.filtersWereApplied() )
            return VCFConstants.PASSES_FILTERS_v4;
        else
            return VCFConstants.UNFILTERED;
    }

    private static final String QUAL_FORMAT_STRING = "%.2f";
    private static final String QUAL_FORMAT_EXTENSION_TO_TRIM = ".00";

    private String formatQualValue(double qual) {
        String s = String.format(QUAL_FORMAT_STRING, qual);
        if ( s.endsWith(QUAL_FORMAT_EXTENSION_TO_TRIM) )
            s = s.substring(0, s.length() - QUAL_FORMAT_EXTENSION_TO_TRIM.length());
        return s;
    }

    /**
     * create the info string; assumes that no values are null
     *
     * @param infoFields a map of info fields
     * @throws IOException for writer
     */
    private void writeInfoString(Map<String, String> infoFields) throws IOException {
        if ( infoFields.isEmpty() ) {
            write(VCFConstants.EMPTY_INFO_FIELD);
            return;
        }

        boolean isFirst = true;
        for ( Map.Entry<String, String> entry : infoFields.entrySet() ) {
            if ( isFirst )
                isFirst = false;
            else
                write(VCFConstants.INFO_FIELD_SEPARATOR);

            String key = entry.getKey();
            write(key);

            if ( !entry.getValue().equals("") ) {
                VCFInfoHeaderLine metaData = mHeader.getInfoHeaderLine(key);
                if ( metaData == null || metaData.getCountType() != VCFHeaderLineCount.INTEGER || metaData.getCount() != 0 ) {
                    write("=");
                    write(entry.getValue());
                }
            }
        }
    }

    /**
     * add the genotype data
     *
     * @param vc                     the variant context
     * @param genotypeFormatKeys  Genotype formatting string
     * @param alleleMap              alleles for this context
     * @throws IOException for writer
     */
    private void addGenotypeData(VariantContext vc, Map<Allele, String> alleleMap, List<String> genotypeFormatKeys)
    throws IOException {
        final int ploidy = vc.getMaxPloidy(2);

        for ( String sample : mHeader.getGenotypeSamples() ) {
            write(VCFConstants.FIELD_SEPARATOR);

            Genotype g = vc.getGenotype(sample);
            if ( g == null ) g = GenotypeBuilder.createMissing(sample, ploidy);

            final List<String> attrs = new ArrayList<String>(genotypeFormatKeys.size());
            for ( String field : genotypeFormatKeys ) {
                if ( field.equals(VCFConstants.GENOTYPE_KEY) ) {
                    if ( !g.isAvailable() ) {
                        throw new IllegalStateException("GTs cannot be missing for some samples if they are available for others in the record");
                    }

                    writeAllele(g.getAllele(0), alleleMap);
                    for (int i = 1; i < g.getPloidy(); i++) {
                        write(g.isPhased() ? VCFConstants.PHASED : VCFConstants.UNPHASED);
                        writeAllele(g.getAllele(i), alleleMap);
                    }

                    continue;
                } else {
                    String outputValue;
                    if ( field.equals(VCFConstants.GENOTYPE_FILTER_KEY ) ) {
                        outputValue = g.isFiltered() ? g.getFilters() : VCFConstants.PASSES_FILTERS_v4;
                    } else {
                        final IntGenotypeFieldAccessors.Accessor accessor = intGenotypeFieldAccessors.getAccessor(field);
                        if ( accessor != null ) {
                            final int[] intValues = accessor.getValues(g);
                            if ( intValues == null )
                                outputValue = VCFConstants.MISSING_VALUE_v4;
                            else if ( intValues.length == 1 ) // fast path
                                outputValue = Integer.toString(intValues[0]);
                            else {
                                StringBuilder sb = new StringBuilder();
                                sb.append(intValues[0]);
                                for ( int i = 1; i < intValues.length; i++) {
                                    sb.append(",");
                                    sb.append(intValues[i]);
                                }
                                outputValue = sb.toString();
                            }
                        } else {
                            Object val = g.hasExtendedAttribute(field) ? g.getExtendedAttribute(field) : VCFConstants.MISSING_VALUE_v4;

                            VCFFormatHeaderLine metaData = mHeader.getFormatHeaderLine(field);
                            if ( metaData != null ) {
                                int numInFormatField = metaData.getCount(vc);
                                if ( numInFormatField > 1 && val.equals(VCFConstants.MISSING_VALUE_v4) ) {
                                    // If we have a missing field but multiple values are expected, we need to construct a new string with all fields.
                                    // For example, if Number=2, the string has to be ".,."
                                    StringBuilder sb = new StringBuilder(VCFConstants.MISSING_VALUE_v4);
                                    for ( int i = 1; i < numInFormatField; i++ ) {
                                        sb.append(",");
                                        sb.append(VCFConstants.MISSING_VALUE_v4);
                                    }
                                    val = sb.toString();
                                }
                            }

                            // assume that if key is absent, then the given string encoding suffices
                            outputValue = formatVCFField(val);
                        }
                    }

                    if ( outputValue != null )
                        attrs.add(outputValue);
                }
            }

            // strip off trailing missing values
            for (int i = attrs.size()-1; i >= 0; i--) {
                if ( isMissingValue(attrs.get(i)) )
                    attrs.remove(i);
                else
                    break;
            }

            for (int i = 0; i < attrs.size(); i++) {
                if ( i > 0 || genotypeFormatKeys.contains(VCFConstants.GENOTYPE_KEY) )
                    write(VCFConstants.GENOTYPE_FIELD_SEPARATOR);
                write(attrs.get(i));
            }
        }
    }

    private boolean isMissingValue(String s) {
        // we need to deal with the case that it's a list of missing values
        return (countOccurrences(VCFConstants.MISSING_VALUE_v4.charAt(0), s) + countOccurrences(',', s) == s.length());
    }

    private void writeAllele(Allele allele, Map<Allele, String> alleleMap) throws IOException {
        String encoding = alleleMap.get(allele);
        if ( encoding == null )
            throw new TribbleException.InternalCodecException("Allele " + allele + " is not an allele in the variant context");
        write(encoding);
    }

    /**
     * Takes a double value and pretty prints it to a String for display
     *
     * Large doubles => gets %.2f style formatting
     * Doubles < 1 / 10 but > 1/100 </>=> get %.3f style formatting
     * Double < 1/100 => %.3e formatting
     * @param d
     * @return
     */
    public static final String formatVCFDouble(final double d) {
        String format;
        if ( d < 1 ) {
            if ( d < 0.01 ) {
                if ( Math.abs(d) >= 1e-20 )
                    format = "%.3e";
                else {
                    // return a zero format
                    return "0.00";
                }
            } else {
                format = "%.3f";
            }
        } else {
            format = "%.2f";
        }

        return String.format(format, d);
    }

    public static String formatVCFField(Object val) {
        String result;
        if ( val == null )
            result = VCFConstants.MISSING_VALUE_v4;
        else if ( val instanceof Double )
            result = formatVCFDouble((Double) val);
        else if ( val instanceof Boolean )
            result = (Boolean)val ? "" : null; // empty string for true, null for false
        else if ( val instanceof List ) {
            result = formatVCFField(((List)val).toArray());
        } else if ( val.getClass().isArray() ) {
            final int length = Array.getLength(val);
            if ( length == 0 )
                return formatVCFField(null);
            final StringBuilder sb = new StringBuilder(formatVCFField(Array.get(val, 0)));
            for ( int i = 1; i < length; i++) {
                sb.append(",");
                sb.append(formatVCFField(Array.get(val, i)));
            }
            result = sb.toString();
        } else
            result = val.toString();

        return result;
    }

    /**
     * Determine which genotype fields are in use in the genotypes in VC
     * @param vc
     * @return an ordered list of genotype fields in use in VC.  If vc has genotypes this will always include GT first
     */
    public static List<String> calcVCFGenotypeKeys(final VariantContext vc, final VCFHeader header) {
        Set<String> keys = new HashSet<String>();

        boolean sawGoodGT = false;
        boolean sawGoodQual = false;
        boolean sawGenotypeFilter = false;
        boolean sawDP = false;
        boolean sawAD = false;
        boolean sawPL = false;
        for ( final Genotype g : vc.getGenotypes() ) {
            keys.addAll(g.getExtendedAttributes().keySet());
            if ( g.isAvailable() ) sawGoodGT = true;
            if ( g.hasGQ() ) sawGoodQual = true;
            if ( g.hasDP() ) sawDP = true;
            if ( g.hasAD() ) sawAD = true;
            if ( g.hasPL() ) sawPL = true;
            if (g.isFiltered()) sawGenotypeFilter = true;
        }

        if ( sawGoodQual ) keys.add(VCFConstants.GENOTYPE_QUALITY_KEY);
        if ( sawDP ) keys.add(VCFConstants.DEPTH_KEY);
        if ( sawAD ) keys.add(VCFConstants.GENOTYPE_ALLELE_DEPTHS);
        if ( sawPL ) keys.add(VCFConstants.GENOTYPE_PL_KEY);
        if ( sawGenotypeFilter ) keys.add(VCFConstants.GENOTYPE_FILTER_KEY);

        List<String> sortedList = ParsingUtils.sortList(new ArrayList<String>(keys));

        // make sure the GT is first
        if ( sawGoodGT ) {
            List<String> newList = new ArrayList<String>(sortedList.size()+1);
            newList.add(VCFConstants.GENOTYPE_KEY);
            newList.addAll(sortedList);
            sortedList = newList;
        }

        if ( sortedList.isEmpty() && header.hasGenotypingData() ) {
            // this needs to be done in case all samples are no-calls
            return Collections.singletonList(VCFConstants.GENOTYPE_KEY);
        } else {
            return sortedList;
        }
    }


    private static int countOccurrences(char c, String s) {
           int count = 0;
           for (int i = 0; i < s.length(); i++) {
               count += s.charAt(i) == c ? 1 : 0;
           }
           return count;
    }

    private final void fieldIsMissingFromHeaderError(final VariantContext vc, final String id, final String field) {
        if ( !allowMissingFieldsInHeader)
            throw new IllegalStateException("Key " + id + " found in VariantContext field " + field
                    + " at " + vc.getChr() + ":" + vc.getStart()
                    + " but this key isn't defined in the VCFHeader.  We require all VCFs to have"
                    + " complete VCF headers by default.");
    }
}
