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

import net.sf.samtools.util.BlockCompressedInputStream;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.NameAwareCodec;
import org.broad.tribble.TribbleException;
import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.variant.utils.GeneralUtils;
import org.broadinstitute.variant.variantcontext.*;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.*;
import java.util.zip.GZIPInputStream;


public abstract class AbstractVCFCodec extends AsciiFeatureCodec<VariantContext> implements NameAwareCodec {
    public final static int MAX_ALLELE_SIZE_BEFORE_WARNING = (int)Math.pow(2, 20);

    protected final static int NUM_STANDARD_FIELDS = 8;  // INFO is the 8th column

    // we have to store the list of strings that make up the header until they're needed
    protected VCFHeader header = null;
    protected VCFHeaderVersion version = null;

    // a mapping of the allele
    protected Map<String, List<Allele>> alleleMap = new HashMap<String, List<Allele>>(3);

    // for ParsingUtils.split
    protected String[] GTValueArray = new String[100];
    protected String[] genotypeKeyArray = new String[100];
    protected String[] infoFieldArray = new String[1000];
    protected String[] infoValueArray = new String[1000];

    // for performance testing purposes
    public static boolean validate = true;

    // a key optimization -- we need a per thread string parts array, so we don't allocate a big array over and over
    // todo: make this thread safe?
    protected String[] parts = null;
    protected String[] genotypeParts = null;
    protected final String[] locParts = new String[6];

    // for performance we cache the hashmap of filter encodings for quick lookup
    protected HashMap<String,List<String>> filterHash = new HashMap<String,List<String>>();

    // we store a name to give to each of the variant contexts we emit
    protected String name = "Unknown";

    protected int lineNo = 0;

    protected Map<String, String> stringCache = new HashMap<String, String>();

    protected boolean warnedAboutNoEqualsForNonFlag = false;

    /**
     * If true, then we'll magically fix up VCF headers on the fly when we read them in
     */
    protected boolean doOnTheFlyModifications = true;

    protected AbstractVCFCodec() {
        super(VariantContext.class);
    }

    /**
     * Creates a LazyParser for a LazyGenotypesContext to use to decode
     * our genotypes only when necessary.  We do this instead of eagarly
     * decoding the genotypes just to turn around and reencode in the frequent
     * case where we don't actually want to manipulate the genotypes
     */
    class LazyVCFGenotypesParser implements LazyGenotypesContext.LazyParser {
        final List<Allele> alleles;
        final String contig;
        final int start;

        LazyVCFGenotypesParser(final List<Allele> alleles, final String contig, final int start) {
            this.alleles = alleles;
            this.contig = contig;
            this.start = start;
        }

        @Override
        public LazyGenotypesContext.LazyData parse(final Object data) {
            //System.out.printf("Loading genotypes... %s:%d%n", contig, start);
            return createGenotypeMap((String) data, alleles, contig, start);
        }
    }

    /**
     * parse the filter string, first checking to see if we already have parsed it in a previous attempt
     * @param filterString the string to parse
     * @return a set of the filters applied
     */
    protected abstract List<String> parseFilters(String filterString);

    /**
     * create a VCF header from a set of header record lines
     *
     * @param headerStrings a list of strings that represent all the ## and # entries
     * @return a VCFHeader object
     */
    protected VCFHeader parseHeaderFromLines( final List<String> headerStrings, final VCFHeaderVersion version ) {
        this.version = version;

        Set<VCFHeaderLine> metaData = new LinkedHashSet<VCFHeaderLine>();
        Set<String> sampleNames = new LinkedHashSet<String>();
        int contigCounter = 0;
        // iterate over all the passed in strings
        for ( String str : headerStrings ) {
            if ( !str.startsWith(VCFHeader.METADATA_INDICATOR) ) {
                String[] strings = str.substring(1).split(VCFConstants.FIELD_SEPARATOR);
                if ( strings.length < VCFHeader.HEADER_FIELDS.values().length )
                    throw new TribbleException.InvalidHeader("there are not enough columns present in the header line: " + str);

                int arrayIndex = 0;
                for (VCFHeader.HEADER_FIELDS field : VCFHeader.HEADER_FIELDS.values()) {
                    try {
                        if (field != VCFHeader.HEADER_FIELDS.valueOf(strings[arrayIndex]))
                            throw new TribbleException.InvalidHeader("we were expecting column name '" + field + "' but we saw '" + strings[arrayIndex] + "'");
                    } catch (IllegalArgumentException e) {
                        throw new TribbleException.InvalidHeader("unknown column name '" + strings[arrayIndex] + "'; it does not match a legal column header name.");
                    }
                    arrayIndex++;
                }

                boolean sawFormatTag = false;
                if ( arrayIndex < strings.length ) {
                    if ( !strings[arrayIndex].equals("FORMAT") )
                        throw new TribbleException.InvalidHeader("we were expecting column name 'FORMAT' but we saw '" + strings[arrayIndex] + "'");
                    sawFormatTag = true;
                    arrayIndex++;
                }

                while ( arrayIndex < strings.length )
                    sampleNames.add(strings[arrayIndex++]);

                if ( sawFormatTag && sampleNames.size() == 0 )
                    throw new TribbleException.InvalidHeader("The FORMAT field was provided but there is no genotype/sample data");

            } else {
                if ( str.startsWith(VCFConstants.INFO_HEADER_START) ) {
                    final VCFInfoHeaderLine info = new VCFInfoHeaderLine(str.substring(7), version);
                    metaData.add(info);
                } else if ( str.startsWith(VCFConstants.FILTER_HEADER_START) ) {
                    final VCFFilterHeaderLine filter = new VCFFilterHeaderLine(str.substring(9), version);
                    metaData.add(filter);
                } else if ( str.startsWith(VCFConstants.FORMAT_HEADER_START) ) {
                    final VCFFormatHeaderLine format = new VCFFormatHeaderLine(str.substring(9), version);
                    metaData.add(format);
                } else if ( str.startsWith(VCFConstants.CONTIG_HEADER_START) ) {
                    final VCFContigHeaderLine contig = new VCFContigHeaderLine(str.substring(9), version, VCFConstants.CONTIG_HEADER_START.substring(2), contigCounter++);
                    metaData.add(contig);
                } else if ( str.startsWith(VCFConstants.ALT_HEADER_START) ) {
                    final VCFSimpleHeaderLine alt = new VCFSimpleHeaderLine(str.substring(6), version, VCFConstants.ALT_HEADER_START.substring(2), Arrays.asList("ID", "Description"));
                    metaData.add(alt);
                } else {
                    int equals = str.indexOf("=");
                    if ( equals != -1 )
                        metaData.add(new VCFHeaderLine(str.substring(2, equals), str.substring(equals+1)));
                }
            }
        }

        this.header = new VCFHeader(metaData, sampleNames);
        if ( doOnTheFlyModifications )
            this.header = VCFStandardHeaderLines.repairStandardHeaderLines(this.header);
        return this.header;
    }

    /**
     * the fast decode function
     * @param line the line of text for the record
     * @return a feature, (not guaranteed complete) that has the correct start and stop
     */
    public Feature decodeLoc(String line) {
        return decodeLine(line, false);
    }

    /**
     * decode the line into a feature (VariantContext)
     * @param line the line
     * @return a VariantContext
     */
    public VariantContext decode(String line) {
        return decodeLine(line, true);
    }

    private final VariantContext decodeLine(final String line, final boolean includeGenotypes) {
        // the same line reader is not used for parsing the header and parsing lines, if we see a #, we've seen a header line
        if (line.startsWith(VCFHeader.HEADER_INDICATOR)) return null;

        // our header cannot be null, we need the genotype sample names and counts
        if (header == null) throw new TribbleException("VCF Header cannot be null when decoding a record");

        if (parts == null)
            parts = new String[Math.min(header.getColumnCount(), NUM_STANDARD_FIELDS+1)];

        int nParts = ParsingUtils.split(line, parts, VCFConstants.FIELD_SEPARATOR_CHAR, true);

        // if we have don't have a header, or we have a header with no genotyping data check that we have eight columns.  Otherwise check that we have nine (normal colummns + genotyping data)
        if (( (header == null || !header.hasGenotypingData()) && nParts != NUM_STANDARD_FIELDS) ||
                (header != null && header.hasGenotypingData() && nParts != (NUM_STANDARD_FIELDS + 1)) )
            throw new TribbleException("Line " + lineNo + ": there aren't enough columns for line " + line + " (we expected " + (header == null ? NUM_STANDARD_FIELDS : NUM_STANDARD_FIELDS + 1) +
                    " tokens, and saw " + nParts + " )");

        return parseVCFLine(parts, includeGenotypes);
    }

    /**
     * parse out the VCF line
     *
     * @param parts the parts split up
     * @return a variant context object
     */
    private VariantContext parseVCFLine(final String[] parts, final boolean includeGenotypes) {
        VariantContextBuilder builder = new VariantContextBuilder();
        builder.source(getName());

        // increment the line count
        // TODO -- because of the way the engine utilizes Tribble, we can parse a line multiple times (especially when
        // TODO --   the first record is far along the contig) and the line counter can get out of sync
        lineNo++;

        // parse out the required fields
        final String chr = getCachedString(parts[0]);
        builder.chr(chr);
        int pos = -1;
        try {
            pos = Integer.valueOf(parts[1]);
        } catch (NumberFormatException e) {
            generateException(parts[1] + " is not a valid start position in the VCF format");
        }
        builder.start(pos);

        if ( parts[2].length() == 0 )
            generateException("The VCF specification requires a valid ID field");
        else if ( parts[2].equals(VCFConstants.EMPTY_ID_FIELD) )
            builder.noID();
        else
            builder.id(parts[2]);

        final String ref = getCachedString(parts[3].toUpperCase());
        final String alts = getCachedString(parts[4].toUpperCase());
        builder.log10PError(parseQual(parts[5]));

        final List<String> filters = parseFilters(getCachedString(parts[6]));
        if ( filters != null ) builder.filters(new HashSet<String>(filters));
        final Map<String, Object> attrs = parseInfo(parts[7]);
        builder.attributes(attrs);

        if ( attrs.containsKey(VCFConstants.END_KEY) ) {
            // update stop with the end key if provided
            try {
                builder.stop(Integer.valueOf(attrs.get(VCFConstants.END_KEY).toString()));
            } catch (Exception e) {
                generateException("the END value in the INFO field is not valid");
            }
        } else {
            builder.stop(pos + ref.length() - 1);
        }

        // get our alleles, filters, and setup an attribute map
        final List<Allele> alleles = parseAlleles(ref, alts, lineNo);
        builder.alleles(alleles);

        // do we have genotyping data
        if (parts.length > NUM_STANDARD_FIELDS && includeGenotypes) {
            final LazyGenotypesContext.LazyParser lazyParser = new LazyVCFGenotypesParser(alleles, chr, pos);
            final int nGenotypes = header.getNGenotypeSamples();
            LazyGenotypesContext lazy = new LazyGenotypesContext(lazyParser, parts[8], nGenotypes);

            // did we resort the sample names?  If so, we need to load the genotype data
            if ( !header.samplesWereAlreadySorted() )
                lazy.decode();

            builder.genotypesNoValidation(lazy);
        }

        VariantContext vc = null;
        try {
            vc = builder.make();
        } catch (Exception e) {
            generateException(e.getMessage());
        }

        return vc;
    }

    /**
     * get the name of this codec
     * @return our set name
     */
    public String getName() {
        return name;
    }

    /**
     * set the name of this codec
     * @param name new name
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * Return a cached copy of the supplied string.
     *
     * @param str string
     * @return interned string
     */
    protected String getCachedString(String str) {
        String internedString = stringCache.get(str);
        if ( internedString == null ) {
            internedString = new String(str);
            stringCache.put(internedString, internedString);
        }
        return internedString;
    }

    /**
     * parse out the info fields
     * @param infoField the fields
     * @return a mapping of keys to objects
     */
    private Map<String, Object> parseInfo(String infoField) {
        Map<String, Object> attributes = new HashMap<String, Object>();

        if ( infoField.length() == 0 )
            generateException("The VCF specification requires a valid info field");

        if ( !infoField.equals(VCFConstants.EMPTY_INFO_FIELD) ) {
            if ( infoField.indexOf("\t") != -1 || infoField.indexOf(" ") != -1 )
                generateException("The VCF specification does not allow for whitespace in the INFO field");

            int infoFieldSplitSize = ParsingUtils.split(infoField, infoFieldArray, VCFConstants.INFO_FIELD_SEPARATOR_CHAR, false);
            for (int i = 0; i < infoFieldSplitSize; i++) {
                String key;
                Object value;

                int eqI = infoFieldArray[i].indexOf("=");
                if ( eqI != -1 ) {
                    key = infoFieldArray[i].substring(0, eqI);
                    String valueString = infoFieldArray[i].substring(eqI+1);

                    // split on the INFO field separator
                    int infoValueSplitSize = ParsingUtils.split(valueString, infoValueArray, VCFConstants.INFO_FIELD_ARRAY_SEPARATOR_CHAR, false);
                    if ( infoValueSplitSize == 1 ) {
                        value = infoValueArray[0];
                        final VCFInfoHeaderLine headerLine = header.getInfoHeaderLine(key);
                        if ( headerLine != null && headerLine.getType() == VCFHeaderLineType.Flag && value.equals("0") ) {
                            // deal with the case where a flag field has =0, such as DB=0, by skipping the add
                            continue;
                        }
                    } else {
                        ArrayList<String> valueList = new ArrayList<String>(infoValueSplitSize);
                        for ( int j = 0; j < infoValueSplitSize; j++ )
                            valueList.add(infoValueArray[j]);
                        value = valueList;
                    }
                } else {
                    key = infoFieldArray[i];
                    final VCFInfoHeaderLine headerLine = header.getInfoHeaderLine(key);
                    if ( headerLine != null && headerLine.getType() != VCFHeaderLineType.Flag ) {
                        if ( GeneralUtils.DEBUG_MODE_ENABLED && ! warnedAboutNoEqualsForNonFlag ) {
                            System.err.println("Found info key " + key + " without a = value, but the header says the field is of type "
                                               + headerLine.getType() + " but this construct is only value for FLAG type fields");
                            warnedAboutNoEqualsForNonFlag = true;
                        }

                        value = VCFConstants.MISSING_VALUE_v4;
                    } else {
                        value = true;
                    }
                }

                // this line ensures that key/value pairs that look like key=; are parsed correctly as MISSING
                if ( "".equals(value) ) value = VCFConstants.MISSING_VALUE_v4;

                attributes.put(key, value);
            }
        }

        return attributes;
    }

    /**
     * create a an allele from an index and an array of alleles
     * @param index the index
     * @param alleles the alleles
     * @return an Allele
     */
    protected static Allele oneAllele(String index, List<Allele> alleles) {
        if ( index.equals(VCFConstants.EMPTY_ALLELE) )
            return Allele.NO_CALL;
        final int i;
        try {
            i = Integer.valueOf(index);
        } catch ( NumberFormatException e ) {
            throw new TribbleException.InternalCodecException("The following invalid GT allele index was encountered in the file: " + index);
        }
        if ( i >= alleles.size() )
            throw new TribbleException.InternalCodecException("The allele with index " + index + " is not defined in the REF/ALT columns in the record");
        return alleles.get(i);
    }


    /**
     * parse genotype alleles from the genotype string
     * @param GT         GT string
     * @param alleles    list of possible alleles
     * @param cache      cache of alleles for GT
     * @return the allele list for the GT string
     */
    protected static List<Allele> parseGenotypeAlleles(String GT, List<Allele> alleles, Map<String, List<Allele>> cache) {
        // cache results [since they are immutable] and return a single object for each genotype
        List<Allele> GTAlleles = cache.get(GT);

        if ( GTAlleles == null ) {
            StringTokenizer st = new StringTokenizer(GT, VCFConstants.PHASING_TOKENS);
            GTAlleles = new ArrayList<Allele>(st.countTokens());
            while ( st.hasMoreTokens() ) {
                String genotype = st.nextToken();
                GTAlleles.add(oneAllele(genotype, alleles));
            }
            cache.put(GT, GTAlleles);
        }

        return GTAlleles;
    }

    /**
     * parse out the qual value
     * @param qualString the quality string
     * @return return a double
     */
    protected static Double parseQual(String qualString) {
        // if we're the VCF 4 missing char, return immediately
        if ( qualString.equals(VCFConstants.MISSING_VALUE_v4))
            return VariantContext.NO_LOG10_PERROR;

        Double val = Double.valueOf(qualString);

        // check to see if they encoded the missing qual score in VCF 3 style, with either the -1 or -1.0.  check for val < 0 to save some CPU cycles
        if ((val < 0) && (Math.abs(val - VCFConstants.MISSING_QUALITY_v3_DOUBLE) < VCFConstants.VCF_ENCODING_EPSILON))
            return VariantContext.NO_LOG10_PERROR;

        // scale and return the value
        return val / -10.0;
    }

    /**
     * parse out the alleles
     * @param ref the reference base
     * @param alts a string of alternates to break into alleles
     * @param lineNo  the line number for this record
     * @return a list of alleles, and a pair of the shortest and longest sequence
     */
    protected static List<Allele> parseAlleles(String ref, String alts, int lineNo) {
        List<Allele> alleles = new ArrayList<Allele>(2); // we are almost always biallelic
        // ref
        checkAllele(ref, true, lineNo);
        Allele refAllele = Allele.create(ref, true);
        alleles.add(refAllele);

        if ( alts.indexOf(",") == -1 ) // only 1 alternatives, don't call string split
            parseSingleAltAllele(alleles, alts, lineNo);
        else
            for ( String alt : alts.split(",") )
                parseSingleAltAllele(alleles, alt, lineNo);

        return alleles;
    }

    /**
     * check to make sure the allele is an acceptable allele
     * @param allele the allele to check
     * @param isRef are we the reference allele?
     * @param lineNo  the line number for this record
     */
    private static void checkAllele(String allele, boolean isRef, int lineNo) {
        if ( allele == null || allele.length() == 0 )
            generateException("Empty alleles are not permitted in VCF records", lineNo);

        if ( GeneralUtils.DEBUG_MODE_ENABLED && MAX_ALLELE_SIZE_BEFORE_WARNING != -1 && allele.length() > MAX_ALLELE_SIZE_BEFORE_WARNING ) {
            System.err.println(String.format("Allele detected with length %d exceeding max size %d at approximately line %d, likely resulting in degraded VCF processing performance", allele.length(), MAX_ALLELE_SIZE_BEFORE_WARNING, lineNo));
        }

        if ( isSymbolicAllele(allele) ) {
            if ( isRef ) {
                generateException("Symbolic alleles not allowed as reference allele: " + allele, lineNo);
            }
        } else {
            // check for VCF3 insertions or deletions
            if ( (allele.charAt(0) == VCFConstants.DELETION_ALLELE_v3) || (allele.charAt(0) == VCFConstants.INSERTION_ALLELE_v3) )
                generateException("Insertions/Deletions are not supported when reading 3.x VCF's. Please" +
                        " convert your file to VCF4 using VCFTools, available at http://vcftools.sourceforge.net/index.html", lineNo);

            if (!Allele.acceptableAlleleBases(allele))
                generateException("Unparsable vcf record with allele " + allele, lineNo);

            if ( isRef && allele.equals(VCFConstants.EMPTY_ALLELE) )
                generateException("The reference allele cannot be missing", lineNo);
        }
    }

    /**
     * return true if this is a symbolic allele (e.g. <SOMETAG>) or
     * structural variation breakend (with [ or ]), otherwise false
     * @param allele the allele to check
     * @return true if the allele is a symbolic allele, otherwise false
     */
    private static boolean isSymbolicAllele(String allele) {
        return (allele != null && allele.length() > 2 &&
                ((allele.startsWith("<") && allele.endsWith(">")) ||
                        (allele.contains("[") || allele.contains("]"))));
    }

    /**
     * parse a single allele, given the allele list
     * @param alleles the alleles available
     * @param alt the allele to parse
     * @param lineNo  the line number for this record
     */
    private static void parseSingleAltAllele(List<Allele> alleles, String alt, int lineNo) {
        checkAllele(alt, false, lineNo);

        Allele allele = Allele.create(alt, false);
        if ( ! allele.isNoCall() )
            alleles.add(allele);
    }

    public final static boolean canDecodeFile(final String potentialInput, final String MAGIC_HEADER_LINE) {
        try {
            return isVCFStream(new FileInputStream(potentialInput), MAGIC_HEADER_LINE) ||
                    isVCFStream(new GZIPInputStream(new FileInputStream(potentialInput)), MAGIC_HEADER_LINE) ||
                    isVCFStream(new BlockCompressedInputStream(new FileInputStream(potentialInput)), MAGIC_HEADER_LINE);
        } catch ( FileNotFoundException e ) {
            return false;
        } catch ( IOException e ) {
            return false;
        }
    }

    private final static boolean isVCFStream(final InputStream stream, final String MAGIC_HEADER_LINE) {
        try {
            byte[] buff = new byte[MAGIC_HEADER_LINE.length()];
            int nread = stream.read(buff, 0, MAGIC_HEADER_LINE.length());
            boolean eq = Arrays.equals(buff, MAGIC_HEADER_LINE.getBytes());
            return eq;
//            String firstLine = new String(buff);
//            return firstLine.startsWith(MAGIC_HEADER_LINE);
        } catch ( IOException e ) {
            return false;
        } catch ( RuntimeException e ) {
            return false;
        } finally {
            try { stream.close(); } catch ( IOException e ) {}
        }
    }


    /**
     * create a genotype map
     *
     * @param str the string
     * @param alleles the list of alleles
     * @return a mapping of sample name to genotype object
     */
    public LazyGenotypesContext.LazyData createGenotypeMap(final String str,
                                                              final List<Allele> alleles,
                                                              final String chr,
                                                              final int pos) {
        if (genotypeParts == null)
            genotypeParts = new String[header.getColumnCount() - NUM_STANDARD_FIELDS];

        int nParts = ParsingUtils.split(str, genotypeParts, VCFConstants.FIELD_SEPARATOR_CHAR);
        if ( nParts != genotypeParts.length )
            generateException("there are " + (nParts-1) + " genotypes while the header requires that " + (genotypeParts.length-1) + " genotypes be present for all records at " + chr + ":" + pos, lineNo);

        ArrayList<Genotype> genotypes = new ArrayList<Genotype>(nParts);

        // get the format keys
        int nGTKeys = ParsingUtils.split(genotypeParts[0], genotypeKeyArray, VCFConstants.GENOTYPE_FIELD_SEPARATOR_CHAR);

        // cycle through the sample names
        Iterator<String> sampleNameIterator = header.getGenotypeSamples().iterator();

        // clear out our allele mapping
        alleleMap.clear();

        // cycle through the genotype strings
        for (int genotypeOffset = 1; genotypeOffset < nParts; genotypeOffset++) {
            int GTValueSplitSize = ParsingUtils.split(genotypeParts[genotypeOffset], GTValueArray, VCFConstants.GENOTYPE_FIELD_SEPARATOR_CHAR);

            final String sampleName = sampleNameIterator.next();
            final GenotypeBuilder gb = new GenotypeBuilder(sampleName);

            // check to see if the value list is longer than the key list, which is a problem
            if (nGTKeys < GTValueSplitSize)
                generateException("There are too many keys for the sample " + sampleName + ", keys = " + parts[8] + ", values = " + parts[genotypeOffset]);

            int genotypeAlleleLocation = -1;
            if (nGTKeys >= 1) {
                gb.maxAttributes(nGTKeys - 1);

                for (int i = 0; i < nGTKeys; i++) {
                    final String gtKey = genotypeKeyArray[i];
                    boolean missing = i >= GTValueSplitSize;

                    // todo -- all of these on the fly parsing of the missing value should be static constants
                    if (gtKey.equals(VCFConstants.GENOTYPE_KEY)) {
                        genotypeAlleleLocation = i;
                    } else if ( missing ) {
                        // if its truly missing (there no provided value) skip adding it to the attributes
                    } else if (gtKey.equals(VCFConstants.GENOTYPE_FILTER_KEY)) {
                        final List<String> filters = parseFilters(getCachedString(GTValueArray[i]));
                        if ( filters != null ) gb.filters(filters);
                    } else if ( GTValueArray[i].equals(VCFConstants.MISSING_VALUE_v4) ) {
                        // don't add missing values to the map
                    } else {
                        if (gtKey.equals(VCFConstants.GENOTYPE_QUALITY_KEY)) {
                            if ( GTValueArray[i].equals(VCFConstants.MISSING_GENOTYPE_QUALITY_v3) )
                                gb.noGQ();
                            else
                                gb.GQ((int)Math.round(Double.valueOf(GTValueArray[i])));
                        } else if (gtKey.equals(VCFConstants.GENOTYPE_ALLELE_DEPTHS)) {
                            gb.AD(decodeInts(GTValueArray[i]));
                        } else if (gtKey.equals(VCFConstants.GENOTYPE_PL_KEY)) {
                            gb.PL(decodeInts(GTValueArray[i]));
                        } else if (gtKey.equals(VCFConstants.GENOTYPE_LIKELIHOODS_KEY)) {
                            gb.PL(GenotypeLikelihoods.fromGLField(GTValueArray[i]).getAsPLs());
                        } else if (gtKey.equals(VCFConstants.DEPTH_KEY)) {
                            gb.DP(Integer.valueOf(GTValueArray[i]));
                        } else {
                            gb.attribute(gtKey, GTValueArray[i]);
                        }
                    }
                }
            }

            // check to make sure we found a genotype field if our version is less than 4.1 file
            if ( version != VCFHeaderVersion.VCF4_1 && genotypeAlleleLocation == -1 )
                generateException("Unable to find the GT field for the record; the GT field is required in VCF4.0");
            if ( genotypeAlleleLocation > 0 )
                generateException("Saw GT field at position " + genotypeAlleleLocation + ", but it must be at the first position for genotypes when present");

            final List<Allele> GTalleles = (genotypeAlleleLocation == -1 ? new ArrayList<Allele>(0) : parseGenotypeAlleles(GTValueArray[genotypeAlleleLocation], alleles, alleleMap));
            gb.alleles(GTalleles);
            gb.phased(genotypeAlleleLocation != -1 && GTValueArray[genotypeAlleleLocation].indexOf(VCFConstants.PHASED) != -1);

            // add it to the list
            try {
                genotypes.add(gb.make());
            } catch (TribbleException e) {
                throw new TribbleException.InternalCodecException(e.getMessage() + ", at position " + chr+":"+pos);
            }
        }

        return new LazyGenotypesContext.LazyData(genotypes, header.getSampleNamesInOrder(), header.getSampleNameToOffset());
    }


    private final static String[] INT_DECODE_ARRAY = new String[10000];
    private final static int[] decodeInts(final String string) {
        final int nValues = ParsingUtils.split(string, INT_DECODE_ARRAY, ',');
        final int[] values = new int[nValues];
        for ( int i = 0; i < nValues; i++ )
            values[i] = Integer.valueOf(INT_DECODE_ARRAY[i]);
        return values;
    }

    /**
     * Forces all VCFCodecs to not perform any on the fly modifications to the VCF header
     * of VCF records.  Useful primarily for raw comparisons such as when comparing
     * raw VCF records
     */
    public final void disableOnTheFlyModifications() {
        doOnTheFlyModifications = false;
    }


    protected void generateException(String message) {
        throw new TribbleException(String.format("The provided VCF file is malformed at approximately line number %d: %s", lineNo, message));
    }

    protected static void generateException(String message, int lineNo) {
        throw new TribbleException(String.format("The provided VCF file is malformed at approximately line number %d: %s", lineNo, message));
    }
}
