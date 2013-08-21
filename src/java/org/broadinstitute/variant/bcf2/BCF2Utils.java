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

package org.broadinstitute.variant.bcf2;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broad.tribble.TribbleException;
import org.broadinstitute.variant.vcf.*;

import java.io.*;
import java.util.*;

/**
 * Common utilities for working with BCF2 files
 *
 * Includes convenience methods for encoding, decoding BCF2 type descriptors (size + type)
 *
 * @author depristo
 * @since 5/12
 */
public final class BCF2Utils {
    public static final int MAX_ALLELES_IN_GENOTYPES = 127;

    public static final int OVERFLOW_ELEMENT_MARKER = 15;
    public static final int MAX_INLINE_ELEMENTS = 14;

    public final static BCF2Type[] INTEGER_TYPES_BY_SIZE = new BCF2Type[]{BCF2Type.INT8, BCF2Type.INT16, BCF2Type.INT32};
    public final static BCF2Type[] ID_TO_ENUM;

    static {
        int maxID = -1;
        for ( BCF2Type v : BCF2Type.values() ) maxID = Math.max(v.getID(), maxID);
        ID_TO_ENUM = new BCF2Type[maxID+1];
        for ( BCF2Type v : BCF2Type.values() ) ID_TO_ENUM[v.getID()] = v;
    }

    private BCF2Utils() {}

    /**
     * Create a strings dictionary from the VCF header
     *
     * The dictionary is an ordered list of common VCF identifers (FILTER, INFO, and FORMAT)
     * fields.
     *
     * Note that its critical that the list be dedupped and sorted in a consistent manner each time,
     * as the BCF2 offsets are encoded relative to this dictionary, and if it isn't determined exactly
     * the same way as in the header each time it's very bad
     *
     * @param header the VCFHeader from which to build the dictionary
     * @return a non-null dictionary of elements, may be empty
     */
    @Requires("header != null")
    @Ensures({"result != null", "new HashSet(result).size() == result.size()"})
    public static ArrayList<String> makeDictionary(final VCFHeader header) {
        final Set<String> seen = new HashSet<String>();
        final ArrayList<String> dict = new ArrayList<String>();

        // special case the special PASS field which doesn't show up in the FILTER field definitions
        seen.add(VCFConstants.PASSES_FILTERS_v4);
        dict.add(VCFConstants.PASSES_FILTERS_v4);

        // set up the strings dictionary
        for ( VCFHeaderLine line : header.getMetaDataInInputOrder() ) {
            if ( line.shouldBeAddedToDictionary() ) {
                final VCFIDHeaderLine idLine = (VCFIDHeaderLine)line;
                if ( ! seen.contains(idLine.getID())) {
                    dict.add(idLine.getID());
                    seen.add(idLine.getID());
                }
            }
        }

        return dict;
    }

    @Requires({"nElements >= 0", "nElements <= OVERFLOW_ELEMENT_MARKER", "type != null"})
    public static byte encodeTypeDescriptor(final int nElements, final BCF2Type type ) {
        return (byte)((0x0F & nElements) << 4 | (type.getID() & 0x0F));
    }

    @Ensures("result >= 0")
    public static int decodeSize(final byte typeDescriptor) {
        return (0xF0 & typeDescriptor) >> 4;
    }

    @Ensures("result >= 0")
    public static int decodeTypeID(final byte typeDescriptor) {
        return typeDescriptor & 0x0F;
    }

    @Ensures("result != null")
    public static BCF2Type decodeType(final byte typeDescriptor) {
        return ID_TO_ENUM[decodeTypeID(typeDescriptor)];
    }

    public static boolean sizeIsOverflow(final byte typeDescriptor) {
        return decodeSize(typeDescriptor) == OVERFLOW_ELEMENT_MARKER;
    }

    public static byte readByte(final InputStream stream) throws IOException {
        return (byte)(stream.read() & 0xFF);
    }

    /**
     * Collapse multiple strings into a comma separated list
     *
     * ["s1", "s2", "s3"] => ",s1,s2,s3"
     *
     * @param strings size > 1 list of strings
     * @return
     */
    @Requires({"strings != null"})
    @Ensures("result != null")
    public static String collapseStringList(final List<String> strings) {
        if ( strings.isEmpty() ) return "";
        else if ( strings.size() == 1 ) return strings.get(0);
        else {
            final StringBuilder b = new StringBuilder();
            for ( final String s : strings ) {
                if ( s != null ) {
                    assert s.indexOf(",") == -1; // no commas in individual strings
                    b.append(",").append(s);
                }
            }
            return b.toString();
        }
    }

    /**
     * Inverse operation of collapseStringList.
     *
     * ",s1,s2,s3" => ["s1", "s2", "s3"]
     *
     *
     * @param collapsed
     * @return
     */
    @Requires({"collapsed != null", "isCollapsedString(collapsed)"})
    @Ensures("result != null")
    public static List<String> explodeStringList(final String collapsed) {
        assert isCollapsedString(collapsed);
        final String[] exploded = collapsed.substring(1).split(",");
        return Arrays.asList(exploded);
    }

    @Requires("s != null")
    public static boolean isCollapsedString(final String s) {
        return s.length() > 0 && s.charAt(0) == ',';
    }

    /**
     * Returns a good name for a shadow BCF file for vcfFile.
     *
     * foo.vcf => foo.bcf
     * foo.xxx => foo.xxx.bcf
     *
     * If the resulting BCF file cannot be written, return null.  Happens
     * when vcfFile = /dev/null for example
     *
     * @param vcfFile
     * @return the BCF
     */
    @Requires("vcfFile != null")
    public static final File shadowBCF(final File vcfFile) {
        final String path = vcfFile.getAbsolutePath();
        if ( path.contains(".vcf") )
            return new File(path.replace(".vcf", ".bcf"));
        else {
            final File bcf = new File( path + ".bcf" );
            if ( bcf.canRead() )
                return bcf;
            else {
                try {
                    // this is the only way to robustly decide if we could actually write to BCF
                    final FileOutputStream o = new FileOutputStream(bcf);
                    o.close();
                    bcf.delete();
                    return bcf;
                } catch ( FileNotFoundException e ) {
                    return null;
                } catch ( IOException e ) {
                    return null;
                }
            }
        }
    }

    @Ensures("result.isIntegerType()")
    public static BCF2Type determineIntegerType(final int value) {
        for ( final BCF2Type potentialType : INTEGER_TYPES_BY_SIZE) {
            if ( potentialType.withinRange(value) )
                return potentialType;
        }

        throw new TribbleException("Integer cannot be encoded in allowable range of even INT32: " + value);
    }

    @Ensures("result.isIntegerType()")
    public static BCF2Type determineIntegerType(final int[] values) {
        // find the min and max values in the array
        int max = 0, min = 0;
        for ( final int v : values ) {
            if ( v > max ) max = v;
            if ( v < min ) min = v;
        }

        final BCF2Type maxType = determineIntegerType(max);
        final BCF2Type minType = determineIntegerType(min);

        // INT8 < INT16 < INT32 so this returns the larger of the two
        return maxType.compareTo(minType) >= 0 ? maxType : minType;
    }

    /**
     * Returns the maximum BCF2 integer size of t1 and t2
     *
     * For example, if t1 == INT8 and t2 == INT16 returns INT16
     *
     * @param t1
     * @param t2
     * @return
     */
    @Requires({"t1.isIntegerType()","t2.isIntegerType()"})
    @Ensures("result.isIntegerType()")
    public static BCF2Type maxIntegerType(final BCF2Type t1, final BCF2Type t2) {
        switch ( t1 ) {
            case INT8: return t2;
            case INT16: return t2 == BCF2Type.INT32 ? t2 : t1;
            case INT32: return t1;
            default: throw new TribbleException("BUG: unexpected BCF2Type " + t1);
        }
    }

    @Ensures("result.isIntegerType()")
    public static BCF2Type determineIntegerType(final List<Integer> values) {
        BCF2Type maxType = BCF2Type.INT8;
        for ( final int value : values ) {
            final BCF2Type type1 = determineIntegerType(value);
            switch ( type1 ) {
                case INT8: break;
                case INT16: maxType = BCF2Type.INT16; break;
                case INT32: return BCF2Type.INT32; // fast path for largest possible value
                default: throw new TribbleException("Unexpected integer type " + type1 );
            }
        }
        return maxType;
    }

    /**
     * Helper function that takes an object and returns a list representation
     * of it:
     *
     * o == null => []
     * o is a list => o
     * else => [o]
     *
     * @param o
     * @return
     */
    public static List<Object> toList(final Object o) {
        if ( o == null ) return Collections.emptyList();
        else if ( o instanceof List ) return (List<Object>)o;
        else return Collections.singletonList(o);
    }

    /**
     * Are the elements and their order in the output and input headers consistent so that
     * we can write out the raw genotypes block without decoding and recoding it?
     *
     * If the order of INFO, FILTER, or contrig elements in the output header is different than
     * in the input header we must decode the blocks using the input header and then recode them
     * based on the new output order.
     *
     * If they are consistent, we can simply pass through the raw genotypes block bytes, which is
     * a *huge* performance win for large blocks.
     *
     * Many common operations on BCF2 files (merging them for -nt, selecting a subset of records, etc)
     * don't modify the ordering of the header fields and so can safely pass through the genotypes
     * undecoded.  Some operations -- those at add filters or info fields -- can change the ordering
     * of the header fields and so produce invalid BCF2 files if the genotypes aren't decoded
     */
    public static boolean headerLinesAreOrderedConsistently(final VCFHeader outputHeader, final VCFHeader genotypesBlockHeader) {
        // first, we have to have the same samples in the same order
        if ( ! nullAsEmpty(outputHeader.getSampleNamesInOrder()).equals(nullAsEmpty(genotypesBlockHeader.getSampleNamesInOrder())) )
            return false;

        final Iterator<? extends VCFIDHeaderLine> outputLinesIt = outputHeader.getIDHeaderLines().iterator();
        final Iterator<? extends VCFIDHeaderLine> inputLinesIt = genotypesBlockHeader.getIDHeaderLines().iterator();

        while ( inputLinesIt.hasNext() ) {
            if ( ! outputLinesIt.hasNext() ) // missing lines in output
                return false;

            final VCFIDHeaderLine outputLine = outputLinesIt.next();
            final VCFIDHeaderLine inputLine = inputLinesIt.next();

            if ( ! inputLine.getClass().equals(outputLine.getClass()) || ! inputLine.getID().equals(outputLine.getID()) )
                return false;
        }

        return true;
    }

    private static <T> List<T> nullAsEmpty(List<T> l) {
        if ( l == null )
            return Collections.emptyList();
        else
            return l;
    }
}
