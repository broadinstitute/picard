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
import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.readers.LineReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/**
 * A feature codec for the VCF3 specification, to read older VCF files.  VCF3 has been
 * depreciated in favor of VCF4 (See VCF codec for the latest information)
 *
 * <p>
 * Reads historical VCF3 encoded files (1000 Genomes Pilot results, for example)
 * </p>
 *
 * <p>
 * See also: @see <a href="http://vcftools.sourceforge.net/specs.html">VCF specification</a><br>
 * See also: @see <a href="http://www.ncbi.nlm.nih.gov/pubmed/21653522">VCF spec. publication</a>
 * </p>
 *
 * @author Mark DePristo
 * @since 2010
 */
public class VCF3Codec extends AbstractVCFCodec {
    public final static String VCF3_MAGIC_HEADER = "##fileformat=VCFv3";

    /**
     * @param reader the line reader to take header lines from
     * @return the number of header lines
     */
    public Object readActualHeader(final LineIterator reader) {
        final List<String> headerStrings = new ArrayList<String>();

        VCFHeaderVersion version = null;
        boolean foundHeaderVersion = false;
        while (reader.hasNext()) {
            lineNo++;
            final String line = reader.peek();
            if (line.startsWith(VCFHeader.METADATA_INDICATOR)) {
                final String[] lineFields = line.substring(2).split("=");
                if (lineFields.length == 2 && VCFHeaderVersion.isFormatString(lineFields[0]) ) {
                    if ( !VCFHeaderVersion.isVersionString(lineFields[1]) )
                        throw new TribbleException.InvalidHeader(lineFields[1] + " is not a supported version");
                    foundHeaderVersion = true;
                    version = VCFHeaderVersion.toHeaderVersion(lineFields[1]);
                    if ( version != VCFHeaderVersion.VCF3_3 && version != VCFHeaderVersion.VCF3_2 )
                        throw new TribbleException.InvalidHeader("This codec is strictly for VCFv3 and does not support " + lineFields[1]);
                }
                headerStrings.add(reader.next());
            }
            else if (line.startsWith(VCFHeader.HEADER_INDICATOR)) {
                if (!foundHeaderVersion) {
                    throw new TribbleException.InvalidHeader("We never saw a header line specifying VCF version");
                }
                headerStrings.add(reader.next());
                return super.parseHeaderFromLines(headerStrings, version);
            }
            else {
                throw new TribbleException.InvalidHeader("We never saw the required CHROM header line (starting with one #) for the input VCF file");
            }

        }
        throw new TribbleException.InvalidHeader("We never saw the required CHROM header line (starting with one #) for the input VCF file");
    }


    /**
     * parse the filter string, first checking to see if we already have parsed it in a previous attempt
     * @param filterString the string to parse
     * @return a set of the filters applied
     */
    protected List<String> parseFilters(String filterString) {

        // null for unfiltered
        if ( filterString.equals(VCFConstants.UNFILTERED) )
            return null;

        // empty set for passes filters
        List<String> fFields = new ArrayList<String>();

        if ( filterString.equals(VCFConstants.PASSES_FILTERS_v3) )
            return new ArrayList<String>(fFields);

        if ( filterString.length() == 0 )
            generateException("The VCF specification requires a valid filter status");

        // do we have the filter string cached?
        if ( filterHash.containsKey(filterString) )
            return new ArrayList<String>(filterHash.get(filterString));

        // otherwise we have to parse and cache the value
        if ( filterString.indexOf(VCFConstants.FILTER_CODE_SEPARATOR) == -1 )
            fFields.add(filterString);
        else
            fFields.addAll(Arrays.asList(filterString.split(VCFConstants.FILTER_CODE_SEPARATOR)));

        filterHash.put(filterString, fFields);

        return fFields;
    }

    @Override
    public boolean canDecode(final String potentialInput) {
        return canDecodeFile(potentialInput, VCF3_MAGIC_HEADER);
    }
}
