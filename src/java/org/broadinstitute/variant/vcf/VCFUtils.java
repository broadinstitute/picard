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

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import org.broadinstitute.variant.utils.GeneralUtils;

import java.io.File;
import java.util.*;

public class VCFUtils {

    public static Set<VCFHeaderLine> smartMergeHeaders(Collection<VCFHeader> headers, boolean emitWarnings) throws IllegalStateException {
        HashMap<String, VCFHeaderLine> map = new HashMap<String, VCFHeaderLine>(); // from KEY.NAME -> line
        HeaderConflictWarner conflictWarner = new HeaderConflictWarner(emitWarnings);

        // todo -- needs to remove all version headers from sources and add its own VCF version line
        for ( VCFHeader source : headers ) {
            //System.out.printf("Merging in header %s%n", source);
            for ( VCFHeaderLine line : source.getMetaDataInSortedOrder()) {

                String key = line.getKey();
                if ( line instanceof VCFIDHeaderLine )
                    key = key + "-" + ((VCFIDHeaderLine)line).getID();

                if ( map.containsKey(key) ) {
                    VCFHeaderLine other = map.get(key);
                    if ( line.equals(other) ) {
                        // continue;
                    } else if ( ! line.getClass().equals(other.getClass()) ) {
                        throw new IllegalStateException("Incompatible header types: " + line + " " + other );
                    } else if ( line instanceof VCFFilterHeaderLine ) {
                        String lineName = ((VCFFilterHeaderLine) line).getID();
                        String otherName = ((VCFFilterHeaderLine) other).getID();
                        if ( ! lineName.equals(otherName) )
                            throw new IllegalStateException("Incompatible header types: " + line + " " + other );
                    } else if ( line instanceof VCFCompoundHeaderLine ) {
                        VCFCompoundHeaderLine compLine = (VCFCompoundHeaderLine)line;
                        VCFCompoundHeaderLine compOther = (VCFCompoundHeaderLine)other;

                        // if the names are the same, but the values are different, we need to quit
                        if (! (compLine).equalsExcludingDescription(compOther) ) {
                            if ( compLine.getType().equals(compOther.getType()) ) {
                                // The Number entry is an Integer that describes the number of values that can be
                                // included with the INFO field. For example, if the INFO field contains a single
                                // number, then this value should be 1. However, if the INFO field describes a pair
                                // of numbers, then this value should be 2 and so on. If the number of possible
                                // values varies, is unknown, or is unbounded, then this value should be '.'.
                                conflictWarner.warn(line, "Promoting header field Number to . due to number differences in header lines: " + line + " " + other);
                                compOther.setNumberToUnbounded();
                            } else if ( compLine.getType() == VCFHeaderLineType.Integer && compOther.getType() == VCFHeaderLineType.Float ) {
                                // promote key to Float
                                conflictWarner.warn(line, "Promoting Integer to Float in header: " + compOther);
                                map.put(key, compOther);
                            } else if ( compLine.getType() == VCFHeaderLineType.Float && compOther.getType() == VCFHeaderLineType.Integer ) {
                                // promote key to Float
                                conflictWarner.warn(line, "Promoting Integer to Float in header: " + compOther);
                            } else {
                                throw new IllegalStateException("Incompatible header types, collision between these two types: " + line + " " + other );
                            }
                        }
                        if ( ! compLine.getDescription().equals(compOther.getDescription()) )
                            conflictWarner.warn(line, "Allowing unequal description fields through: keeping " + compOther + " excluding " + compLine);
                    } else {
                        // we are not equal, but we're not anything special either
                        conflictWarner.warn(line, "Ignoring header line already in map: this header line = " + line + " already present header = " + other);
                    }
                } else {
                    map.put(key, line);
                    //System.out.printf("Adding header line %s%n", line);
                }
            }
        }

        return new HashSet<VCFHeaderLine>(map.values());
    }

    /**
     * Add / replace the contig header lines in the VCFHeader with the in the reference file and master reference dictionary
     *
     * @param oldHeader the header to update
     * @param referenceFile the file path to the reference sequence used to generate this vcf
     * @param refDict the SAM formatted reference sequence dictionary
     */
    public static VCFHeader withUpdatedContigs(final VCFHeader oldHeader, final File referenceFile, final SAMSequenceDictionary refDict) {
        return new VCFHeader(withUpdatedContigsAsLines(oldHeader.getMetaDataInInputOrder(), referenceFile, refDict), oldHeader.getGenotypeSamples());
    }

    public static Set<VCFHeaderLine> withUpdatedContigsAsLines(final Set<VCFHeaderLine> oldLines, final File referenceFile, final SAMSequenceDictionary refDict) {
        return withUpdatedContigsAsLines(oldLines, referenceFile, refDict, false);
    }

    public static Set<VCFHeaderLine> withUpdatedContigsAsLines(final Set<VCFHeaderLine> oldLines, final File referenceFile, final SAMSequenceDictionary refDict, boolean referenceNameOnly) {
        final Set<VCFHeaderLine> lines = new LinkedHashSet<VCFHeaderLine>(oldLines.size());

        for ( final VCFHeaderLine line : oldLines ) {
            if ( line instanceof VCFContigHeaderLine )
                continue; // skip old contig lines
            if ( line.getKey().equals(VCFHeader.REFERENCE_KEY) )
                continue; // skip the old reference key
            lines.add(line);
        }

        for ( final VCFHeaderLine contigLine : makeContigHeaderLines(refDict, referenceFile) )
            lines.add(contigLine);

        String referenceValue;
        if (referenceFile != null) {
            if (referenceNameOnly) {
                int extensionStart = referenceFile.getName().lastIndexOf(".");
                referenceValue = extensionStart == -1 ? referenceFile.getName() : referenceFile.getName().substring(0, extensionStart);
            }
            else {
                referenceValue = "file://" + referenceFile.getAbsolutePath();
            }
            lines.add(new VCFHeaderLine(VCFHeader.REFERENCE_KEY, referenceValue));
        }
        return lines;
    }

    /**
     * Create VCFHeaderLines for each refDict entry, and optionally the assembly if referenceFile != null
     * @param refDict reference dictionary
     * @param referenceFile for assembly name.  May be null
     * @return list of vcf contig header lines
     */
    public static List<VCFContigHeaderLine> makeContigHeaderLines(final SAMSequenceDictionary refDict,
                                                                  final File referenceFile) {
        final List<VCFContigHeaderLine> lines = new ArrayList<VCFContigHeaderLine>();
        final String assembly = referenceFile != null ? getReferenceAssembly(referenceFile.getName()) : null;
        for ( SAMSequenceRecord contig : refDict.getSequences() )
            lines.add(makeContigHeaderLine(contig, assembly));
        return lines;
    }

    private static VCFContigHeaderLine makeContigHeaderLine(final SAMSequenceRecord contig, final String assembly) {
        final Map<String, String> map = new LinkedHashMap<String, String>(3);
        map.put("ID", contig.getSequenceName());
        map.put("length", String.valueOf(contig.getSequenceLength()));
        if ( assembly != null ) map.put("assembly", assembly);
        return new VCFContigHeaderLine(map, contig.getSequenceIndex());
    }

    private static String getReferenceAssembly(final String refPath) {
        // This doesn't need to be perfect as it's not a required VCF header line, but we might as well give it a shot
        String assembly = null;
        if (refPath.contains("b37") || refPath.contains("v37"))
            assembly = "b37";
        else if (refPath.contains("b36"))
            assembly = "b36";
        else if (refPath.contains("hg18"))
            assembly = "hg18";
        else if (refPath.contains("hg19"))
            assembly = "hg19";
        return assembly;
    }

    /** Only displays a warning if warnings are enabled and an identical warning hasn't been already issued */
    private static final class HeaderConflictWarner {
        boolean emitWarnings;
        Set<String> alreadyIssued = new HashSet<String>();

        private HeaderConflictWarner( final boolean emitWarnings ) {
            this.emitWarnings = emitWarnings;
        }

        public void warn(final VCFHeaderLine line, final String msg) {
            if ( GeneralUtils.DEBUG_MODE_ENABLED && emitWarnings && ! alreadyIssued.contains(line.getKey()) ) {
                alreadyIssued.add(line.getKey());
                System.err.println(msg);
            }
        }
    }
}
