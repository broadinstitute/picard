/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

package picard.analysis;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/** Utilities to calculate GC Bias
 * Created by kbergin on 9/23/15.
 */
public class GcBiasUtils {

    /////////////////////////////////////////////////////////////////////////////
    // Calculates GC as a number from 0 to 100 in the specified window.
    // If the window includes more than five no-calls then -1 is returned.
    /////////////////////////////////////////////////////////////////////////////
    public static int calculateGc(final byte[] bases, final int startIndex, final int endIndex, final CalculateGcState state) {
        if (state.init) {
            state.init = false;
            state.gcCount = 0;
            state.nCount = 0;
            for (int i = startIndex; i < endIndex; ++i) {
                final byte base = bases[i];
                if (SequenceUtil.basesEqual(base, (byte)'G') || SequenceUtil.basesEqual(base, (byte)'C')) ++state.gcCount;
                else if (SequenceUtil.basesEqual(base, (byte)'N')) ++state.nCount;
            }
        } else {
            final byte newBase = bases[endIndex - 1];
            if (SequenceUtil.basesEqual(newBase, (byte)'G') || SequenceUtil.basesEqual(newBase, (byte)'C')) ++state.gcCount;
            else if (newBase == 'N') ++state.nCount;

            if (SequenceUtil.basesEqual(state.priorBase, (byte)'G') || SequenceUtil.basesEqual(state.priorBase, (byte)'C')) --state.gcCount;
            else if (SequenceUtil.basesEqual(state.priorBase, (byte)'N')) --state.nCount;
        }
        state.priorBase = bases[startIndex];
        if (state.nCount > 4) return -1;
        else return (state.gcCount * 100) / (endIndex - startIndex);
    }

    /////////////////////////////////////////////////////////////////////////////
    // Calculate number of 100bp windows in the refBases passed in that fall into
    // each gc content bin (0-100% gc)
    /////////////////////////////////////////////////////////////////////////////
    public static int[] calculateRefWindowsByGc(final int windows, final File referenceSequence, final int windowSize) {
        final ReferenceSequenceFile refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceSequence);
        ReferenceSequence ref;

        final int [] windowsByGc = new int [windows];

        while ((ref = refFile.nextSequence()) != null) {
            accumulateWindowsByGc(ref.getBases(), windowSize, windowsByGc);
        }

        return windowsByGc;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Calculate the number of windows at each gc content bin for each named group of contigs
    // in a single pass over the reference. contigToGroupNames maps a contig name
    // to the groups it belongs to, as returned by parseContigGroups. Contigs that are not
    // listed in any group are skipped entirely (their bases are never scanned). Every contig
    // named in contigToGroupNames must actually be present in the reference; this can only
    // fail to hold if REFERENCE_SEQUENCE doesn't match the BAM's sequence dictionary that
    // CONTIG_GROUP was validated against, which is a misconfiguration worth failing loudly on
    // rather than silently under-counting the affected group's windows.
    /////////////////////////////////////////////////////////////////////////////
    public static Map<String, int[]> calculateRefWindowsByGcForGroups(final int windows, final File referenceSequence,
                                                                     final int windowSize,
                                                                     final Map<String, List<String>> contigToGroupNames) {
        final Map<String, int[]> windowsByGcByGroup = new LinkedHashMap<>();
        final Set<String> missingContigs = new LinkedHashSet<>(contigToGroupNames.keySet());

        final ReferenceSequenceFile refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceSequence);
        ReferenceSequence ref;

        while ((ref = refFile.nextSequence()) != null) {
            final List<String> groupNames = contigToGroupNames.get(ref.getName());
            if (groupNames == null) continue;

            missingContigs.remove(ref.getName());
            for (final String groupName : groupNames) {
                accumulateWindowsByGc(ref.getBases(), windowSize, windowsByGcByGroup.computeIfAbsent(groupName, k -> new int[windows]));
            }
        }

        if (!missingContigs.isEmpty()) {
            throw new PicardException("Contig(s) named in CONTIG_GROUP were not found in REFERENCE_SEQUENCE (" +
                    referenceSequence + "): " + missingContigs + ". REFERENCE_SEQUENCE must contain every contig " +
                    "named in CONTIG_GROUP.");
        }

        return windowsByGcByGroup;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Adds the window/gc counts for a single contig's bases into target.
    /////////////////////////////////////////////////////////////////////////////
    private static void accumulateWindowsByGc(final byte[] bases, final int windowSize, final int[] target) {
        StringUtil.toUpperCase(bases);
        final int refLength = bases.length;
        final int lastWindowStart = refLength - windowSize;

        final CalculateGcState state = new GcBiasUtils().new CalculateGcState();

        for (int i = 1; i < lastWindowStart; ++i) {
            final int windowEnd = i + windowSize;
            final int gcBin = calculateGc(bases, i, windowEnd, state);
            if (gcBin != -1) target[gcBin]++;
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    // Parses CONTIG_GROUP arguments, each of the form 'group_name:contig1,contig2,...',
    // into the contig name -> group names map used to bin reads and reference windows.
    // A contig may appear in more than one group. Every contig named is validated against
    // dict, when dict is non-null.
    /////////////////////////////////////////////////////////////////////////////
    public static Map<String, List<String>> parseContigGroups(final List<String> contigGroupSpecs, final SAMSequenceDictionary dict) {
        final Map<String, List<String>> contigToGroupNames = new LinkedHashMap<>();
        final Set<String> seenGroupNames = new LinkedHashSet<>();

        for (final String spec : contigGroupSpecs) {
            final int colon = spec.indexOf(':');
            if (colon < 0) {
                throw new PicardException("Malformed CONTIG_GROUP '" + spec + "'. Expected 'group_name:contig1,contig2,...'.");
            }

            final String groupName = spec.substring(0, colon).trim();
            if (groupName.isEmpty()) {
                throw new PicardException("Missing group name in CONTIG_GROUP '" + spec + "'.");
            }
            if (!seenGroupNames.add(groupName)) {
                throw new PicardException("CONTIG_GROUP name '" + groupName + "' was given more than once.");
            }

            int contigCount = 0;
            for (final String token : spec.substring(colon + 1).split(",")) {
                final String contig = token.trim();
                if (contig.isEmpty()) continue;
                if (dict != null && dict.getSequence(contig) == null) {
                    throw new PicardException("Contig '" + contig + "' in CONTIG_GROUP '" + groupName +
                            "' was not found in the sequence dictionary.");
                }
                final List<String> groupsForContig = contigToGroupNames.computeIfAbsent(contig, k -> new ArrayList<>());
                if (!groupsForContig.contains(groupName)) {
                    groupsForContig.add(groupName);
                    contigCount++;
                }
            }
            if (contigCount == 0) {
                throw new PicardException("CONTIG_GROUP '" + groupName + "' lists no contigs.");
            }
        }

        return contigToGroupNames;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Calculate all the GC values for all windows
    /////////////////////////////////////////////////////////////////////////////
    public static byte [] calculateAllGcs(final byte[] refBases, final int lastWindowStart, final int windowSize) {

        final CalculateGcState state = new GcBiasUtils().new CalculateGcState();

        final int refLength = refBases.length;
        final byte[] gc = new byte[refLength + 1];

        for (int i = 1; i < lastWindowStart; ++i) {
            final int windowEnd = i + windowSize;
            final int windowGc = calculateGc(refBases, i, windowEnd, state);
            gc[i] = (byte) windowGc;
        }

        return gc;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Keeps track of current GC calculation state
    /////////////////////////////////////////////////////////////////////////////
    class CalculateGcState {
        boolean init = true;
        int nCount;
        int gcCount;
        byte priorBase;
    }
}
