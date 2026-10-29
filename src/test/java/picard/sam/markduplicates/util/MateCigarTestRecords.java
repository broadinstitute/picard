/*
 * The MIT License
 *
 * Copyright (c) 2026 The Broad Institute
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

package picard.sam.markduplicates.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;

import java.util.Arrays;

/**
 * Builds SAM records with mate cigars for the tests of {@link ReadEndsForMateCigar}, {@link MarkQueue} and the
 * mate-cigar location sets. The header has two sequences and three read groups: rgA and rgC in library libA,
 * rgB in library libB.
 */
final class MateCigarTestRecords {

    final SAMFileHeader header = new SAMFileHeader();

    MateCigarTestRecords() {
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        header.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.addSequence(new SAMSequenceRecord("chr2", 100000));
        header.addReadGroup(readGroup("rgA", "libA"));
        header.addReadGroup(readGroup("rgB", "libB"));
        header.addReadGroup(readGroup("rgC", "libA"));
    }

    private static SAMReadGroupRecord readGroup(final String id, final String library) {
        final SAMReadGroupRecord rg = new SAMReadGroupRecord(id);
        rg.setLibrary(library);
        rg.setSample("sample");
        return rg;
    }

    /** An unpaired mapped read whose bases all have the given quality. */
    SAMRecord fragment(final String name, final String readGroup, final int referenceIndex, final int start,
                       final boolean negativeStrand, final String cigar, final int quality) {
        final SAMRecord rec = new SAMRecord(header);
        rec.setReadName(name);
        if (readGroup != null) rec.setAttribute("RG", readGroup);
        rec.setReferenceIndex(referenceIndex);
        rec.setAlignmentStart(start);
        rec.setReadNegativeStrandFlag(negativeStrand);
        rec.setCigarString(cigar);
        rec.setMappingQuality(60);
        final int length = TextCigarCodec.decode(cigar).getReadLength();
        final byte[] bases = new byte[length];
        Arrays.fill(bases, (byte) 'A');
        rec.setReadBases(bases);
        final byte[] qualities = new byte[length];
        Arrays.fill(qualities, (byte) quality);
        rec.setBaseQualities(qualities);
        return rec;
    }

    /** One end of a pair whose mate is mapped, with the MC tag set from the mate's cigar. */
    SAMRecord pairEnd(final String name, final String readGroup, final int referenceIndex, final int start,
                      final boolean negativeStrand, final String cigar, final boolean firstOfPair,
                      final int mateReferenceIndex, final int mateStart, final boolean mateNegativeStrand,
                      final String mateCigar, final int quality) {
        final SAMRecord rec = fragment(name, readGroup, referenceIndex, start, negativeStrand, cigar, quality);
        rec.setReadPairedFlag(true);
        rec.setFirstOfPairFlag(firstOfPair);
        rec.setSecondOfPairFlag(!firstOfPair);
        rec.setMateReferenceIndex(mateReferenceIndex);
        rec.setMateAlignmentStart(mateStart);
        rec.setMateNegativeStrandFlag(mateNegativeStrand);
        rec.setAttribute("MC", mateCigar);
        return rec;
    }

    /** A paired read whose mate is unmapped and placed at the read's own position. */
    SAMRecord mateUnmapped(final String name, final String readGroup, final int referenceIndex, final int start,
                           final boolean negativeStrand, final String cigar, final boolean firstOfPair, final int quality) {
        final SAMRecord rec = fragment(name, readGroup, referenceIndex, start, negativeStrand, cigar, quality);
        rec.setReadPairedFlag(true);
        rec.setFirstOfPairFlag(firstOfPair);
        rec.setSecondOfPairFlag(!firstOfPair);
        rec.setMateUnmappedFlag(true);
        rec.setMateReferenceIndex(referenceIndex);
        rec.setMateAlignmentStart(start);
        return rec;
    }

    /** Both ends of a pair, first end first. */
    SAMRecord[] pair(final String name, final String readGroup, final int referenceIndex,
                     final int start1, final boolean negative1, final String cigar1,
                     final int start2, final boolean negative2, final String cigar2, final int quality) {
        return new SAMRecord[]{
                pairEnd(name, readGroup, referenceIndex, start1, negative1, cigar1, true, referenceIndex, start2, negative2, cigar2, quality),
                pairEnd(name, readGroup, referenceIndex, start2, negative2, cigar2, false, referenceIndex, start1, negative1, cigar1, quality)
        };
    }
}
