/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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
package picard.sam;

import com.google.common.cache.Cache;
import com.google.common.cache.CacheBuilder;
import htsjdk.samtools.*;

import java.io.File;
import java.io.Serializable;
import java.util.concurrent.ExecutionException;

/**
 * Comparator for sorting SAMRecords by chr + location + queryname.
 *
 */
public class YossiSortComparator implements SAMRecordComparator, Serializable {

    protected class Coordinate {
        public String chr;
        public int pos;
        public Coordinate(String chr, int pos) {
            this.chr = chr; this.pos = pos;
        }
    }

    private static final long serialVersionUID = 1L;


    protected Cache<String, Coordinate> primaryCache = CacheBuilder.newBuilder()
            .maximumSize(1000)
            .build();

    protected ThreadLocal<SamReader> queryReader;

    public YossiSortComparator(SamReaderFactory readerFactory, File refSeq, File bam) {
         queryReader = ThreadLocal.withInitial(() -> readerFactory.referenceSequence(refSeq).open(bam));
    }

    protected String getRefFromSATag(String saTag) {
        //SA:Z:11,6482131,-,44S107M,60,0;<next>
        return saTag.split(";")[0].split(",")[0];
    }

    protected int getPosFromSATag(String saTag) {
        return Integer.parseInt(saTag.split(";")[0].split(",")[1]);
    }

    protected boolean recIsWeird(SAMRecord rec) {
        return rec.isSecondaryOrSupplementary()
            && rec.getMateUnmappedFlag()
            && rec.getReferenceName().equals(rec.getMateReferenceName())
            && rec.getAlignmentStart() == rec.getMateAlignmentStart();
    }

    //SA:Z:11,6482131,-,44S107M,60,0

    protected String getFirstRefName(final SAMRecord rec) {
        if (rec.getFirstOfPairFlag()) {
            if(rec.isSecondaryOrSupplementary()) {
                return getRefFromSATag(rec.getStringAttribute("SA"));
            } else {
                return rec.getReferenceName();
            }
        } else {
            return rec.getMateReferenceName();
        }
    }

    protected int getFirstPos(final SAMRecord rec) {
        if (rec.getFirstOfPairFlag()) {
            if(rec.isSecondaryOrSupplementary()) {
                return getPosFromSATag(rec.getStringAttribute("SA"));
            }
            return rec.getAlignmentStart();
        } else {
            return rec.getMateAlignmentStart();
        }
    }

    protected Coordinate getPrimary1(SAMRecord rec) {
        try {
            return primaryCache.get(rec.getReadName(), () -> getPrimary1Internal(rec));
        } catch (ExecutionException e) {
            e.printStackTrace();
            System.exit(1);
        }
        assert(false); //shouldn't get here
        return null;
    }

    protected Coordinate getPrimary1Internal(SAMRecord rec) {
        /*
         * I found mysterious pair2 records with mate unmapped.
         * Its RNEXT and PNEXT pointed to itself, though it did have an SA tag.
         * Detect this, and switcheroo us with
         */
        //its mate info points to itself(?!), though it does have an SA tag.
        //look it up with a query
        if(recIsWeird(rec)) {
            String primaryRef = getRefFromSATag(rec.getStringAttribute("SA"));
            long primaryPos = getPosFromSATag(rec.getStringAttribute("SA"));
            try (
                    SAMRecordIterator iter = queryReader.get().queryAlignmentStart(primaryRef, (int) primaryPos)
            ) {
                while(iter.hasNext()) {
                    SAMRecord pRec = iter.next();
                    if (pRec.getReadName().equals(rec.getReadName())) {
                        return new Coordinate(getFirstRefName(pRec), getFirstPos(pRec));
                    }
                }
            }
            throw new RuntimeException("weird record couldn't find its primary mate");
        } else {
            return new Coordinate(getFirstRefName(rec), getFirstPos(rec));
        }
    }

    @Override
    public int compare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        Coordinate s1Coord = getPrimary1(samRecord1);
        Coordinate s2Coord = getPrimary1(samRecord2);

        if (!s1Coord.chr.equals(s2Coord.chr)) {
            return s1Coord.chr.compareTo(s2Coord.chr);
        } else if (s1Coord.pos != s2Coord.pos) {
            return s1Coord.pos - s2Coord.pos;
        } else {
            return samRecord1.getReadName().compareTo(samRecord2.getReadName());
        }
    }

    /**
     * Less stringent compare method than the regular compare.  If the two records
     * are equal enough that their ordering in a sorted SAM file would be arbitrary,
     * this method returns 0.
     *
     * @return negative if samRecord1 < samRecord2,  0 if equal, else positive
     */
    @Override
    public int fileOrderCompare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        return compare(samRecord1, samRecord2);
    }
}
