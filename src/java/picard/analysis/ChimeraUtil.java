/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.SamPairUtil.PairOrientation;

import java.util.EnumSet;
import java.util.Set;

public class ChimeraUtil {
    public static int DEFAULT_INSERT_SIZE_LIMIT = 100000;
    public static Set<PairOrientation> DEFAULT_EXPECTED_ORIENTATIONS = EnumSet.of(PairOrientation.FR);

    /**
     * Checks whether the given read is part of a chimeric pair.
     * Note that this method returns false if the read is unpaired or if either end of the pair is unmapped.
     *
     * @param rec                    the read
     * @param maxInsertSize          max insert size to be considered non-chimeric
     * @param expectedOrientations   set of orientations that are not chimeric; must not ne null
     * @return true if this record is part of a chimeric read pair, false otherwise
     */
    public static boolean isChimeric(final SAMRecord rec, final int maxInsertSize, final Set<PairOrientation> expectedOrientations) {
        return !(rec.getReadUnmappedFlag() || !rec.getReadPairedFlag() || rec.getMateUnmappedFlag()) &&  // the read pair needs to be mapped
                (Math.abs(rec.getInferredInsertSize()) > maxInsertSize ||
                !rec.getReferenceIndex().equals(rec.getMateReferenceIndex()) ||
                expectedOrientations.contains(SamPairUtil.getPairOrientation(rec)) ||
                rec.getAttribute("SA") != null);
    }

    /**
     * Checks whether the given read is part of a chimeric pair.
     * Note that this method returns false if either end of the pair is unmapped.
     *
     * @param r1                     first read of the pair
     * @param r2                     second read of the pair
     * @param maxInsertSize          max insert size to be considered non-chimeric
     * @param expectedOrientations   set of orientations that are not chimeric
     * @return true if this pair is chimeric, false otherwise
     */
    public static boolean isChimeric(final SAMRecord r1, final SAMRecord r2, final int maxInsertSize, final Set<PairOrientation> expectedOrientations) {
        return !(r1.getReadUnmappedFlag() || r2.getMateUnmappedFlag()) &&  // the read pair needs to be mapped
                (Math.abs(r1.getInferredInsertSize()) > maxInsertSize ||
                !r1.getReferenceIndex().equals(r2.getReferenceIndex()) ||
                expectedOrientations.contains(SamPairUtil.getPairOrientation(r1)) ||
                r1.getAttribute("SA") != null ||
                r2.getAttribute("SA") != null);
    }
}
