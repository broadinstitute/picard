/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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
package net.sf.picard.util;

import net.sf.picard.PicardException;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMException;
import net.sf.samtools.util.CoordMath;

import java.util.LinkedList;
import java.util.List;

/**
 * @author alecw@broadinstitute.org
 */
public class CigarUtil {
    /** adjust the cigar based on adapter clipping.
     * TODO: If there is hard clipping at the end of the input CIGAR, it is lost.  It should not be. 
     * *
     * @param clipFrom       1-based position where the clipping starts
     * @param oldCigar       The existing unclipped cigar
     * @return               New adjusted list of cigar elements
     */
    // package visible so can be unit-tested
    public static List<CigarElement> softClipEndOfRead(final int clipFrom, final List<CigarElement> oldCigar) {
        final int clippedBases = (int)CoordMath.getLength(clipFrom, Cigar.getReadLength(oldCigar));
        List<CigarElement> newCigar = new LinkedList<CigarElement>();
        int pos = 1;

        for (CigarElement c : oldCigar) {
            // Distinguish two cases:
            //	c occurs before the clipped region
            //	c is adjacent to or straddles the boundary between clipped and unclipped region.
            //  c never occurs after the clipped region; clipped region is always at the end

            final CigarOperator op = c.getOperator();
            final int length = op.consumesReadBases()? c.getLength() : 0;
            final int endPos = pos + length - 1;  // same as pos on next iteration

            if (endPos < (clipFrom - 1)) {
                // handle elements before clip position (just copy them)
                newCigar.add(c);

            } else if (endPos >= (clipFrom - 1)) {
                // handle adjacent or straddling element
                elementStraddlesClippedRead(newCigar, c,
                        (clipFrom -1) - (pos -1) , clippedBases);
                break;
            }

            pos = endPos + 1;      // update pos for next iteration
        } // end loop over cigar elements
        return newCigar;
    }

    // a cigar element occurs in the middle of an adapter clipping
    static private void elementStraddlesClippedRead(List<CigarElement> newCigar, CigarElement c,
                                                    int relativeClippedPosition,
                                                    int clippedBases){
        final CigarOperator op = c.getOperator();
        int clipAmount = clippedBases;
        if (op.consumesReadBases()){
            if (op.consumesReferenceBases() & relativeClippedPosition > 0){
               newCigar.add(new CigarElement(relativeClippedPosition, op));
            }
            if (!op.consumesReferenceBases()){
                clipAmount = clippedBases + relativeClippedPosition;
            }
        } else if (relativeClippedPosition != 0){
                throw new SAMException("Unexpected non-0 relativeClippedPosition " + relativeClippedPosition);
        }
        newCigar.add(new CigarElement(clipAmount, CigarOperator.S));  // S is always last element
    }

}
