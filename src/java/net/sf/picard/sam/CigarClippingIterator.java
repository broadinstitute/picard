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
package net.sf.picard.sam;

import net.sf.picard.util.CigarUtil;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloseableIterator;

import java.util.List;

/**
 * @author ktibbett@broadinsitute.org
 */
public class CigarClippingIterator implements CloseableIterator<SAMRecord> {

    private final CloseableIterator<SAMRecord> iterator;

    public CigarClippingIterator(CloseableIterator<SAMRecord> iterator) {
        this.iterator = iterator;
    }

    @Override
    public SAMRecord next() {
        SAMRecord rec = iterator.next();
        if (!rec.getReadUnmappedFlag()) {
            SAMSequenceRecord refseq = rec.getHeader().getSequence(rec.getReferenceIndex());
            if (rec.getAlignmentEnd() > refseq.getSequenceLength()) {
                // 1-based index of first base in read to clip.
                int clipFrom = refseq.getSequenceLength() - rec.getAlignmentStart() + 1;
                List<CigarElement> newCigarElements  = CigarUtil.softClipEndOfRead(clipFrom, rec.getCigar().getCigarElements());
                rec.setCigar(new Cigar(newCigarElements));
            }
        }
        return rec;
    }

    @Override
    public boolean hasNext() {
        return iterator.hasNext();
    }

    @Override
    public void remove() {
        iterator.remove();
    }

    @Override
    public void close() {
        iterator.close();
    }
}
