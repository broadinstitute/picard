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
package net.sf.samtools;

import net.sf.samtools.util.CloseableIterator;

/**
 * One crisp, informative sentence or noun phrase that explains
 * the concept modeled by the class.
 * <p/>
 * This class is [<em>not</em>] thread safe [because it is immutable].
 */
public class BamIndexValidator {

    public static int exhaustivelyTestIndex(SAMFileReader reader) { // throws Exception {
        // look at all chunk offsets in a linear index to make sure they are valid

        if (reader.hasBrowseableIndex()) {

            // content is from an existing bai file
            final CachingBAMFileIndex existingIndex = (CachingBAMFileIndex) reader.getBrowseableIndex(); // new CachingBAMFileIndex(inputBai, null);
            final int n_ref = existingIndex.getNumberOfReferences();

            int chunkCount = 0;
            int indexCount = 0;
            for (int i = 0; i < n_ref; i++) {
                BAMIndexContent content = existingIndex.getQueryResults(i);
                for (Chunk c : content.getAllChunks()) {
                    final CloseableIterator<SAMRecord> iter = reader.iterator(new BAMFileSpan(c));
                    chunkCount++;
                    BAMRecord b = null;
                    try {
                        // if (iter.hasNext()) {   // not needed since there should be something there
                        b = (BAMRecord) iter.next();
                        // }
                        iter.close();
                    } catch (Exception e) {
                        throw new SAMException("Exception in BamIndexValidator. Last good record " + b + " in chunk " + c + " chunkCount=" + chunkCount, e);
                    }
                }
                // also seek to every position in the linear index
                // final BAMRecordCodec bamRecordCodec = new BAMRecordCodec(reader.getFileHeader());
                // bamRecordCodec.setInputStream(reader.getInputStream());

                LinearIndex linearIndex = content.getLinearIndex();
                for (long l : linearIndex.getIndexEntries()) {
                    try {
                        if (l != 0) {
                            final CloseableIterator<SAMRecord> iter = reader.iterator(new BAMFileSpan(new Chunk(l, l + 1)));
                            BAMRecord b = (BAMRecord) iter.next();   // read the first record identified by the linear index
                            indexCount++;
                            iter.close();
                        }
                    } catch (Exception e) {
                        throw new SAMException("Exception in BamIndexValidator. Linear index access failure " + l + " indexCount=" + indexCount, e);
                    }

                }
            }
            return chunkCount;
            // System.out.println("Found " chunkCount + " chunks in test " + inputBai +
            // " linearIndex positions = " + indexCount);
        } // else  not a bam file with a browseable index
        //    System.err.println("No browseableIndex for reader");
        return 0;
    }

}
