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

package picard.sam.markduplicates;

/**
 * This is the testing class for the QuerySortedMarkDuplicates (it's the same CLP, only it
 * takes a query-sorted bam as input.) The difference in the tests themselves is that
 * supplementary and secondary reads will be marked the same as their primary alignments,
 * and that unmapped reads with mapped mates will be marked according to their mapped mate.
 */
public class QuerySortedMarkDuplicatesTest extends MarkDuplicatesTest {

    @Override
    protected boolean markUnmappedRecordsLikeTheirMates() {
        return true;
    }

    @Override
    protected boolean markSecondaryAndSupplementaryRecordsLikeTheCanonical() { return true; }

    @Override
    protected AbstractMarkDuplicatesCommandLineProgramTester getTester() { return new QuerySortedMarkDuplicatesTester(); }
}
