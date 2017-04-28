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

import htsjdk.samtools.util.Histogram;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;

import java.util.ArrayList;

import static org.testng.Assert.assertTrue;
import static picard.sam.markduplicates.EstimateLibraryComplexity.PairedReadSequence;

public class ElcIdenticalBasesDuplicatesFinderTest {

    static final double MAX_DIFF_RATE = 0.03;
    static final int MAX_READ_LENGTH = Integer.MAX_VALUE;
    static final int MIN_IDENTICAL_BASES = 5;
    static final boolean USE_BARCODES = false;
    static final OpticalDuplicateFinder OPTICAL_DUPLICATE_FINDER = new OpticalDuplicateFinder();

    private ElcDuplicatesFinder duplicatesFinder = new ElcIdenticalBasesDuplicatesFinder(
            MAX_DIFF_RATE,
            MAX_READ_LENGTH,
            MIN_IDENTICAL_BASES,
            USE_BARCODES,
            OPTICAL_DUPLICATE_FINDER
    );

    @DataProvider(name = "fillHistogramDataProvider")
    public Object[][] fillHistogramDataProvider() {
        return new Object[][]{
                // empty dups, that mean increment first bin
                {
                        duplicatesFinder,
                        new Histogram<>(),
                        new Histogram<>(),
                        generatePairedReadSequence(false),
                        generatePairedReadSequences(0, false),
                        1,
                        0,
                        0
                },
                // 10 dups + 1 identical paired-read, that mean increment 11 bin
                {
                        duplicatesFinder,
                        new Histogram<>(),
                        new Histogram<>(),
                        generatePairedReadSequence(false),
                        generatePairedReadSequences(10, false),
                        11,
                        0,
                        0
                },
                // 10 dups + 1 identical paired-read, that mean increment 11 bin and isOpticalDuplicates is true, then
                // increment opticalHisto 11 bin 10 times
                {
                        duplicatesFinder,
                        new Histogram<>(),
                        new Histogram<>(),
                        generatePairedReadSequence(true),
                        generatePairedReadSequences(10, true),
                        11,
                        11,
                        10
                }
        };
    }

    @DataProvider(name = "searchDuplicatesDataProvider")
    public Object[][] searchDuplicatesDataProvider() {
        return new Object[][]{
                // empty dups, that mean increment first bin
                {
                        duplicatesFinder,
                        new Histogram<>(),
                        new Histogram<>(),
                        generatePairedReadSequences(1, false),
                        1,
                        0,
                        0
                },
                // 10 dups, that mean increment 10 bin
                {
                        duplicatesFinder,
                        new Histogram<>(),
                        new Histogram<>(),
                        generatePairedReadSequences(10, false),
                        10,
                        0,
                        0
                },
                // 10 dups, that mean increment 10 bin and isOpticalDuplicates is true, then
                // increment opticalHisto 10 bin 9 times
                {
                        duplicatesFinder,
                        new Histogram<>(),
                        new Histogram<>(),
                        generatePairedReadSequences(10, true),
                        10,
                        10,
                        9
                },
                // 10 paired-reads but only 9 dups, that mean increment 9 bin and isOpticalDuplicates is true, then
                // increment opticalHisto 9 bin 8 times
                {
                        duplicatesFinder,
                        new Histogram<>(),
                        new Histogram<>(),
                        generateSeqsWithNoDup(10, true),
                        9,
                        9,
                        8
                }
        };
    }

    @Test(dataProvider = "fillHistogramDataProvider")
    public void testFillHistogram(ElcDuplicatesFinder duplicatesFinder,
                                  Histogram<Integer> duplicationHisto,
                                  Histogram<Integer> opticalHisto,
                                  PairedReadSequence prs,
                                  ArrayList<PairedReadSequence> dupes,
                                  int dupHistoIndex,
                                  int optHistoIndex,
                                  int optHistoValue) throws Exception {
        duplicatesFinder.fillHistogram(duplicationHisto, opticalHisto, prs, dupes);
        assertTrue(duplicationHisto.get(dupHistoIndex).getValue() == 1);
        if (optHistoIndex > 0) {
            assertTrue(opticalHisto.get(optHistoIndex).getValue() == optHistoValue);
        }
    }

    @Test(dataProvider = "searchDuplicatesDataProvider")
    public void testSearchDuplicates(ElcDuplicatesFinder duplicatesFinder,
                                     Histogram<Integer> duplicationHisto,
                                     Histogram<Integer> opticalHisto,
                                     ArrayList<PairedReadSequence> dupes,
                                     int dupHistoIndex,
                                     int optHistoIndex,
                                     int optHistoValue) throws Exception {
        duplicatesFinder.searchDuplicates(dupes, duplicationHisto, opticalHisto);
        assertTrue(duplicationHisto.get(dupHistoIndex).getValue() == 1);
        if (optHistoIndex > 0) {
            assertTrue(opticalHisto.get(optHistoIndex).getValue() == optHistoValue);
        }
    }

    protected ArrayList<PairedReadSequence> generatePairedReadSequences(int seqsSize, boolean isOpticalDuplicates) {
        ArrayList<PairedReadSequence> seq = new ArrayList<>(seqsSize);
        for (int i = 0; i < seqsSize; i++) {
            seq.add(generatePairedReadSequence(isOpticalDuplicates));
        }
        return seq;
    }

    protected ArrayList<PairedReadSequence> generateSeqsWithNoDup(int seqsSize, boolean isOpticalDuplicates) {
        ArrayList<PairedReadSequence> seq = generatePairedReadSequences(seqsSize, isOpticalDuplicates);
        PairedReadSequence prs = seq.get(seq.size() - 1);
        changeReadContent(prs.read1);
        changeReadContent(prs.read2);
        return seq;
    }

    protected PairedReadSequence generatePairedReadSequence(boolean isOpticalDuplicates) {
        PairedReadSequence prs = new PairedReadSequence();
        prs.read1 = new byte[100];
        prs.read2 = new byte[100];
        if (isOpticalDuplicates) {
            short readGroup = 1;
            short tile = Short.MAX_VALUE;
            prs.setReadGroup(readGroup);
            prs.setTile(tile);
            prs.setX(127);
            prs.setY(255);
        }
        return prs;
    }

    private void changeReadContent(byte[] read) {
        for (int i = MIN_IDENTICAL_BASES; i <= MIN_IDENTICAL_BASES + read.length * MAX_DIFF_RATE; i++) {
            read[i] += 1;
        }
    }
}