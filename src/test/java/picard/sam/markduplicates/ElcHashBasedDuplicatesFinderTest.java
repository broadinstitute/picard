package picard.sam.markduplicates;

import htsjdk.samtools.util.Histogram;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;

import static picard.sam.markduplicates.EstimateLibraryComplexity.PairedReadSequence;

public class ElcHashBasedDuplicatesFinderTest extends ElcIdenticalBasesDuplicatesFinderTest {

    private ElcDuplicatesFinder duplicatesFinder = new ElcHashBasedDuplicatesFinder(
            MAX_DIFF_RATE,
            MAX_READ_LENGTH,
            MIN_IDENTICAL_BASES,
            OPTICAL_DUPLICATE_FINDER
    );

    @DataProvider(name = "fillHistogramDataProvider")
    public Object[][] fillHistogramDataProvider() {
        return new Object[][]{
                // empty dups + 1 identical paired-read, that mean increment first bin
                {
                        duplicatesFinder,
                        new Histogram<>(),
                        new Histogram<>(),
                        generatePrs(false),
                        generateSeqs(0, false),
                        1,
                        0,
                        0
                },
                // 10 dups + 1 identical paired-read, that mean increment 11 bin
                {
                        duplicatesFinder,
                        new Histogram<>(),
                        new Histogram<>(),
                        generatePrs(true),
                        generateSeqs(10, false),
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
                        generatePrs(true),
                        generateSeqs(10, true),
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
                        generateSeqs(1, false),
                        1,
                        0,
                        0
                },
                // 10 dups, that mean increment 10 bin
                {
                        duplicatesFinder,
                        new Histogram<>(),
                        new Histogram<>(),
                        generateSeqs(10, false),
                        10,
                        0,
                        0
                },
                // 10 paired-reads but only 9 dupes, that mean increment 9 bin
                {
                        duplicatesFinder,
                        new Histogram<>(),
                        new Histogram<>(),
                        generateSeqsWithNoDup(10, true),
                        9,
                        0,
                        0
                },
                // 10 dups, that mean increment 10 bin and isOpticalDuplicates is true, then
                // increment opticalHisto 10 bin 9 times
                {
                        duplicatesFinder,
                        new Histogram<>(),
                        new Histogram<>(),
                        generateSeqs(10, true),
                        10,
                        10,
                        9
                },
                // 10 paired-reads but only 9 dupes, that mean increment 9 bin and isOpticalDuplicates is true, then
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
        super.testFillHistogram(
                duplicatesFinder,
                duplicationHisto,
                opticalHisto,
                prs,
                dupes,
                dupHistoIndex,
                optHistoIndex,
                optHistoValue
        );
    }

    @Test(dataProvider = "searchDuplicatesDataProvider")
    public void testSearchDuplicates(ElcDuplicatesFinder duplicatesFinder,
                                     Histogram<Integer> duplicationHisto,
                                     Histogram<Integer> opticalHisto,
                                     ArrayList<PairedReadSequence> dupes,
                                     int dupHistoIndex,
                                     int optHistoIndex,
                                     int optHistoValue) throws Exception {
        super.testSearchDuplicates(
                duplicatesFinder,
                duplicationHisto,
                opticalHisto,
                dupes,
                dupHistoIndex,
                optHistoIndex,
                optHistoValue
        );
    }

}