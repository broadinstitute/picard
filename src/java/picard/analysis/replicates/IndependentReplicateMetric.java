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

package picard.analysis.replicates;

/**
 * A class to store information relevant for biological rate estimation
 *
 * @author Yossi Farjoun
 */
public class IndependentReplicateMetric extends MergeableMetricBase {

    // the count of sites used
    @MergeByAdding
    public Integer nSites = 0;
    // the count of sites in which a third allele was found
    @MergeByAdding
    public Integer nThreeAllelesSites = 0;
    // the total number of reads over the het sites
    @MergeByAdding
    public Integer nTotalReads = 0;
    // the number of duplicate sets examined
    @MergeByAdding
    public Integer nDuplicateSets = 0;
    // the number of sets of size exactly 3 found
    @MergeByAdding
    public Integer nExactlyTriple = 0;
    // the number of sets of size exactly 2 found
    @MergeByAdding
    public Integer nExactlyDouble = 0;
    // the number of reads in duplicate of sizes greater than 3
    @MergeByAdding
    public Integer nReadsInBigSets = 0;
    // the number of doubletons where the two reads had different bases in the locus
    @MergeByAdding
    public Integer nDifferentAllelesBiDups = 0;
    // the number of doubletons where the two reads matched the reference
    @MergeByAdding
    public Integer nReferenceAllelesBiDups = 0;
    // the number of doubletons where the two reads matched the alternate
    @MergeByAdding
    public Integer nAlternateAllelesBiDups = 0;
    // the number of tripletons where at least one of the reads didn't match either allele of the het site
    @MergeByAdding
    public Integer nDifferentAllelesTriDups = 0;
    // the number of tripletons where the two reads had different bases in the locus
    @MergeByAdding
    public Integer nMismatchingAllelesBiDups = 0;
    // the number of tripletons where the two reads matched the reference
    @MergeByAdding
    public Integer nReferenceAllelesTriDups = 0;
    // the number of tripletons where the two reads matched the alternate
    @MergeByAdding
    public Integer nAlternateAllelesTriDups = 0;
    // the number of tripletons where at least one of the reads didn't match either allele of the het site
    @MergeByAdding
    public Integer nMismatchingAllelesTriDups = 0;
    // the number of reference alleles in the reads;
    @MergeByAdding
    public Integer nReferenceReads = 0;
    // the number of alternate alleles in the reads;
    @MergeByAdding
    public Integer nAlternateReads = 0;
    // the number of UMIs that are different within Bi-sets that come from different Alleles
    @MergeByAdding
    public Integer nMismatchingUMIsInDiffBiDups = 0;
    // the number of UMIs that are match within Bi-sets that come from different Alleles
    @MergeByAdding
    public Integer nMatchingUMIsInDiffBiDups = 0;
    // the number of UMIs that are different within Bi-sets that come from the same Alleles
    @MergeByAdding
    public Integer nMismatchingUMIsInSameBiDups = 0;
    // the number of UMIs that are match within Bi-sets that come from the same Alleles
    @MergeByAdding
    public Integer nMatchingUMIsInSameBiDups = 0;
    // the number of bi-sets with mismatching UMIs and same orientation
    @MergeByAdding
    public Integer nMismatchingUMIsInCoOrientedBiDups = 0;
    // the number of bi-sets with mismatching UMIs and opposite orientation
    @MergeByAdding
    public Integer nMismatchingUMIsInContraOrientedBiDups = 0;
    // the number of sets where the UMIs had poor quality bases and were not used for any comparisons.
    @MergeByAdding
    public Integer nBadBarcodes = 0;
    // the number of sets where the UMIs had good quality bases and were used for any comparisons.
    @MergeByAdding
    public Integer nGoodBarcodes = 0;
    // the rate of heterogeneity within doubleton sets
    @NoMergingIsDerived
    public Double biSiteHeterogeneityRate = 0.0;
    // the rate of heterogeneity within tripleton sets
    @NoMergingIsDerived
    public Double triSiteHeterogeneityRate = 0.0;
    // the rate of homogeneity within doubleton sets
    @NoMergingIsDerived
    public Double biSiteHomogeneityRate = 0.0;
    // the rate of homogeneity within tripleton sets
    @NoMergingIsDerived
    public Double triSiteHomogeneityRate = 0.0;
    //The biological duplication rate calculated from doublton sets
    @NoMergingIsDerived
    public Double independentReplicationRateFromBiDups = 0.0;
    //The biological duplication rate calculated from tripleton sets
    @NoMergingIsDerived
    public Double independentReplicationRateFromTriDups = 0.0;
    // when the alleles are different, we know that this is a biological duplication, thus we expect nearly all
    // the UMIs to be different (allowing for equality due to chance). So we expect this to be near 1.
    @NoMergingIsDerived
    public Double pSameUmiInIndependentBiDup = 0.0;
    //when the UMIs mismatch, we expect about the same number of different alleles as the same (assuming
    //that different UMI implied biological duplicate. thus, this value should be near 0.5
    @NoMergingIsDerived
    public Double pSameAlleleWhenMismatchingUmi = 0.0;
    // given the UMIs one can estimate the rate of biological duplication directly, as this would be the
    // rate of having different UMIs in all duplicate sets. This is only a good estimate if the assumptions hold, for example if pSameUmiInIndependentBiDup is near 1.
    @NoMergingIsDerived
    public Double independentReplicationRateFromUmi = 0.0;
    //an estimate of the duplication rate that is based on the duplicate sets we observed
    @NoMergingIsDerived
    public Double replicationRateFromReplicateSets = 0.0;

    @Override
    public void calculateDerivedFields() {
        // In doubleton sets, the rate of different alleles over het sites is half the replication rate,
        biSiteHeterogeneityRate = nDifferentAllelesBiDups / (double) (nDifferentAllelesBiDups + nAlternateAllelesBiDups + nReferenceAllelesBiDups);
        biSiteHomogeneityRate = 1 - biSiteHeterogeneityRate;

        this.independentReplicationRateFromBiDups = 2 * biSiteHeterogeneityRate;

        // in tripleton sets, the calculation is a little more complicated....see below
        triSiteHeterogeneityRate = nDifferentAllelesTriDups / (double) (nDifferentAllelesTriDups + nAlternateAllelesTriDups + nReferenceAllelesTriDups);
        triSiteHomogeneityRate = 1 - triSiteHeterogeneityRate;
        independentReplicationRateFromTriDups = 2 * (1 - Math.sqrt(triSiteHomogeneityRate));

        // Some more metric collection here:
        pSameUmiInIndependentBiDup = nMatchingUMIsInDiffBiDups / (double) (nMismatchingUMIsInDiffBiDups + nMatchingUMIsInDiffBiDups);
        pSameAlleleWhenMismatchingUmi = nMismatchingUMIsInSameBiDups / (double) (nMismatchingUMIsInSameBiDups + nMismatchingUMIsInDiffBiDups);
        independentReplicationRateFromUmi = (nMismatchingUMIsInDiffBiDups + nMismatchingUMIsInSameBiDups) / (double) nExactlyDouble;

        final int numberOfBigSets = nDuplicateSets - nExactlyDouble - nExactlyDouble;
        replicationRateFromReplicateSets = (nExactlyDouble + nExactlyTriple * 2 + nReadsInBigSets - numberOfBigSets) / (double) nTotalReads;
    }
}

/*
Explanation of calculation of independent replication rate from heterozygosity rate (within triplicate sets):

We assume the there are two types of replication events:

- those that are "independent", such that we just happen to get 2 fragments from the exact same region
(These get a random allele so effectively change the allele with probability 0.5), and
- those that are "artifactual" (do not change the allele)

We skip sets that have unexpected alleles as they do not fit our model. In the following we use the term "duplicates" to indicate that
the read-pairs would be marked as duplicates, not that they actually are technical duplicates.

To reach a triplicate set, 2 replication events are required so there are the following options, assuming
that the independent replication rate is x (thus the artifactual is 1-x). We assume a diploid organism with no bias towards either allele
in heterozygous sites, so an idependent replication will result in the other alleles in half the cases:

0 -> 0, 1 happens with probability x/2, therefore
0 -> 0, 0 happens with probability 1-x/2

(This is the explanation that is required for calculating the independent replication rate from doubleton sets...quite simpler)

Without loss of generality we assume that we "start" with allele 0 and that 1 is the other allele.

Each of the resulting alleles ("First" or "second") can replicate (each with probability 0.5) so we get:

from first row:
0=>1 1 (0 replicate to a 1) with probability  x/2 * 0.5 * x/2
0=>0 1 (0 replicate to a 0) with probability  x/2 * 0.5 * (1-x/2)
=====subtotal = x/4

0 1=>0 (1 replicate to a 0) with probability x/2 * 0.5 * x/2
0 1=>1 (1 replicate to a 1) with probability x/2 * 0.5 * (1-x/2)
=====subtotal = x/4

from second row:
0=>1 0 (first 0 replicate to a 1) with probability (1-x/2) * 0.5 * x/2
0=>0 0 (first 0 replicate to a 0) with probability (1-x/2) * 0.5 * (1-x/2)   <======= Homogeneous set
=====subtotal = (1-x/2)/2

0 0=>1 (second 0 replicate to a 1) with probability (1-x/2) * 0.5 * x/2
0 0=>0 (second 0 replicate to a 0) with probability (1-x/2) * 0.5 * (1-x/2)  <======= Homogeneous set
=====subtotal = (1-x/2)/2


total is x/2 (from first two sub-totals) + (1-x/2) (from last two sub-totals) = 1

We differentiate between a heterogeneous result (with 0's and 1's in the set) and homogeneous results
(all zeros, since we assumed WLOG that we start with 0)
The probability of a homogeneous set is therefore the sum of the two only homogeneous results

P(hom) = (1-x/2) * 0.5 * (1-x/2) + (1-x/2) * 0.5 * (1-x/2) = (1-x/2)^2 = y

where y = P(hom) is the rate of homogeneity within triplicate sets.

Solving for x we find that 2 * ( 1 - sqrt(y) ) = x.

*/
