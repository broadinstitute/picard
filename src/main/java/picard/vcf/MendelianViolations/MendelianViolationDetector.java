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

package picard.vcf.MendelianViolations;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import picard.pedigree.Sex;
import picard.vcf.processor.VariantProcessor;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * @author mccowan
 */

class MendelianViolationDetector implements VariantProcessor.Accumulator<MendelianViolationDetector.Result> {
    private final Set<String> SKIP_CHROMS, MALE_CHROMS, FEMALE_CHROMS;
    private final List<MendelianViolationMetrics> trios;
    private final List<Interval> parIntervals;
    private final double MIN_HET_FRACTION;
    private final int MIN_GQ;
    private final int MIN_DP;
    private final ProgressLogger logger;
    private final MendelianViolationsByFamily familyToViolations;

    MendelianViolationDetector(final Set<String> skip_chroms, final Set<String> male_chroms, final Set<String> female_chroms,
                               final double min_het_fraction, final int min_gq, final int min_dp, final List<MendelianViolationMetrics> trios,
                               final List<Interval> parIntervals, final ProgressLogger logger) {
        SKIP_CHROMS = skip_chroms;
        MALE_CHROMS = male_chroms;
        FEMALE_CHROMS = female_chroms;
        MIN_HET_FRACTION = min_het_fraction;
        MIN_GQ = min_gq;
        MIN_DP = min_dp;
        this.trios = trios;
        this.parIntervals = parIntervals;
        this.logger = logger;
        familyToViolations = new MendelianViolationsByFamily();
        }

    /** A little enum to describe the different type of Mendelian violations possible at a bi-allelic site. */
    enum MendelianViolation {
        Diploid_Denovo("Parents are both Homozygous Ref and offspring is Heterozygous."),
        HomVar_HomVar_Het("Parents are both Homozygous Variant and offspring is Heterozygous."),
        HomRef_HomVar_Hom("One parent is Homozygous Ref, the other Homozygous Variant and the offspring is Homozygous."),
        Hom_Het_Hom("One parent is Homozygous (Ref or Var), the other is Heterozygous and the offspring is incompatibly Homozygous."),
        Haploid_Denovo("Offspring variant genotype is expected to be haploid and does not match an expected parental allele."),
        Haploid_Other("Offspring reference genotype is expected to be haploid and does not match an expected parental allele."),
        Other("Other unclassified violation of allele transmission.");

        private final String description;

        MendelianViolation(final String desc) {
            this.description = desc;
        }

        public String getDescription() {
            return description;
        }
    }

    public final static String MENDELIAN_VIOLATION_KEY = "MV";
    public final static String ORIGINAL_AC = "AC_Orig";
    public final static String ORIGINAL_AF = "AF_Orig";
    public final static String ORIGINAL_AN = "AN_Orig";

    @Override
    public void accumulate(final VariantContext ctx) {
        logger.record(ctx.getContig(), ctx.getStart());

        final String variantChrom = ctx.getContig();
        final int variantPos = ctx.getStart();

        // Skip anything a little too funky
        if (ctx.isFiltered()) return;
        if (!ctx.isVariant()) return;
        if (SKIP_CHROMS.contains(variantChrom)) return;

        for (final MendelianViolationMetrics trio : trios) {
            final Genotype momGt = ctx.getGenotype(trio.MOTHER);
            final Genotype dadGt = ctx.getGenotype(trio.FATHER);
            final Genotype kidGt = ctx.getGenotype(trio.OFFSPRING);

            // if any genotype:
            // - has a non-snp allele; or
            // - lacks a reference allele
            //
            // then ignore this trio
            if (CollectionUtil.makeList(momGt, dadGt, kidGt).stream().anyMatch(gt ->
                    gt.isHetNonRef() ||
                            Stream.concat(Stream.of(ctx.getReference()), gt.getAlleles().stream()).anyMatch(a -> a.length() != 1 || a.isSymbolic()))) {
                continue;
            }

            // if between the trio there are more than 2 alleles including the reference, continue
            if (Stream.concat(
                    Collections.singleton(ctx.getReference()).stream(),
                    CollectionUtil.makeList(momGt, dadGt, kidGt)
                            .stream()
                            .flatMap(gt -> gt.getAlleles().stream()))
                    .collect(Collectors.toSet()).size() > 2) continue;

            // Test to make sure:
            //   1) That the site is in fact variant in the trio
            //   2) that the offspring doesn't have a really wacky het allele balance
            if (!isVariant(momGt, dadGt, kidGt)) continue;
            if (kidGt.isHet()) {
                final int[] ad = kidGt.getAD();
                if (ad == null) continue;

                final List<Integer> adOfAlleles = kidGt.getAlleles().stream()
                        .map(a->ad[ctx.getAlleleIndex(a)]).collect(Collectors.toList());
                final double minAlleleFraction = Math.min(adOfAlleles.get(0), adOfAlleles.get(1)) / (double) (adOfAlleles.get(0) + adOfAlleles.get(1));
                if (minAlleleFraction < MIN_HET_FRACTION) continue;
            }

            ///////////////////////////////////////////////////////////////
            // Determine whether the offspring should be haploid at this
            // locus and which is the parental donor of the haploid genotype
            ///////////////////////////////////////////////////////////////
            boolean haploid = false;
            Genotype haploidParentalGenotype = null;

            if (FEMALE_CHROMS.contains(variantChrom) && trio.OFFSPRING_SEX != Sex.Unknown) {
                if (trio.OFFSPRING_SEX == Sex.Female) {
                    // famale
                    haploid = false;
                } else if (isInPseudoAutosomalRegion(variantChrom, variantPos)) {
                    // male but in PAR on X, so diploid
                    haploid = false;
                } else {
                    // male, out of PAR on X, haploid
                    haploid = true;
                    haploidParentalGenotype = momGt;
                }
            }

            // the PAR on the male chromosome should be masked so that reads
            // align to the female chromosomes instead, so there's no point
            // of worrying about that here.

            if (MALE_CHROMS.contains(variantChrom)) {
                if (trio.OFFSPRING_SEX == Sex.Male) {
                    haploid = true;
                    haploidParentalGenotype = dadGt;
                } else {
                    continue;
                }
            }

            // We only want to look at sites where we have high enough confidence that the genotypes we are looking at are
            // interesting.  We want to ensure that parents are always GQ>=MIN_GQ, and that the kid is either GQ>=MIN_GQ or in the
            // case where kid is het that the phred-scaled-likelihood of being reference is >=MIN_GQ.
            if (haploid && (haploidParentalGenotype.isNoCall() || haploidParentalGenotype.getGQ() < MIN_GQ)) continue;
            if (!haploid && (momGt.isNoCall() || momGt.getGQ() < MIN_GQ || dadGt.isNoCall() || dadGt.getGQ() < MIN_GQ))
                continue;
            if (kidGt.isNoCall()) continue;
            if (momGt.isHomRef() && dadGt.isHomRef() && !kidGt.isHomRef()) {
                if (kidGt.getPL()[0] < MIN_GQ) continue;
            } else if (kidGt.getGQ() < MIN_GQ) continue;

            // Also filter on the DP for each of the samples - it's possible to miss hets when DP is too low
            if (haploid && (kidGt.getDP() < MIN_DP || haploidParentalGenotype.getDP() < MIN_DP)) continue;
            if (!haploid && (kidGt.getDP() < MIN_DP || momGt.getDP() < MIN_DP || dadGt.getDP() < MIN_DP)) continue;

            trio.NUM_VARIANT_SITES++;

            ///////////////////////////////////////////////////////////////
            // First test for haploid violations
            ///////////////////////////////////////////////////////////////
            MendelianViolation type = null;
            if (haploid) {
                if (kidGt.isHet()) continue; // Should not see heterozygous calls at haploid regions

                if (!haploidParentalGenotype.getAlleles().contains(kidGt.getAllele(0))) {
                    if (kidGt.isHomRef()) {
                        type = MendelianViolation.Haploid_Other;
                        trio.NUM_HAPLOID_OTHER++;
                    } else {
                        type = MendelianViolation.Haploid_Denovo;
                        trio.NUM_HAPLOID_DENOVO++;
                    }
                }
            }
            ///////////////////////////////////////////////////////////////
            // Then test for diploid mendelian violations
            ///////////////////////////////////////////////////////////////
            else if (isMendelianViolation(momGt, dadGt, kidGt)) {
                if (momGt.isHomRef() && dadGt.isHomRef() && !kidGt.isHomRef()) {
                    trio.NUM_DIPLOID_DENOVO++;
                    type = MendelianViolation.Diploid_Denovo;
                } else if (momGt.isHomVar() && dadGt.isHomVar() && kidGt.isHet()) {
                    trio.NUM_HOMVAR_HOMVAR_HET++;
                    type = MendelianViolation.HomVar_HomVar_Het;
                } else if (kidGt.isHom() && ((momGt.isHomRef() && dadGt.isHomVar()) || (momGt.isHomVar() && dadGt.isHomRef()))) {
                    trio.NUM_HOMREF_HOMVAR_HOM++;
                    type = MendelianViolation.HomRef_HomVar_Hom;
                } else if (kidGt.isHom() && ((momGt.isHom() && dadGt.isHet()) || (momGt.isHet() && dadGt.isHom()))) {
                    trio.NUM_HOM_HET_HOM++;
                    type = MendelianViolation.Hom_Het_Hom;
                } else {
                    trio.NUM_OTHER++;
                    type = MendelianViolation.Other;
                }
            }

            // Output a record into the family's violation VCF
            if (type != null) {
                // Create a new Context subsetted to the three samples
                final VariantContextBuilder builder = new VariantContextBuilder(ctx);
                builder.genotypes(ctx.getGenotypes().subsetToSamples(CollectionUtil.makeSet(trio.MOTHER, trio.FATHER, trio.OFFSPRING)));
                builder.attribute(MENDELIAN_VIOLATION_KEY, type.name());

                // Copy over some useful attributes from the full context
                if (ctx.hasAttribute(VCFConstants.ALLELE_COUNT_KEY))
                    builder.attribute(ORIGINAL_AC, ctx.getAttribute(VCFConstants.ALLELE_COUNT_KEY));
                if (ctx.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY))
                    builder.attribute(ORIGINAL_AF, ctx.getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY));
                if (ctx.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY))
                    builder.attribute(ORIGINAL_AN, ctx.getAttribute(VCFConstants.ALLELE_NUMBER_KEY));

                // Write out the variant record
                familyToViolations.get(trio.FAMILY_ID).add(builder.make());
            }
        }
    }

    @Override
    public Result result() {
        return new Result(trios, familyToViolations);
    }

    /** Tests to ensure that a locus is bi-allelic within the given set of samples. */
    private boolean isVariant(final Genotype... gts) {
        for (final Genotype gt : gts) {
            if (gt.isCalled() && !gt.isHomRef()) return true;
        }

        return false;
    }

    /** Tests whether the alleles of the offspring are possible given the alleles of the parents. */
    private boolean isMendelianViolation(final Genotype p1, final Genotype p2, final Genotype off) {
        final Allele offAllele1 = off.getAllele(0);
        final Allele offAllele2 = off.getAllele(1);

        if(p1.getAlleles().contains(offAllele1) && p2.getAlleles().contains(offAllele2)) return false;
        if(p2.getAlleles().contains(offAllele1) && p1.getAlleles().contains(offAllele2)) return false;

        return true;
    }

    /**
     * Tests whether the variant is within one of the pseudo-autosomal regions
     */
    private boolean isInPseudoAutosomalRegion(final String chr, final int pos) {
        for (final Interval par : parIntervals) {
            if (par.getContig().equals(chr) && pos >= par.getStart() && pos <= par.getEnd()) return true;
        }
        return false;
    }

    /** Represents the result of the work this class does. */
    static class Result {
        private final Collection<MendelianViolationMetrics> metrics;
        private final MendelianViolationsByFamily violations;

        Result(final Collection<MendelianViolationMetrics> metrics, final MendelianViolationsByFamily violations) {
            this.metrics = metrics;
            this.violations = violations;
        }
        
        public Collection<MendelianViolationMetrics> metrics() {
            return metrics;
        }


        MendelianViolationsByFamily violations() {
            return violations;
        }

        public static MendelianViolationDetector.Result merge(final Collection<MendelianViolationDetector.Result> results) {
            final Collection<Collection<MendelianViolationMetrics>> metricCollections = new ArrayList<>();
            final Collection<MendelianViolationsByFamily> violationCollections = new ArrayList<>();
            for (final MendelianViolationDetector.Result result : results) {
                metricCollections.add(result.metrics());
                violationCollections.add(result.violations());
            }

            return new MendelianViolationDetector.Result(mergeMetrics(metricCollections), mergeViolations(violationCollections));
        }

        private static MendelianViolationsByFamily mergeViolations(final Collection<MendelianViolationsByFamily> resultsToReduce) {
            final MendelianViolationsByFamily masterFamilyViolationsMap = new MendelianViolationsByFamily();
            
            for (final Map<String, Collection<VariantContext>> childFamilyViolationsMap : resultsToReduce) {
                for (final String childFamily : childFamilyViolationsMap.keySet()) {
                    masterFamilyViolationsMap.get(childFamily).addAll(childFamilyViolationsMap.get(childFamily));
                }
            }
            
            return masterFamilyViolationsMap;
        }

        /** Flattens out the provided metrics, collates them by "sample", and merges those collations. */
        private static  Collection<MendelianViolationMetrics> mergeMetrics(final Collection<Collection<MendelianViolationMetrics>> resultsToReduce) {
            final Collection<MendelianViolationMetrics> allMetrics = new ArrayList<>();
            resultsToReduce.forEach(allMetrics::addAll);

            final Map<String, List<MendelianViolationMetrics>> sampleToMetricsMap =
                    allMetrics
                            .stream()
                            .collect(Collectors
                                    .groupingBy(m -> String.format("%s|%s|%s|%s", m.FAMILY_ID, m.FATHER, m.MOTHER, m.OFFSPRING)));

            return sampleToMetricsMap
                    .values()
                    .stream()
                    .map(a-> (MendelianViolationMetrics) new MendelianViolationMetrics().merge(a))
                    .collect(Collectors.<MendelianViolationMetrics, List<MendelianViolationMetrics>>toCollection(ArrayList<MendelianViolationMetrics>::new));
        }
    }
}
