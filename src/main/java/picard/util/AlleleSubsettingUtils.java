package picard.util;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.*;
import picard.fingerprint.Snp;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Utilities class containing methods for restricting {@link VariantContext} and {@link GenotypesContext} objects to a
 * reduced set of alleles, as well as for choosing the best set of alleles to keep and for cleaning up annotations and
 * genotypes after subsetting.
 *
 * @author David Benjamin davidben@broadinstitute.org;
 * @author Yossi Farjoun farjoun@broadinstitute.org;
 */

public final class AlleleSubsettingUtils {

    static final Allele NON_REF_ALLELE = Allele.create("<NON_REF>");
    static List<Allele> DIPLOID_NO_CALL = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);

    /**
     * Method to subset the alleles in the VariantContext to those in the input snp.
     * Will also subset appropriately the AD and PL fields in all the genotypes, and in
     * the case of a monomorphic variant with a NON_REF allele, will replace the NON_REF
     * allele with the non-reference allele from snp.
     * <p>
     * returns null if it is not possible to subset ctx to the alleles in snp.
     *
     * @param ctx The VariantContext to subset.
     * @param snp Snp whose alleles are used for subsetting ctx.
     * @return a VariantContext subsetted to the alleles in snp, or null if unable.
     */
    public static VariantContext subsetVCToMatchSnp(final VariantContext ctx, final Snp snp) {
        if (ctx.isFiltered()) return null;

        // we can use this VC if it contains both alleles in snp as alleles, one of which is the reference,
        // or if it contains the NON_REF_ALLELE and one of the snp's alleles, and its reference allele is length 1.

        //TODO: figure out a way to identify the allele corresponding to the snp in the presence of a deletion

        if (ctx.getReference().length() != 1) return null;

        // one of the alleles in snp must match the reference allele
        final Optional<Byte> referenceAlleleMaybe = Stream.of(snp.getAllele1(), snp.getAllele2())
                .filter(b -> StringUtil.toUpperCase(b) == StringUtil.toUpperCase(ctx.getReference().getBases()[0]))
                .findAny();
        if (!referenceAlleleMaybe.isPresent()) return null;

        final byte refAllele = referenceAlleleMaybe.get();
        final byte otherAllele = snp.getAllele1() == refAllele ? snp.getAllele2() : snp.getAllele1();

        // do we have otherAllele in ctx?
        final Optional<Allele> altAlleleMaybe = ctx.getAlternateAlleles()
                .stream()
                .filter(a -> a.length() == 1 && StringUtil.toUpperCase(a.getBases()[0]) == StringUtil.toUpperCase(otherAllele))
                .findAny();

        if (altAlleleMaybe.isPresent()) {
            if (ctx.isBiallelic()) return ctx;

            return AlleleSubsettingUtils.subsetAlleles(ctx, Arrays.asList(ctx.getReference(), altAlleleMaybe.get()));
        }

        // if not, perhaps we have NON_REF_ALLELE
        final Optional<Allele> nonRefAlleleMaybe = ctx.getAlternateAlleles().stream().filter(a -> a.equals(NON_REF_ALLELE)).findAny();
        if (nonRefAlleleMaybe.isPresent()) {

            final VariantContext vcSubsetted = ctx.isBiallelic() ? ctx : AlleleSubsettingUtils.subsetAlleles(ctx, Arrays.asList(ctx.getReference(), nonRefAlleleMaybe.get()));
            return AlleleSubsettingUtils.swapAlleles(vcSubsetted, NON_REF_ALLELE, Allele.create(otherAllele));
        }
        return null;
    }

    public static VariantContext subsetAlleles(final VariantContext originalVc, final List<Allele> allelesToKeep) {
        VariantContextBuilder vcBuilder = new VariantContextBuilder(originalVc).alleles(allelesToKeep);

        GenotypesContext newGenotypes = subsetAlleles(originalVc.getGenotypes(), originalVc.getAlleles(), allelesToKeep);
        vcBuilder.genotypes(newGenotypes);
        return vcBuilder.make();
    }

    /**
     * Swaps one of the alleles in a VC (and its genotypes) with another.
     *
     * @param originalVc The {@link VariantContext} whose oldAllele will be swapped for newAllele.
     * @param oldAllele The {@link Allele} in originalVc that needs to be replaced.
     * @param newAllele The new {@link Allele} to use instead of oldAllele.
     *
     * @return A new {@link VariantContext} with newAllele swapped in for oldAllele.
     * @throws IllegalArgumentException if originalVc doesn't contain oldAllele.
     */
    public static VariantContext swapAlleles(final VariantContext originalVc, final Allele oldAllele, final Allele newAllele) throws IllegalArgumentException{
        if (!originalVc.getAlleles().contains(oldAllele)) throw new IllegalArgumentException("Couldn't find allele " +
                oldAllele + " in VariantContext " + originalVc);

        final List<Allele> alleles = new ArrayList<>(originalVc.getAlleles());
        alleles.set(alleles.indexOf(oldAllele), newAllele);

        VariantContextBuilder vcBuilder = new VariantContextBuilder(originalVc).alleles(alleles);
        GenotypesContext newGTs = GenotypesContext.create(originalVc.getGenotypes().size());

        for (final Genotype g : originalVc.getGenotypes()) {
            if (!g.getAlleles().contains(oldAllele)) {
                newGTs.add(g);
            } else {
                final GenotypeBuilder gb = new GenotypeBuilder(g);
                gb.alleles(g.getAlleles().stream().map(a -> a.equals(oldAllele) ? newAllele : a).collect(Collectors.toList()));
                newGTs.add(gb.make());
            }
        }
        vcBuilder.genotypes(newGTs);
        return vcBuilder.make();
    }

    /**
     * Create the new GenotypesContext with the subsetted PLs and ADs
     * Expects allelesToKeep to be in the same order
     * in which they are in originalAlleles.
     * <p>
     *
     * @param originalGs      the original GenotypesContext
     * @param originalAlleles the original alleles
     * @param allelesToKeep   the subset of alleles to use with the new Genotypes
     * @return a new non-null GenotypesContext
     */
    public static GenotypesContext subsetAlleles(final GenotypesContext originalGs,
                                                 final List<Allele> originalAlleles,
                                                 final List<Allele> allelesToKeep) {
        nonNull(originalGs, "original GenotypesContext must not be null.");
        nonNull(allelesToKeep, "allelesToKeep is null.");
        nonEmpty(allelesToKeep, "must keep at least one allele.");
        validateTrue(allelesToKeep.get(0).isReference(), "First allele must be the reference allele.");
        validateTrue(allelesToKeep.stream().allMatch(originalAlleles::contains), "OriginalAlleles must contain allelesToKeep.");

        int indexOfLast = -1;
        for (Allele a : allelesToKeep) {
            validateTrue(indexOfLast < originalAlleles.indexOf(a), "alleles to keep must maintain the order of the original alleles.");
            indexOfLast = originalAlleles.indexOf(a);
        }

        final int[] allelesIndex = allelesToKeep.stream().mapToInt(originalAlleles::indexOf).toArray();
        final GenotypesContext newGTs = GenotypesContext.create(originalGs.size());
        int[] subsettedLikelihoodIndices = subsettedPLIndices(originalAlleles, allelesToKeep);

        for (final Genotype g : originalGs) {
            validateTrue(g.getPloidy() == 2, "only implemented for ploidy 2 for now.");

            final int expectedNumLikelihoods = GenotypeLikelihoods.numLikelihoods(allelesToKeep.size(), 2);
            // create the new likelihoods array from the alleles we are to use
            int[] newPLs = null;
            double newLog10GQ = -1;
            if (g.hasLikelihoods()) {

                int[] originalPLs = g.getPL();

                if (originalPLs.length != expectedNumLikelihoods) {
                    newPLs = Arrays.stream(subsettedLikelihoodIndices)
                            .map(idx -> originalPLs[idx]).toArray();
                    final int minLikelihood = MathUtil.min(newPLs);
                    for (int i = 0; i < expectedNumLikelihoods; i++) {
                        newPLs[i] = newPLs[i] - minLikelihood;
                    }
                    final int indexOfMostLikely = MathUtil.indexOfMin(newPLs);

                    newLog10GQ = GenotypeLikelihoods.getGQLog10FromLikelihoods(indexOfMostLikely, GenotypeLikelihoods.fromPLs(newPLs).getAsVector());
                } else {
                    newPLs = null;
                }
            }

            final GenotypeBuilder gb;
            if (newPLs == null) {
                gb = new GenotypeBuilder(g).noPL().noGQ().alleles(DIPLOID_NO_CALL);
            } else {
                gb = new GenotypeBuilder(g).PL(newPLs).log10PError(newLog10GQ);
                final List<Integer> originalDiploidAlleles = GenotypeLikelihoods.getAlleles(MathUtil.indexOfMin(newPLs), 2);
                gb.alleles(originalDiploidAlleles.stream().map(allelesToKeep::get).collect(Collectors.toList()));
            }

            // restrict AD to the new allele subset
            if (g.hasAD()) {
                final int[] oldAD = g.getAD();
                final int[] newAD = IntStream.range(0, allelesToKeep.size()).map(n -> oldAD[allelesIndex[n]]).toArray();
                gb.AD(newAD);
            }
            newGTs.add(gb.make());
        }
        return newGTs;
    }

    /**
     * Given a list of original alleles and a subset of new alleles to retain, find the array of old PL indices that correspond
     * to new PL indices i.e. result[7] = old PL index of genotype containing same alleles as the new genotype with PL index 7.
     * <p>
     * This method is written in terms f indices rather than subsetting PLs directly in order to produce output that can be
     * recycled from sample to sample, provided that the ploidy is the same.
     *
     * @param originalAlleles List of original alleles
     * @param newAlleles      New alleles -- must be a subset of {@code originalAlleles}
     * @return old PL indices of new genotypes
     */
    public static int[] subsettedPLIndices(final List<Allele> originalAlleles, final List<Allele> newAlleles) {
        final int[] result = new int[GenotypeLikelihoods.numLikelihoods(newAlleles.size(), 2)];

        for (int oldPLIndex = 0; oldPLIndex < GenotypeLikelihoods.numLikelihoods(originalAlleles.size(), 2); oldPLIndex++) {
            final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair allelePairFromPLIndex = GenotypeLikelihoods.getAllelePair(oldPLIndex);

            final Allele allele1 = originalAlleles.get(allelePairFromPLIndex.alleleIndex1);
            final Allele allele2 = originalAlleles.get(allelePairFromPLIndex.alleleIndex2);

            final boolean containsOnlyNewAlleles = newAlleles.contains(allele1) && newAlleles.contains(allele2);

            if (containsOnlyNewAlleles) {
                // thus we want this PL...but we need to figure out where its new index is....
                final int newPLIndex = GenotypeLikelihoods.calculatePLindex(newAlleles.indexOf(allele1), newAlleles.indexOf(allele2));
                result[newPLIndex] = oldPLIndex;
            }
        }
        return result;
    }

    /**
     * Checks that a boolean is true and returns the same object or throws an {@link IllegalArgumentException}
     *
     * @param condition any Object
     * @param msg       the text message that would be passed to the exception thrown when {@code !condition}.
     * @throws IllegalArgumentException if a {@code !condition}
     */
    private static void validateTrue(final boolean condition, final String msg) {
        if (!condition) {
            throw new IllegalArgumentException(msg);
        }
    }

    /**
     * Checks that an {@link Object} is not {@code null} and returns the same object or throws an {@link IllegalArgumentException}
     *
     * @param object  any Object
     * @param message the text message that would be passed to the exception thrown when {@code o == null}.
     * @return the same object
     * @throws IllegalArgumentException if a {@code o == null}
     */
    private static <T> T nonNull(final T object, String message) {
        if (object == null) {
            throw new IllegalArgumentException(message);
        }
        return object;
    }

    /**
     * Checks that a {@link Collection} is not {@code null} and that it is not empty.
     * If it's non-null and non-empty it returns the input, otherwise it throws an {@link IllegalArgumentException}
     *
     * @param collection any Collection
     * @param message    a message to include in the output
     * @return the original collection
     * @throws IllegalArgumentException if collection is null or empty
     */
    private static <I, T extends Collection<I>> T nonEmpty(T collection, String message) {
        nonNull(collection, "The collection is null: " + message);
        if (collection.isEmpty()) {
            throw new IllegalArgumentException("The collection is empty: " + message);
        } else {
            return collection;
        }
    }
}
