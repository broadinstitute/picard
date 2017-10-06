package picard.util;

import htsjdk.variant.variantcontext.*;

import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
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
    private AlleleSubsettingUtils() {
    }  // prevent instantiation

    public static VariantContext subsetAlleles(final VariantContext originalVc, final List<Allele> allelesToKeep) {
        VariantContextBuilder vcBuilder = new VariantContextBuilder(originalVc).alleles(allelesToKeep);

        GenotypesContext newGenotypes = subsetAlleles(originalVc.getGenotypes(), originalVc.getAlleles(), allelesToKeep);
        vcBuilder.genotypes(newGenotypes);
        return vcBuilder.make();
    }

    /**
     * swaps one of the alleles in a VC (and its genotypes) with another.
     */
    public static VariantContext swapAlleles(final VariantContext originalVc, final Allele oldAllele, final Allele newAllele) {
        if (!originalVc.getAlleles().contains(oldAllele)) return originalVc;

        final LinkedList<Allele> alleles = new LinkedList<>(originalVc.getAlleles());
        alleles.set(alleles.indexOf(oldAllele), newAllele);

        VariantContextBuilder vcBuilder = new VariantContextBuilder(originalVc).alleles(alleles);
        GenotypesContext newGTs = GenotypesContext.create(originalVc.getGenotypes().size());

        for (final Genotype g : originalVc.getGenotypes()) {
            if (!g.getAlleles().contains(oldAllele)) {
                newGTs.add(g);
            } else {
                final GenotypeBuilder gb = new GenotypeBuilder(g);
                final LinkedList<Allele> genotypeAlleles = new LinkedList<>(g.getAlleles());
                genotypeAlleles.set(genotypeAlleles.indexOf(oldAllele), newAllele);
                gb.alleles(genotypeAlleles);
                newGTs.add(gb.make());
            }
        }
        vcBuilder.genotypes(newGTs);
        return vcBuilder.make();
    }

    /**
     * Create the new GenotypesContext with the subsetted PLs and ADs
     * <p>
     * Will reorder subsetted alleles according to the ordering provided by the list allelesToKeep
     *
     * @param originalGs      the original GenotypesContext
     * @param originalAlleles the original alleles
     * @param allelesToKeep   the subset of alleles to use with the new Genotypes
     * @return a new non-null GenotypesContext
     */
    public static GenotypesContext subsetAlleles(final GenotypesContext originalGs,
                                                 final List<Allele> originalAlleles,
                                                 final List<Allele> allelesToKeep) {
        nonNull(originalGs, "original GenotypesContext must not be null");
        nonNull(allelesToKeep, "allelesToKeep is null");
        nonEmpty(allelesToKeep, "must keep at least one allele");
        validateTrue(allelesToKeep.get(0).isReference(), "First allele must be the reference allele");
        validateTrue(allelesToKeep.stream().allMatch(originalAlleles::contains), "OriginalAlleles must contain allelesToKeep");

        int indexOfLast = -1;
        for (Allele a : allelesToKeep) {
            validateTrue(indexOfLast < originalAlleles.indexOf(a), "alleles to keep must maintain the order of the original alleles");
            indexOfLast = originalAlleles.indexOf(a);
        }

        final int allelesIndex[] = allelesToKeep.stream().mapToInt(originalAlleles::indexOf).toArray();
        final GenotypesContext newGTs = GenotypesContext.create(originalGs.size());
        int[] subsettedLikelihoodIndices = subsettedPLIndices(originalAlleles, allelesToKeep);

        for (final Genotype g : originalGs) {
            validateTrue(g.getPloidy() == 2, "only implemented for ploidy 2 for now");

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
                    final int PLindex = MathUtil.indexOfMin(newPLs);

                    newLog10GQ = GenotypeLikelihoods.getGQLog10FromLikelihoods(PLindex, GenotypeLikelihoods.fromPLs(newPLs).getAsVector());
                } else {
                    newPLs = null;
                }
            }

            final boolean useNewLikelihoods = newPLs != null;
            final GenotypeBuilder gb = useNewLikelihoods ? new GenotypeBuilder(g).PL(newPLs).log10PError(newLog10GQ) : new GenotypeBuilder(g).noPL().noGQ();
            GenotypeLikelihoods.GenotypeLikelihoodsAllelePair allelePair = GenotypeLikelihoods.getAllelePair(MathUtil.indexOfMin(newPLs));
            gb.alleles(Stream.of(allelePair.alleleIndex1, allelePair.alleleIndex2).map(allelesToKeep::get).collect(Collectors.toList()));
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
     * Checks that a boolea is true and returns the same object or throws an {@link IllegalArgumentException}
     *
     * @param condition  any Object
     * @param msg the text message that would be passed to the exception thrown when {@code !condition}.
     *
     * @throws IllegalArgumentException if a {@code !condition}
     */
    public static void validateTrue(final boolean condition, final String msg) {
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
    public static <T> T nonNull(final T object, String message) {
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
    public static <I, T extends Collection<I>> T nonEmpty(T collection, String message) {
        nonNull(collection, "The collection is null: " + message);
        if (collection.isEmpty()) {
            throw new IllegalArgumentException("The collection is empty: " + message);
        } else {
            return collection;
        }
    }
}
