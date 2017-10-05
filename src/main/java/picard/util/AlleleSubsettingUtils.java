package picard.util;

import htsjdk.variant.variantcontext.*;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Utilities class containing methods for restricting {@link VariantContext} and {@link GenotypesContext} objects to a
 * reduced set of alleles, as well as for choosing the best set of alleles to keep and for cleaning up annotations and
 * genotypes after subsetting.
 *
 * @author David Benjamin davidben@broadinstitute.org;
 */

public final class AlleleSubsettingUtils {
    private AlleleSubsettingUtils() {}  // prevent instantiation

    /**
     * Create the new GenotypesContext with the subsetted PLs and ADs
     *
     * Will reorder subsetted alleles according to the ordering provided by the list allelesToKeep
     *
     * @param originalGs               the original GenotypesContext
     * @param originalAlleles          the original alleles
     * @param allelesToKeep            the subset of alleles to use with the new Genotypes
     * @param depth                    the original variant DP or 0 if there was no DP
     * @return                         a new non-null GenotypesContext
     */
    public static GenotypesContext subsetAlleles(final GenotypesContext originalGs,
                                                 final List<Allele> originalAlleles,
                                                 final List<Allele> allelesToKeep,
                                                 final int depth) {
        nonNull(originalGs, "original GenotypesContext must not be null");
        nonNull(allelesToKeep, "allelesToKeep is null");
        nonEmpty(allelesToKeep, "must keep at least one allele");
        validateTrue(allelesToKeep.get(0).isReference(), "First allele must be the reference allele");
        validateTrue(allelesToKeep.stream().allMatch(originalAlleles::contains),"OriginalAlleles must contain allelesToKeep");

        final int allelesIndex[] = allelesToKeep.stream().mapToInt(originalAlleles::indexOf).toArray();

        final GenotypesContext newGTs = GenotypesContext.create(originalGs.size());

        int[] subsettedLikelihoodIndices = subsettedPLIndices(originalAlleles, allelesToKeep);
        for (final Genotype g : originalGs) {
            validateTrue(g.getPloidy()==2,"only implemented for ploidy 2 for now");

            final int expectedNumLikelihoods = GenotypeLikelihoods.numLikelihoods(originalAlleles.size(), 2);
            // create the new likelihoods array from the alleles we are allowed to use
            double[] newLikelihoods = null;
            double newLog10GQ = -1;
            if (g.hasLikelihoods()) {

                double[] originalLikelihoods = g.getLikelihoods().getAsVector();

                if (originalLikelihoods.length != expectedNumLikelihoods) {
                    newLikelihoods = Arrays.stream(subsettedLikelihoodIndices)
                            .mapToDouble(idx -> originalLikelihoods[idx]).toArray();
                    final double minPl = MathUtil.min(newLikelihoods);
                    for (int i = 0; i < expectedNumLikelihoods; i++) {
                        newLikelihoods[i] = newLikelihoods[i] - minPl;
                    }
                    final int PLindex = MathUtil.indexOfMax(newLikelihoods);
                    newLog10GQ = GenotypeLikelihoods.getGQLog10FromLikelihoods(PLindex, newLikelihoods);
                } else {
                    newLikelihoods = null;
                }
            }

            final boolean useNewLikelihoods = newLikelihoods != null;
            final GenotypeBuilder gb = useNewLikelihoods ? new GenotypeBuilder(g).PL(newLikelihoods).log10PError(newLog10GQ) : new GenotypeBuilder(g).noPL().noGQ();
            GenotypeLikelihoods.GenotypeLikelihoodsAllelePair allelePair = GenotypeLikelihoods.getAllelePair(MathUtil.indexOfMin(newLikelihoods));
            gb.alleles(Stream.of(allelePair.alleleIndex1, allelePair.alleleIndex2).map(allelesToKeep::get).collect(Collectors.toList()));
            // restrict AD to the new allele subset
            if(g.hasAD()) {
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
     *
     * This method is written in terms f indices rather than subsetting PLs directly in order to produce output that can be
     * recycled from sample to sample, provided that the ploidy is the same.
     *
     * @param originalAlleles       List of original alleles
     * @param newAlleles            New alleles -- must be a subset of {@code originalAlleles}
     * @return                      old PL indices of new genotypes
     */
    public static int[] subsettedPLIndices(final List<Allele> originalAlleles, final List<Allele> newAlleles) {
        final int[] result = new int[GenotypeLikelihoods.numLikelihoods(newAlleles.size(), 2)];

        for (int oldPLIndex = 0; oldPLIndex < GenotypeLikelihoods.numLikelihoods(originalAlleles.size(),2); oldPLIndex++) {
            final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair allelePairFromPLIndex = GenotypeLikelihoods.getAllelePair(oldPLIndex);

            final Allele allele1 = originalAlleles.get(allelePairFromPLIndex.alleleIndex1);
            final Allele allele2 = originalAlleles.get(allelePairFromPLIndex.alleleIndex2);

            final boolean containsOnlyNewAlleles = newAlleles.contains(allele1) && newAlleles.contains(allele2);

            if (containsOnlyNewAlleles) {
                // thus we want this PL...but we need to figure out where it's new index is....

                final int newPLIndex = GenotypeLikelihoods.calculatePLindex(newAlleles.indexOf(allele1), newAlleles.indexOf(allele2));
                result[newPLIndex] = oldPLIndex;
            }
        }
        return  result;
    }


    public static void validateTrue(final boolean condition, final String msg){
        if (!condition){
            throw new IllegalArgumentException(msg);
        }
    }

    /**
     * Checks that an {@link Object} is not {@code null} and returns the same object or throws an {@link IllegalArgumentException}
     * @param object any Object
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
     * @param collection any Collection
     * @param message a message to include in the output
     * @return the original collection
     * @throws IllegalArgumentException if collection is null or empty
     */
    public static <I, T extends Collection<I>> T nonEmpty(T collection, String message){
        nonNull(collection, "The collection is null: " + message);
        if(collection.isEmpty()){
            throw new IllegalArgumentException("The collection is empty: " + message);
        } else {
            return collection;
        }
    }
}
