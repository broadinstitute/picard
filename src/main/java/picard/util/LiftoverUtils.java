package picard.util;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.*;
import picard.vcf.LiftoverVcf;

import java.util.*;
import java.util.stream.Collectors;


public class LiftoverUtils {
    /**
     * This will take an input VariantContext and lift to the provided interval.
     * If this interval is in the opposite orientation, all alleles and genotypes will be reverse complemented
     * and indels will be left-aligned.  Currently this is only able to invert biallelic indels, and null will be
     * returned for any unsupported VC.
     * @param source The VariantContext to lift
     * @param target The target interval
     * @param refSeq The reference sequence, which should match the target interval
     * @param writeOriginalPosition If true, INFO field annotations will be added to store the original position and contig
     * @return The lifted VariantContext.  This will be null if the input VariantContext could not be lifted.
     */
    public static VariantContext liftVariant(VariantContext source, Interval target, ReferenceSequence refSeq, boolean writeOriginalPosition){
        if (target == null) {
            return null;
        }

        final VariantContextBuilder builder;
        if (target.isNegativeStrand()){
            builder = reverseComplementVariantContext(source, target, refSeq);
        }
        else {
            builder = liftSimpleVariantContext(source, target);
        }

        if (builder == null){
            return null;
        }

        builder.filters(source.getFilters());
        builder.log10PError(source.getLog10PError());
        builder.attributes(source.getAttributes());
        builder.id(source.getID());

        if (writeOriginalPosition) {
            builder.attribute(LiftoverVcf.ORIGINAL_CONTIG, source.getContig());
            builder.attribute(LiftoverVcf.ORIGINAL_START, source.getStart());
        }

        return builder.make();
    }

    protected static VariantContextBuilder liftSimpleVariantContext(VariantContext source, Interval target){
        if (target == null || source.getReference().length() != target.length()) {
            return null;
        }

        // Build the new variant context
        final VariantContextBuilder builder = new VariantContextBuilder(source);
        builder.chr(target.getContig());
        builder.start(target.getStart());
        builder.stop(target.getEnd());

        return builder;
    }

    protected static VariantContextBuilder reverseComplementVariantContext(VariantContext source, Interval target, ReferenceSequence refSeq){
        if (target.isPositiveStrand()){
            throw new IllegalArgumentException("This should only be called for negative strand liftovers");
        }

        //not currently supported
        if (source.isIndel() && !source.isBiallelic()){
            return null;
        }

        List<Allele> origAlleles = new ArrayList<>(source.getAlleles());
        VariantContextBuilder vcb = new VariantContextBuilder(source);
        vcb.chr(target.getContig());

        final boolean addToStart = source.isIndel() && target.getStart() > 1;

        final int start = target.getStart() - (addToStart ? 1 : 0);
        vcb.start(start);

        final int stop = target.getEnd() - (addToStart ? 1 : 0);
        vcb.stop(stop);

        vcb.alleles(reverseComplementAlleles(origAlleles, target, refSeq, source.isIndel(), addToStart));
        if (source.isIndel()){
            leftAlignVariant(vcb, start, stop, vcb.getAlleles(), refSeq);
        }

        vcb.genotypes(fixGenotypes(source.getGenotypes(), origAlleles, vcb.getAlleles()));

        return vcb;
    }

    private static List<Allele> reverseComplementAlleles(List<Allele> originalAlleles, Interval target, ReferenceSequence refSeq, boolean isBiAllelicIndel, boolean addToStart){
        final List<Allele> alleles = new ArrayList<>();

        for (final Allele oldAllele : originalAlleles) {
            alleles.add(LiftoverUtils.reverseComplement(oldAllele, target, refSeq, isBiAllelicIndel, addToStart));
        }

        return alleles;
    }

    private static Allele reverseComplement(Allele oldAllele, Interval target, ReferenceSequence referenceSequence, boolean isBiAllelicIndel, boolean addToStart){

        if (oldAllele.isSymbolic()){
            return oldAllele;
        }
        else if (isBiAllelicIndel){
            // target.getStart is 1-based, reference bases are 0-based
            final StringBuilder alleleBuilder = new StringBuilder(target.getEnd() - target.getStart() + 1);

            if (addToStart) alleleBuilder.append((char) referenceSequence.getBases()[target.getStart() - 2]);
            alleleBuilder.append(SequenceUtil.reverseComplement(oldAllele.getBaseString().substring(1, oldAllele.length())));
            if (!addToStart) alleleBuilder.append((char) referenceSequence.getBases()[target.getEnd() - 1]);

            return Allele.create(alleleBuilder.toString(), oldAllele.isReference());
        }
        else {
            return Allele.create(SequenceUtil.reverseComplement(oldAllele.getBaseString()), oldAllele.isReference());
        }
    }

    protected static GenotypesContext fixGenotypes(final GenotypesContext originals, List<Allele> originalAlleles, List<Allele> newAlleles) {
        // optimization: if nothing needs to be fixed then don't bother
        if (originalAlleles.equals(newAlleles)) {
            return originals;
        }

        final GenotypesContext fixedGenotypes = GenotypesContext.create(originals.size());
        for (final Genotype genotype : originals) {
            final List<Allele> fixedAlleles = new ArrayList<>();
            for (final Allele allele : genotype.getAlleles()) {
                if (allele.isSymbolic() || allele.isNoCall()){
                    fixedAlleles.add(allele);
                }
                else {
                    int idx = originalAlleles.indexOf(allele);
                    if (idx == -1) {
                        throw new IllegalStateException("Allele not found: " + allele.toString() + ", " + originalAlleles + "/ " + newAlleles);
                    }
                    fixedAlleles.add(newAlleles.get(idx));
                }
            }
            fixedGenotypes.add(new GenotypeBuilder(genotype).alleles(fixedAlleles).make());
        }
        return fixedGenotypes;
    }

    /**
     *    Normalizes and left aligns a {@link VariantContext}.
     *
     *    Based on Adrian Tan, Gonçalo R. Abecasis and Hyun Min Kang. (2015)
     *    Unified Representation of Genetic Variants. Bioinformatics.
     *
     * @return a new {@link VariantContext} which represents the same variation as vc but has
     * been normalized and left-aligned
     */
    protected static void leftAlignVariant(VariantContextBuilder builder, int start, int end, final List<Allele> alleles, final ReferenceSequence referenceSequence) {

        boolean changesInAlleles = true;

        final Map<Allele, byte[]> alleleBasesMap = new HashMap<>();
        alleles.forEach(a -> alleleBasesMap.put(a, a.getBases()));

        // 1. while changes in alleles do
        while (changesInAlleles) {

            changesInAlleles = false;
            // 2. if alleles end with the same nucleotide then
            if (alleleBasesMap.values().stream()
                    .collect(Collectors.groupingBy(a -> a[a.length - 1], Collectors.toSet()))
                    .size() == 1 && end > 1) {
                // 3. truncate rightmost nucleotide of each allele
                for (final Allele allele : alleleBasesMap.keySet()) {
                    alleleBasesMap.put(allele, truncateBase(alleleBasesMap.get(allele), true));
                }
                changesInAlleles = true;
                end--;
                // 4. end if
            }

            // 5. if there exists an empty allele then
            if (alleleBasesMap.values().stream()
                    .map(a -> a.length)
                    .anyMatch(l -> l == 0)) {
                // 6. extend alleles 1 nucleotide to the left
                for (final Allele allele : alleleBasesMap.keySet()) {
                    // the first -1 for zero-base (getBases) versus 1-based (variant position)
                    // another   -1 to get the base prior to the location of the start of the allele
                    final byte extraBase =  (start > 1) ?
                            referenceSequence.getBases()[start - 2] :
                            referenceSequence.getBases()[end];

                    alleleBasesMap.put(allele, extendOneBase(alleleBasesMap.get(allele), extraBase));
                }
                changesInAlleles = true;
                start--;

                // 7. end if
            }
        }

        // 8. while leftmost nucleotide of each allele are the same and all alleles have length 2 or more do
        while (alleleBasesMap.values().stream()
                .allMatch(a -> a.length >= 2) &&

                alleleBasesMap.values().stream()
                        .collect(Collectors.groupingBy(a -> a[0], Collectors.toSet()))
                        .size() == 1
                ) {

            //9. truncate the leftmost base of the alleles
            for (final Allele allele : alleleBasesMap.keySet()) {
                alleleBasesMap.put(allele, truncateBase(alleleBasesMap.get(allele), false));
            }
            start++;
        }

        builder.start(start);
        builder.stop(end);

        final Map<Allele, Allele> fixedAlleleMap = alleleBasesMap.entrySet().stream()
                .collect(Collectors.toMap(Map.Entry::getKey, me -> Allele.create(me.getValue(), me.getKey().isReference())));

        //retain original order:
        List<Allele> fixedAlleles = alleles.stream().map(a -> fixedAlleleMap.get(a)).collect(Collectors.toList());
        //Collections.reverse(fixedAlleles);
        builder.alleles(fixedAlleles);
    }

    private static byte[] truncateBase(final byte[] allele, final boolean truncateRightmost) {
        return Arrays.copyOfRange(allele, truncateRightmost ? 0 : 1, truncateRightmost ?
                allele.length - 1 :
                allele.length);
    }

    //creates a new byte array with the base added at the begining
    private static byte[] extendOneBase(final byte[] bases, final byte base) {

        final byte[] newBases = new byte[bases.length + 1];

        System.arraycopy(bases, 0, newBases, 1, bases.length);
        newBases[0] = base;

        return newBases;
    }
}
