package picard.util;

/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.ArrayUtils;
import picard.vcf.LiftoverVcf;

import java.util.*;
import java.util.stream.Collectors;

public class LiftoverUtils {

    /**
     * Attribute used to store the fact that the alt and ref alleles of the variant have been swapped, while all the INFO annotations have not.
     */
    public static final String SWAPPED_ALLELES = "SwappedAlleles";
    public static final String REV_COMPED_ALLELES = "ReverseComplementedAlleles";

    /**
     * Default list of attributes that need to be reversed or dropped from the INFO field when alleles have been swapped.
     */
    public static final Collection<String> DEFAULT_TAGS_TO_REVERSE = Arrays.asList("AF");
    public static final Collection<String> DEFAULT_TAGS_TO_DROP = Arrays.asList("MAX_AF");

    public static final Log log = Log.getInstance(LiftoverUtils.class);

    /**
     * This will take an input VariantContext and lift to the provided interval.
     * If this interval is in the opposite orientation, all alleles and genotypes will be reverse complemented
     * and indels will be left-aligned.  Currently this is only able to invert biallelic indels, and null will be
     * returned for any unsupported VC.
     *
     * @param source                The VariantContext to lift
     * @param target                The target interval
     * @param refSeq                The reference sequence, which should match the target interval
     * @param writeOriginalPosition If true, INFO field annotations will be added to store the original position and contig
     * @param writeOriginalAlleles  If true, an INFO field annotation will be added to store the original alleles in order.  This can be useful in the case of more complex liftovers, like reverse complements, left aligned indels or swapped REF/ALT
     * @return The lifted VariantContext.  This will be null if the input VariantContext could not be lifted.
     */
    public static VariantContext liftVariant(final VariantContext source, final Interval target, final ReferenceSequence refSeq, final boolean writeOriginalPosition, final boolean writeOriginalAlleles) {
        if (target == null) {
            return null;
        }

        final VariantContextBuilder builder;
        if (target.isNegativeStrand()) {
            builder = reverseComplementVariantContext(source, target, refSeq);
        } else {
            builder = liftSimpleVariantContext(source, target);
        }

        if (builder == null) {
            return null;
        }

        builder.filters(source.getFilters());
        builder.log10PError(source.getLog10PError());

        // If any of the source alleles is symbolic, do not populate the END tag (protecting things like <NON_REF> from
        // getting screwed up in reverse-complemented variants
        if (source.hasAttribute(VCFConstants.END_KEY) && builder.getAlleles().stream().noneMatch(Allele::isSymbolic)) {
            // TODO: add start() and stop() methods to the builder in htsjdk and use stop() here.
            builder.attribute(VCFConstants.END_KEY, builder.make().getEnd());
        } else {
            builder.rmAttribute(VCFConstants.END_KEY);
        }

        // make sure that the variant isn't mistakenly set as "SwappedAlleles"
        builder.rmAttribute(SWAPPED_ALLELES);
        if (target.isNegativeStrand()) {
            builder.attribute(REV_COMPED_ALLELES, true);
        } else {
            // make sure that the variant isn't mistakenly set as "ReverseComplementedAlleles" (from a previous liftover, say)
            builder.rmAttribute(REV_COMPED_ALLELES);
        }

        builder.id(source.getID());

        if (writeOriginalPosition) {
            builder.attribute(LiftoverVcf.ORIGINAL_CONTIG, source.getContig());
            builder.attribute(LiftoverVcf.ORIGINAL_START, source.getStart());
        }

        if (writeOriginalAlleles && !source.getAlleles().equals(builder.getAlleles())) {
            builder.attribute(LiftoverVcf.ORIGINAL_ALLELES, allelesToStringList(source.getAlleles()));
        }

        return builder.make();
    }

    /**
     * This is a utility method that will convert a list of alleles into a list of base strings.  Reference status
     * is ignored when creating these strings (i.e. 'A', not 'A*').  These strings should be sufficient
     * to recreate an Allele using Allele.create()
     *
     * @param alleles The list of alleles
     * @return A list of strings representing the bases of the input alleles.
     */
    protected static List<String> allelesToStringList(final List<Allele> alleles) {
        final List<String> ret = new ArrayList<>();
        alleles.forEach(a -> ret.add(a.isNoCall() ? Allele.NO_CALL_STRING : a.getDisplayString()));
        return ret;
    }

    protected static VariantContextBuilder liftSimpleVariantContext(VariantContext source, Interval target) {
        if (target == null || source.getReference().length() != target.length()) {
            return null;
        }

        // Build the new variant context.  Note: this will copy genotypes, annotations and filters
        final VariantContextBuilder builder = new VariantContextBuilder(source);
        builder.chr(target.getContig());
        builder.start(target.getStart());
        builder.stop(target.getEnd());

        return builder;
    }

    protected static VariantContextBuilder reverseComplementVariantContext(final VariantContext source, final Interval target, final ReferenceSequence refSeq) {
        if (target.isPositiveStrand()) {
            throw new IllegalArgumentException("This should only be called for negative strand liftovers");
        }

        final List<Allele> origAlleles = new ArrayList<>(source.getAlleles());
        final VariantContextBuilder vcb = new VariantContextBuilder(source);

        vcb.rmAttribute(VCFConstants.END_KEY);
        vcb.chr(target.getContig());

        // By convention, indels are left aligned and include the base prior to that indel.
        // When reverse complementing, will be necessary to include this additional base.
        // This check prevents the extremely rare situation in which the indel occurs on the
        // first base of the sequence
        final boolean addToStart = source.isIndel() && target.getStart() > 1;

        final int start = target.getStart() - (addToStart ? 1 : 0);
        vcb.start(start);

        final int stop = target.getEnd() - (addToStart ? 1 : 0);
        vcb.stop(stop);

        vcb.alleles(reverseComplementAlleles(origAlleles, target, refSeq, source.isIndel(), addToStart));

        if (!source.isSNP()) {
            // check that the reverse complemented bases match the new reference
            if (!referenceAlleleMatchesReferenceForIndel(vcb.getAlleles(), refSeq, start, stop)) {
                return null;
            }
            leftAlignVariant(vcb, start, stop, vcb.getAlleles(), refSeq);
        }

        vcb.genotypes(fixGenotypes(source.getGenotypes(), origAlleles, vcb.getAlleles()));

        return vcb;
    }

    private static List<Allele> reverseComplementAlleles(final List<Allele> originalAlleles, final Interval target, final ReferenceSequence refSeq, final boolean isIndel, final boolean addToStart) {
        final List<Allele> alleles = new ArrayList<>();

        for (final Allele oldAllele : originalAlleles) {
            alleles.add(LiftoverUtils.reverseComplement(oldAllele, target, refSeq, isIndel, addToStart));
        }

        return alleles;
    }

    private static Allele reverseComplement(final Allele oldAllele, final Interval target, final ReferenceSequence referenceSequence, final boolean isIndel, final boolean addToStart) {

        if (oldAllele.isSymbolic() || oldAllele.isNoCall() || oldAllele.equals(Allele.SPAN_DEL)) {
            return oldAllele;
        } else if (isIndel) {
            // target.getStart is 1-based, reference bases are 0-based
            final StringBuilder alleleBuilder = new StringBuilder(target.getEnd() - target.getStart() + 1);

            if (addToStart) {
                alleleBuilder.append((char) referenceSequence.getBases()[target.getStart() - 2]);
            }
            alleleBuilder.append(SequenceUtil.reverseComplement(oldAllele.getBaseString().substring(1, oldAllele.length())));
            if (!addToStart) {
                alleleBuilder.append((char) referenceSequence.getBases()[target.getEnd() - 1]);
            }

            return Allele.create(alleleBuilder.toString(), oldAllele.isReference());
        } else {
            return Allele.create(SequenceUtil.reverseComplement(oldAllele.getBaseString()), oldAllele.isReference());
        }
    }

    protected static GenotypesContext fixGenotypes(final GenotypesContext originals, final List<Allele> originalAlleles, final List<Allele> newAlleles) {
        // optimization: if nothing needs to be fixed then don't bother
        if (originalAlleles.equals(newAlleles)) {
            return originals;
        }

        if (originalAlleles.size() != newAlleles.size()) {
            throw new IllegalStateException("Error in allele lists: the original and new allele lists are not the same length: " + originalAlleles.toString() + " / " + newAlleles.toString());
        }

        final Map<Allele, Allele> alleleMap = new HashMap<>();
        for (int idx = 0; idx < originalAlleles.size(); idx++) {
            alleleMap.put(originalAlleles.get(idx), newAlleles.get(idx));
        }

        final GenotypesContext fixedGenotypes = GenotypesContext.create(originals.size());
        for (final Genotype genotype : originals) {
            final List<Allele> fixedAlleles = new ArrayList<>();
            for (final Allele allele : genotype.getAlleles()) {
                if (allele.isNoCall()) {
                    fixedAlleles.add(allele);
                } else {
                    Allele newAllele = alleleMap.get(allele);
                    if (newAllele == null) {
                        throw new IllegalStateException("Allele not found: " + allele.toString() + ", " + originalAlleles + "/ " + newAlleles);
                    }
                    fixedAlleles.add(newAllele);
                }
            }
            fixedGenotypes.add(new GenotypeBuilder(genotype).alleles(fixedAlleles).make());
        }
        return fixedGenotypes;
    }

    /**
     * method to swap the reference and alt alleles of a bi-allelic, SNP
     *
     * @param vc                   the {@link VariantContext} (bi-allelic SNP) that needs to have it's REF and ALT alleles swapped.
     * @param annotationsToReverse INFO field annotations (of double value) that will be reversed (x->1-x)
     * @param annotationsToDrop    INFO field annotations that will be dropped from the result since they are invalid when REF and ALT are swapped
     * @return a new {@link VariantContext} with alleles swapped, INFO fields modified and in the genotypes, GT, AD and PL corrected appropriately
     */
    public static VariantContext swapRefAlt(final VariantContext vc, final Collection<String> annotationsToReverse, final Collection<String> annotationsToDrop) {

        if (!vc.isBiallelic() || !vc.isSNP()) {
            throw new IllegalArgumentException("swapRefAlt can only process biallelic, SNPS, found " + vc.toString());
        }

        final VariantContextBuilder swappedBuilder = new VariantContextBuilder(vc);

        swappedBuilder.attribute(SWAPPED_ALLELES, true);

        // Use getBaseString() (rather than the Allele itself) in order to create new Alleles with swapped
        // reference and non-variant attributes
        swappedBuilder.alleles(Arrays.asList(vc.getAlleles().get(1).getBaseString(), vc.getAlleles().get(0).getBaseString()));

        final Map<Allele, Allele> alleleMap = new HashMap<>();

        // A mapping from the old allele to the new allele, to be used when fixing the genotypes
        alleleMap.put(vc.getAlleles().get(0), swappedBuilder.getAlleles().get(1));
        alleleMap.put(vc.getAlleles().get(1), swappedBuilder.getAlleles().get(0));

        final GenotypesContext swappedGenotypes = GenotypesContext.create(vc.getGenotypes().size());
        for (final Genotype genotype : vc.getGenotypes()) {
            final List<Allele> swappedAlleles = new ArrayList<>();
            for (final Allele allele : genotype.getAlleles()) {
                if (allele.isNoCall()) {
                    swappedAlleles.add(allele);
                } else {
                    swappedAlleles.add(alleleMap.get(allele));
                }
            }
            // Flip AD
            final GenotypeBuilder builder = new GenotypeBuilder(genotype).alleles(swappedAlleles);
            if (genotype.hasAD() && genotype.getAD().length == 2) {
                final int[] ad = ArrayUtils.clone(genotype.getAD());
                ArrayUtils.reverse(ad);
                builder.AD(ad);
            } else {
                builder.noAD();
            }

            //Flip PL
            if (genotype.hasPL() && genotype.getPL().length == 3) {
                final int[] pl = ArrayUtils.clone(genotype.getPL());
                ArrayUtils.reverse(pl);
                builder.PL(pl);
            } else {
                builder.noPL();
            }
            swappedGenotypes.add(builder.make());
        }
        swappedBuilder.genotypes(swappedGenotypes);

        for (final String key : vc.getAttributes().keySet()) {
            if (annotationsToDrop.contains(key)) {
                swappedBuilder.rmAttribute(key);
            } else if (annotationsToReverse.contains(key) && !vc.getAttributeAsString(key, "").equals(VCFConstants.MISSING_VALUE_v4)) {
                final double attributeToReverse = vc.getAttributeAsDouble(key, -1);

                if (attributeToReverse < 0 || attributeToReverse > 1) {
                    log.warn("Trying to reverse attribute " + key +
                            " but found value that isn't between 0 and 1: (" + attributeToReverse + ") in variant " + vc + ". Results might be wrong.");
                }
                swappedBuilder.attribute(key, 1 - attributeToReverse);
            }
        }

        return swappedBuilder.make();
    }

    /**
     * Checks whether the reference allele in the provided variant context actually matches the reference sequence
     *
     * @param alleles           list of alleles from which to find the reference allele
     * @param referenceSequence the ref sequence
     * @param start             the start position of the actual indel
     * @param end               the end position of the actual indel
     * @return true if they match, false otherwise
     */
    protected static boolean referenceAlleleMatchesReferenceForIndel(final List<Allele> alleles,
                                                                     final ReferenceSequence referenceSequence,
                                                                     final int start,
                                                                     final int end) {
        final String refString = StringUtil.bytesToString(referenceSequence.getBases(), start - 1, end - start + 1);
        final Allele refAllele = alleles.stream().filter(Allele::isReference).findAny().orElseThrow(() -> new IllegalStateException("Error: no reference allele was present"));
        return (refString.equalsIgnoreCase(refAllele.getBaseString()));
    }

    /**
     * Normalizes and left aligns a {@link VariantContextBuilder}.
     * Note: this will modify the start/stop and alleles of this builder.
     * Also note: if the reference allele does not match the reference sequence, this method will throw an exception
     *
     * Based on Adrian Tan, Gon&ccedil;alo R. Abecasis and Hyun Min Kang. (2015)
     * Unified Representation of Genetic Variants. Bioinformatics.
     */
    protected static void leftAlignVariant(final VariantContextBuilder builder, final int start, final int end, final List<Allele> alleles, final ReferenceSequence referenceSequence) {

        // make sure that referenceAllele matches reference
        if (!referenceAlleleMatchesReferenceForIndel(alleles, referenceSequence, start, end)) {
            throw new IllegalArgumentException(String.format("Reference allele doesn't match reference at %s:%d-%d", referenceSequence.getName(), start, end));
        }

        boolean changesInAlleles = true;

        final Map<Allele, byte[]> alleleBasesMap = new HashMap<>();

        // Put each allele into the alleleBasesMap unless it is a spanning deletion.
        // Spanning deletions are dealt with as a special case later in fixedAlleleMap.
        alleles.stream().filter(a->!a.equals(Allele.SPAN_DEL)).forEach(a -> alleleBasesMap.put(a, a.getBases()));

        int theStart = start;
        int theEnd = end;

        // 1. while changes in alleles do
        while (changesInAlleles) {

            changesInAlleles = false;
            // 2. if alleles end with the same nucleotide then
            if (alleleBasesMap.values().stream()
                    .collect(Collectors.groupingBy(a -> a[a.length - 1], Collectors.toSet()))
                    .size() == 1 && theEnd > 1) {
                // 3. truncate rightmost nucleotide of each allele
                for (final Allele allele : alleleBasesMap.keySet()) {
                    alleleBasesMap.put(allele, truncateBase(alleleBasesMap.get(allele), true));
                }
                changesInAlleles = true;
                theEnd--;
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
                    final byte extraBase = (theStart > 1) ?
                            referenceSequence.getBases()[theStart - 2] :
                            referenceSequence.getBases()[theEnd];

                    alleleBasesMap.put(allele, extendOneBase(alleleBasesMap.get(allele), extraBase));
                }
                changesInAlleles = true;
                theStart--;

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
            theStart++;
        }

        builder.start(theStart);
        builder.stop(theEnd);

        final Map<Allele, Allele> fixedAlleleMap = alleleBasesMap.entrySet().stream()
                .collect(Collectors.toMap(Map.Entry::getKey, me -> Allele.create(me.getValue(), me.getKey().isReference())));

        // A left aligned spanning deletion is itself a spanning deletion
        // Since this is not handled by the code above, it must be handled as
        // a special case here.
        fixedAlleleMap.put(Allele.SPAN_DEL, Allele.SPAN_DEL);

        //retain original order:
        List<Allele> fixedAlleles = alleles.stream()
                .map(fixedAlleleMap::get)
                .collect(Collectors.toList());

        builder.alleles(fixedAlleles);
    }

    private static byte[] truncateBase(final byte[] allele, final boolean truncateRightmost) {
        return Arrays.copyOfRange(allele, truncateRightmost ? 0 : 1, truncateRightmost ?
                allele.length - 1 :
                allele.length);
    }

    //creates a new byte array with the base added at the beginning
    private static byte[] extendOneBase(final byte[] bases, final byte base) {

        final byte[] newBases = new byte[bases.length + 1];

        System.arraycopy(bases, 0, newBases, 1, bases.length);
        newBases[0] = base;

        return newBases;
    }
}
