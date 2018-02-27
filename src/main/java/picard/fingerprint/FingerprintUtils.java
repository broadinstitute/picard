/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
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

package picard.fingerprint;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * A set of utilities used in the fingerprinting environment
 *
 * @author Yossi Farjoun
 */
public class FingerprintUtils {

    /**
     * A function that takes a Fingerprint and writes it as a VCF to a file
     *
     * @param fingerprint               the fingerprint to write
     * @param outputFile                the file to write to
     * @param referenceSequenceFileName the reference sequence (file)
     * @param sample                    the sample name to use in the vcf
     * @param source                    a "source" comment to use in the VCF
     * @throws IOException
     */
    public static void writeFingerPrint(final Fingerprint fingerprint,
                                        final File outputFile,
                                        final File referenceSequenceFileName,
                                        final String sample,
                                        final String source) throws IOException {

        try (final ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceSequenceFileName);
             final VariantContextWriter variantContextWriter = getVariantContextWriter(outputFile, referenceSequenceFileName, sample, source, ref)) {

            createVCSetFromFingerprint(fingerprint, ref, sample).forEach(variantContextWriter::add);
        }
    }

    private static VariantContextWriter getVariantContextWriter(final File outputFile,
                                                                final File referenceSequenceFileName,
                                                                final String sample,
                                                                final String source,
                                                                final ReferenceSequenceFile ref) {
        final VariantContextWriter variantContextWriter = new VariantContextWriterBuilder()
                .setReferenceDictionary(ref.getSequenceDictionary())
                .setOutputFile(outputFile).build();

        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFHeaderLine("reference", referenceSequenceFileName.getAbsolutePath()));
        lines.add(new VCFHeaderLine("source", source));
        lines.add(new VCFHeaderLine("fileDate", new Date().toString()));

        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_PL_KEY));
        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS));
        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.DEPTH_KEY));

        final VCFHeader header = new VCFHeader(lines, Collections.singletonList(sample));
        header.setSequenceDictionary(ref.getSequenceDictionary());
        variantContextWriter.writeHeader(header);
        return variantContextWriter;
    }

    /**
     * A utility function that takes a fingerprint and returns a VariantContextSet with variants representing the haplotypes in the fingerprint
     *
     * @param fingerPrint A fingerprint
     * @param reference   A reference sequence that will be used to create the VariantContexts
     * @param sample      A sample name that will be used for the genotype field
     * @return VariantContextSet with variants representing the haplotypes in the fingerprint
     */
    public static VariantContextSet createVCSetFromFingerprint(final Fingerprint fingerPrint, final ReferenceSequenceFile reference, final String sample) {

        final VariantContextSet variantContexts = new VariantContextSet(reference.getSequenceDictionary());

        // check that the same snp name isn't twice in the fingerprint.
        fingerPrint.values().stream()
                .map(hp -> hp.getRepresentativeSnp().getName())
                .filter(Objects::nonNull)
                .filter(n -> !n.equals(""))
                .collect(Collectors.groupingBy(Function.identity(), Collectors.counting()))
                .entrySet()
                .stream()
                .filter(e -> e.getValue() > 1)
                .findFirst()
                .ifPresent(e -> {
                    throw new IllegalArgumentException("Found same SNP name twice (" + e.getKey() + ") in fingerprint. Cannot create a VCF.");
                });

        // convert all the haplotypes to variant contexts and add them to the set.
        fingerPrint.values().stream()
                .map(hp -> getVariantContext(reference, sample, hp))
                .forEach(variantContexts::add);

        return variantContexts;
    }

    private static VariantContext getVariantContext(final ReferenceSequenceFile reference,
                                                    final String sample,
                                                    final HaplotypeProbabilities haplotypeProbabilities) {
        final Snp snp = haplotypeProbabilities.getRepresentativeSnp();
        final byte refAllele = StringUtil.toUpperCase(reference.getSubsequenceAt(
                snp.getChrom(),
                snp.getPos(),
                snp.getPos()).getBases()[0]);

        final Allele allele1 = Allele.create(snp.getAllele1(), snp.getAllele1() == refAllele);
        final Allele allele2 = Allele.create(snp.getAllele2(), snp.getAllele2() == refAllele);
        final List<Allele> alleles = Arrays.asList(allele1, allele2);

        final Genotype gt = new GenotypeBuilder()
                .DP(haplotypeProbabilities.getTotalObs())
                .noAttributes()
                .PL(haplotypeProbabilities.getLogLikelihoods())
                .AD(new int[]{haplotypeProbabilities.getObsAllele1(), haplotypeProbabilities.getObsAllele2()})
                .name(sample)
                .make();
        try {
            return new VariantContextBuilder(
                    snp.getName(),
                    snp.getChrom(),
                    snp.getPos(),
                    snp.getPos(),
                    alleles)
                    .log10PError(VariantContext.NO_LOG10_PERROR)
                    .genotypes(gt)
                    .unfiltered().make();
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException(String.format("Trouble creating variant at %s-%d", snp.getChrom(), snp.getPos()), e);
        }
    }

    /**
     * A class that holds VariantContexts sorted by genomic position
     */
    public static class VariantContextSet extends TreeSet<VariantContext> {
        VariantContextSet(final SAMSequenceDictionary dict) {
            super((lhs, rhs) -> {
                final int lhsContig = dict.getSequenceIndex(lhs.getContig());
                final int rhsContig = dict.getSequenceIndex(rhs.getContig());

                final int retval = lhsContig - rhsContig;
                if (retval != 0) return retval;

                return lhs.getStart() - rhs.getStart();
            });
        }
    }
}
