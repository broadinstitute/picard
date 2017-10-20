package picard.fingerprint;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import picard.PicardException;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * A set of utilities used in the fingerprinting environment
 *
 * @author Yossi Farjoun
 */
public class FingerprintUtils {

    /**
     * A function that takes a Fingerprint and writes it as a VCF to a file
     *
     * @param fingerprint the fingerprint to write
     * @param outputFile the file to write to
     * @param referenceSequenceFileName the reference sequence (file)
     * @param sample the sample name to use in the vcf
     * @param source a "source" comment to use in the VCF
     * @throws IOException
     */
    public static void writeFingerPrint(final Fingerprint fingerprint, final File outputFile, final File referenceSequenceFileName, final String sample, final String source) throws IOException {

        final ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceSequenceFileName);

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

        variantContextWriter.writeHeader(new VCFHeader(lines, Collections.singletonList(sample)));

        final VariantContextSet variantContexts = createVCSetFromFingerprint(fingerprint, ref, sample);

        ref.close();

        for (final VariantContext vc : variantContexts) {
            variantContextWriter.add(vc);
        }
        variantContextWriter.close();
    }

    /**
     * A utility function that takes a fingerprint and returns a VariantContextSet with variants representing the haplotypes in the fingerprint
     *
     * @param fingerPrint A fingerprint
     * @param reference A reference sequence that will be used to create the VariantContexts
     * @param sample A sample name that will be used for the genotype field
     * @return VariantContextSet with variants representing the haplotypes in the fingerprint
     */
    static public VariantContextSet createVCSetFromFingerprint(final Fingerprint fingerPrint, final ReferenceSequenceFile reference, final String sample) {

        final VariantContextSet variantContexts = new VariantContextSet(reference.getSequenceDictionary());
        final Set<String> snpNames = new HashSet<>();

        for (final HaplotypeProbabilities haplotypeProbabilities : fingerPrint.values()) {

            final Snp snp = haplotypeProbabilities.getRepresentativeSnp();
            final byte refAllele = StringUtil.toUpperCase(reference.getSubsequenceAt(snp.getChrom(),
                    snp.getPos(),
                    snp.getPos()).getBases()[0]);

            final Allele allele1 = Allele.create(snp.getAllele1(), snp.getAllele1() == refAllele);
            final Allele allele2 = Allele.create(snp.getAllele2(), snp.getAllele2() == refAllele);
            final List<Allele> alleles = Arrays.asList(allele1, allele2);

            final VariantContextBuilder builder = new VariantContextBuilder(
                    snp.getName(),
                    snp.getChrom(),
                    snp.getPos(),
                    snp.getPos(),
                    alleles);

            final Genotype gt = new GenotypeBuilder()
                    .DP(haplotypeProbabilities.getTotalObs())
                    .noAttributes()
                    .PL(haplotypeProbabilities.getLogLikelihoods())
                    .AD(new int[]{haplotypeProbabilities.getObsAllele1(), haplotypeProbabilities.getObsAllele2()})
                    .name(sample)
                    .make();

            builder.log10PError(VariantContext.NO_LOG10_PERROR)
                    .genotypes(gt)
                    .unfiltered();

            String snpName = snp.getName();
            if (snpName != null && !snpName.equals("")) {
                if (snpNames.contains(snpName)) {
                    throw new PicardException("Found same SNP name twice (" + snpName + ") in fingerprint. Cannot create a VCF.");
                }
                snpNames.add(snpName);
                builder.id(snpName);
            }

            variantContexts.add(builder.make());
        }

        return variantContexts;
    }

    /**
     * A class that holds VariantContexts sorted by genomic position
     */
    public static class VariantContextSet extends TreeSet<VariantContext> {
        public VariantContextSet(final SAMSequenceDictionary dict) {
            super((lhs, rhs) -> {
                final int lhsContig = dict.getSequenceIndex(lhs.getContig());
                final int rhsContig = dict.getSequenceIndex(rhs.getContig());

                if (lhsContig < rhsContig) return -1;
                else if (rhsContig < lhsContig) return 1;
                else {
                    return lhs.getStart() - rhs.getStart();
                }
            });
        }
    }
}
