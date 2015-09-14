/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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
package picard.vcf.filter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VcfOrBcf;

import java.io.File;
import java.util.List;

/**
 * Applies a set of hard filters to Variants and to Genotypes within a VCF.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage = "Applies one or more hard filters to a VCF file to filter out genotypes and variants.",
        usageShort = "Hard filters a VCF.",
        programGroup = VcfOrBcf.class
)
public class FilterVcf extends CommandLineProgram {
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The INPUT VCF or BCF file.")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output VCF or BCF.")
    public File OUTPUT;

    @Option(doc="The minimum allele balance acceptable before filtering a site. Allele balance is calculated for heterozygotes as " +
            "the number of bases supporting the least-represented allele over the total number of base observations. Different heterozygote " +
            "genotypes at the same locus are measured independently. The locus is filtered if any allele balance is below the limit.")
    public double MIN_AB = 0.0d;

    @Option(doc="The minimum sequencing depth supporting a genotype before the genotype will be filtered out.")
    public int MIN_DP = 0;

    @Option(doc="The minimum genotype quality that must be achieved for a sample otherwise the genotype will be filtered out.")
    public int MIN_GQ = 0;

    @Option(doc="The maximum phred scaled fisher strand value before a site will be filtered out.")
    public double MAX_FS = Double.MAX_VALUE;

    @Option(doc="The minimum QD value to accept or otherwise filter out the variant.")
    public double MIN_QD = 0;

    /** Constructor to default to having index creation on. */
    public FilterVcf() { this.CREATE_INDEX = true; }

    // Stock main method
    public static void main(final String[] args) {
        new FilterVcf().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final List<VariantFilter>  variantFilters = CollectionUtil.makeList(new AlleleBalanceFilter(MIN_AB), new FisherStrandFilter(MAX_FS), new QdFilter(MIN_QD));
        final List<GenotypeFilter> genotypeFilters = CollectionUtil.makeList(new GenotypeQualityFilter(MIN_GQ), new DepthFilter(MIN_DP));
        final VCFFileReader in = new VCFFileReader(INPUT, false);
        final FilterApplyingVariantIterator iterator = new FilterApplyingVariantIterator(in.iterator(), variantFilters, genotypeFilters);

        final VCFHeader header = in.getFileHeader();
        // If the user is writing to a .bcf or .vcf, VariantContextBuilderWriter requires a Sequence Dictionary.  Make sure that the
        // Input VCF has one.
        final VariantContextWriterBuilder variantContextWriterBuilder = new VariantContextWriterBuilder();
        if (isVcfOrBcf(OUTPUT)) {
            final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
            if (sequenceDictionary == null) {
                throw new PicardException("The input vcf must have a sequence dictionary in order to create indexed vcf or bcfs.");
            }
            variantContextWriterBuilder.setReferenceDictionary(sequenceDictionary);
        }
        final VariantContextWriter out = variantContextWriterBuilder.setOutputFile(OUTPUT).build();
        header.addMetaDataLine(new VCFFilterHeaderLine("AllGtsFiltered", "Site filtered out because all genotypes are filtered out."));
        header.addMetaDataLine(new VCFFormatHeaderLine("FT", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Genotype filters."));
        for (final VariantFilter filter : variantFilters) {
            for (final VCFFilterHeaderLine line : filter.headerLines()) {
                header.addMetaDataLine(line);
            }
        }

        out.writeHeader(in.getFileHeader());

        while (iterator.hasNext()) {
            out.add(iterator.next());
        }

        out.close();
        in.close();
        return 0;
    }

    private boolean isVcfOrBcf(final File file) {
        final String fileName = file.getName();
        return fileName.endsWith(".vcf") || fileName.endsWith(".bcf");
    }
}
