/*
 * The MIT License
 *
 * Copyright (c) 2014-2017 The Broad Institute
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
 * 
 */
package picard.vcf.filter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.*;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.JavascriptVariantFilter;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VcfOrBcf;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Applies a set of hard filters to Variants and to Genotypes within a VCF.
 *
 * @author Tim Fennell
 * @author Pierre Lindenbaum added the javascript filter in 2016 
 */
@CommandLineProgramProperties(
        summary = "Applies one or more hard filters to a VCF file to filter out genotypes and variants.",
        oneLineSummary = "Hard filters a VCF.",
        programGroup = VcfOrBcf.class
)
public class FilterVcf extends CommandLineProgram {
    @Argument(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The INPUT VCF or BCF file.")
    public File INPUT;

    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output VCF or BCF.")
    public File OUTPUT;

    @Argument(doc="The minimum allele balance acceptable before filtering a site. Allele balance is calculated for heterozygotes as " +
            "the number of bases supporting the least-represented allele over the total number of base observations. Different heterozygote " +
            "genotypes at the same locus are measured independently. The locus is filtered if any allele balance is below the limit.")
    public double MIN_AB = 0.0d;

    @Argument(doc="The minimum sequencing depth supporting a genotype before the genotype will be filtered out.")
    public int MIN_DP = 0;

    @Argument(doc="The minimum genotype quality that must be achieved for a sample otherwise the genotype will be filtered out.")
    public int MIN_GQ = 0;

    @Argument(doc="The maximum phred scaled fisher strand value before a site will be filtered out.")
    public double MAX_FS = Double.MAX_VALUE;

    @Argument(doc="The minimum QD value to accept or otherwise filter out the variant.")
    public double MIN_QD = 0;

    @Argument(shortName = "JS", doc = "Filters a VCF file with a javascript expression interpreted by the java javascript engine. "
            + " The script puts the following variables in the script context: "
            + " 'variant' a VariantContext ( https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html ) and "
            + " 'header' a VCFHeader ( https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/vcf/VCFHeader.html )."
            + " Last value of the script should be a boolean to tell whether we should accept or reject the record.",
            optional = true)
    public File JAVASCRIPT_FILE = null;


    /** Constructor to default to having index creation on. */
    public FilterVcf() { this.CREATE_INDEX = true; }

    // Stock main method
    public static void main(final String[] args) {
        new FilterVcf().instanceMainWithExit(args);
    }

    final private Log log = Log.getInstance(FilterVcf.class);
    final private ProgressLogger progress = new ProgressLogger(log, 100_000, "Processed", "Variants");

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        VCFFileReader in = null;
        VariantContextWriter out = null;
        try {// try/finally used to close 'in' and 'out'
            in = new VCFFileReader(INPUT, false);
            final List<VariantFilter> variantFilters = new ArrayList<>(4);
            variantFilters.add(new AlleleBalanceFilter(MIN_AB));
            variantFilters.add(new FisherStrandFilter(MAX_FS));
            variantFilters.add(new QdFilter(MIN_QD));
            if (JAVASCRIPT_FILE != null) {
                try {
                    variantFilters.add(new VariantContextJavascriptFilter(JAVASCRIPT_FILE, in.getFileHeader()));
                } catch (final IOException error) {
                    throw new PicardException("javascript-related error", error);
                }
            }
            final List<GenotypeFilter> genotypeFilters = CollectionUtil.makeList(new GenotypeQualityFilter(MIN_GQ), new DepthFilter(MIN_DP));
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
            out = variantContextWriterBuilder.setOutputFile(OUTPUT).build();
            header.addMetaDataLine(new VCFFilterHeaderLine("AllGtsFiltered", "Site filtered out because all genotypes are filtered out."));
            header.addMetaDataLine(new VCFFormatHeaderLine("FT", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Genotype filters."));
            for (final VariantFilter filter : variantFilters) {
                filter.headerLines().forEach(header::addMetaDataLine);
            }

            out.writeHeader(in.getFileHeader());

            while (iterator.hasNext()) {
                final VariantContext vc = iterator.next();
                progress.record(vc.getContig(), vc.getStart());
                out.add(vc);
            }
            return 0;
        } finally {
            CloserUtil.close(out);
            CloserUtil.close(in);
        }
    }

    private boolean isVcfOrBcf(final File file) {
        final String fileName = file.getName();
        return fileName.endsWith(".vcf") || fileName.endsWith(".bcf");
    }
    
    /** Javascript filter implementing VariantFilter */
    private static class VariantContextJavascriptFilter
            extends JavascriptVariantFilter
            implements VariantFilter {
        /** filter name */
        private final String filterName;
        /** script file */
        private final File scriptFile;

        private VariantContextJavascriptFilter(final File scriptFile, final VCFHeader header) throws IOException {
            super(scriptFile, header);
            this.scriptFile = scriptFile;
           /* create filter name using file basename */
            String fname = IOUtil.basename(scriptFile);
            if (fname.isEmpty()) fname = "JSFILTER";
            this.filterName = fname;
        }

        /**
         * returns the filterName if the javascript doesn't accept the variant,
         * null otherwise
         */
        @Override
        public String filter(final VariantContext ctx) {
            return (super.accept(ctx) ? null : this.filterName );
        }

        @Override
        public List<VCFFilterHeaderLine> headerLines() {
            return CollectionUtil.makeList(
                    new VCFFilterHeaderLine(this.filterName, "Variant Filtered by Javascript file " + this.scriptFile));
        }
    }
}
