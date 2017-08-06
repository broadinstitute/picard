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

package picard.vcf;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VcfOrBcf;

import java.io.File;
import java.util.EnumSet;

@CommandLineProgramProperties(
        summary = RenameSampleInVcf.USAGE_SUMMARY + RenameSampleInVcf.USAGE_DETAILS,
        oneLineSummary = RenameSampleInVcf.USAGE_SUMMARY,
        programGroup = VcfOrBcf.class
)
public class RenameSampleInVcf extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Renames a sample within a VCF or BCF.  ";
    static final String USAGE_DETAILS = "This tool enables the user to rename a sample in either a VCF or BCF file.  " +
            "It is intended to change the name of a sample in a VCF prior to merging with VCF files in which one or more samples have " +
            "similar names. Note that the input VCF file must be single-sample VCF and that the NEW_SAMPLE_NAME is required." +
            "<br />" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar RenameSampleInVcf \\<br />" +
            "      I=input.vcf \\<br />" +
            "      O=renamed.vcf \\<br />" +
            "      NEW_SAMPLE_NAME=sample123" +
            "</pre>" +
            "<hr />";
    @Argument(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input single sample VCF.")
    public File INPUT;

    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output single sample VCF.")
    public File OUTPUT;

    @Argument(doc="Existing name of sample in VCF; if provided, asserts that that is the name of the extant sample name", optional = true)
    public String OLD_SAMPLE_NAME = null;

    @Argument(doc="New name to give sample in output VCF.")
    public String NEW_SAMPLE_NAME;


    public static void main(final String[] args) {
        new RenameSampleInVcf().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final VCFFileReader in = new VCFFileReader(INPUT);
        final VCFHeader header = in.getFileHeader();

        if (header.getGenotypeSamples().size() > 1) {
            throw new IllegalArgumentException("Input VCF must be single-sample.");
        }

        if (OLD_SAMPLE_NAME != null && !OLD_SAMPLE_NAME.equals(header.getGenotypeSamples().get(0))) {
            throw new IllegalArgumentException("Input VCF did not contain expected sample. Contained: " + header.getGenotypeSamples().get(0));
        }

        final EnumSet<Options> options = EnumSet.copyOf(VariantContextWriterBuilder.DEFAULT_OPTIONS);
        if (CREATE_INDEX) options.add(Options.INDEX_ON_THE_FLY); else options.remove(Options.INDEX_ON_THE_FLY);

        final VCFHeader outHeader = new VCFHeader(header.getMetaDataInInputOrder(), CollectionUtil.makeList(NEW_SAMPLE_NAME));
        final VariantContextWriter out = new VariantContextWriterBuilder()
                .setOptions(options)
                .setOutputFile(OUTPUT).setReferenceDictionary(outHeader.getSequenceDictionary()).build();
        out.writeHeader(outHeader);

        for (final VariantContext ctx : in) {
            out.add(ctx);
        }

        out.close();
        in.close();

        return 0;
    }
}
