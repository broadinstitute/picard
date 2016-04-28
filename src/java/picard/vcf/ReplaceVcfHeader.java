/*
 * The MIT License
 *
 * Copyright (c) 2016 Nils Homer
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

import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VcfOrBcf;

import java.io.File;
import java.util.ArrayList;

@CommandLineProgramProperties(
        usage =  ReplaceVcfHeader.USAGE_SUMMARY + ReplaceVcfHeader.USAGE_DETAILS,
        usageShort = ReplaceVcfHeader.USAGE_SUMMARY,
        programGroup = VcfOrBcf.class
)
public class ReplaceVcfHeader extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Replaces the VCFHeader in a VCF or BCF file.  ";
    static final String USAGE_DETAILS = "This tool makes it possible to replace the header of a VCF or BCF file with the header of another" +
            "file.  The samples must be the same and in the same order.<br /><br />" +
            "Note that validation is minimal, so it is up to the user to ensure that all the elements referred to in the VCFHeader " +
            "are present in the new header. " +
            "<br />" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar ReplaceVcfHeader \\<br />" +
            "      I=input_1.vcf \\<br />" +
            "      HEADER=input_2.vcf \\<br />" +
            "      O=bam_with_new_head.vcf" +
            "</pre>" +
            "<hr />";
    @Option(doc = "VCF file from which VariantContexts will be read.", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Option(doc = "VCF file from which VCFHeader will be read.")
    public File HEADER;

    @Option(doc = "VCFHeader from HEADER file will be written to this file, followed by VariantContexts from INPUT file",
            shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    public static void main(final String[] argv) {
        new ReplaceVcfHeader().instanceMainWithExit(argv);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(HEADER);
        IOUtil.assertFileIsWritable(OUTPUT);

        final VCFFileReader in = new VCFFileReader(INPUT);
        final VCFHeader header = in.getFileHeader();

        // Check that the samples are the same
        final ArrayList<String> inSamples = in.getFileHeader().getSampleNamesInOrder();
        final ArrayList<String> headerSamples = header.getSampleNamesInOrder();
        if (inSamples.size() != headerSamples.size()) {
            throw new IllegalArgumentException("The same number of samples must be found in the input and header files.");
        }
        for (int i = 0; i < inSamples.size(); i++) {
            if (!inSamples.get(i).equals(headerSamples.get(i))) {
                throw new IllegalArgumentException(String.format("The %dth sample was different: '%s' != '%s'",
                        i + 1, inSamples.get(i), headerSamples.get(i)));
            }
        }

        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                .modifyOption(Options.INDEX_ON_THE_FLY, CREATE_INDEX)
                .setReferenceDictionary(header.getSequenceDictionary())
                .setOutputFile(OUTPUT);
        final VariantContextWriter out = builder.build();

        out.writeHeader(header);

        for (final VariantContext ctx : in) {
            out.add(ctx);
        }

        out.close();
        in.close();

        return 0;
    }
}
