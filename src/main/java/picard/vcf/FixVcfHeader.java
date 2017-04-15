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
 *
 */

package picard.vcf;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VcfOrBcf;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Tool for replacing or fixing up a VCF header.
 *
 * @author Nils Homer
 */
@CommandLineProgramProperties(
        summary = FixVcfHeader.USAGE_SUMMARY + FixVcfHeader.USAGE_DETAILS,
        oneLineSummary = FixVcfHeader.USAGE_SUMMARY,
        programGroup = VcfOrBcf.class
)
public class FixVcfHeader extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Replaces or fixes a VCF header.";
    static final String USAGE_DETAILS =
            "This tool will either replace the header in the input VCF file (INPUT) with the given VCF header (HEADER) or will attempt to fill " +
            "in any field definitions that are missing in the input header by examining the variants in the input VCF file (INPUT).  In the  " +
            " latter case, this tool will perform two passes over the input VCF, and any INFO and FORMAT fields found in the VCF records but " +
            "not found in the input VCF header will be added to the output VCF header with dummy descriptions.<br />" +
            "<h4>Replace header usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar FixVcfHeader \\<br />" +
            "     I=input.vcf \\<br />" +
            "     O=fixed.vcf \\<br />" +
            "     HEADER=header.vcf" +
            "</pre>" +
            "<h4>Fix header usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar FixVcfHeader \\<br />" +
            "     I=input.vcf \\<br />" +
            "     O=fixed.vcf \\<br />" +
            "</pre>" +
            "<hr />";
    @Argument(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The input VCF/BCF file.")
    public File INPUT;

    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output VCF/BCF file.")
    public File OUTPUT;

    @Argument(shortName="N", doc="Check only the first N records when searching for missing INFO and FORMAT fields.", optional=true)
    public int CHECK_FIRST_N_RECORDS = -1;

    @Argument(shortName="H", doc="The replacement VCF header.", optional=true)
    public File HEADER = null;

    @Argument(doc="Enforce that the samples are the same (and in the same order) when replacing the VCF header.", optional=true)
    public boolean ENFORCE_SAME_SAMPLES = true;

    private final Log log = Log.getInstance(FixVcfHeader.class);

    // Stock main method
    public static void main(final String[] args) {
        new FixVcfHeader().instanceMainWithExit(args);
    }

    @Override
    protected String[] customCommandLineValidation() {
        if (HEADER != null && 0 <= CHECK_FIRST_N_RECORDS) return new String[]{"CHECK_FIRST_N_RECORDS should no be specified when HEADER is specified"};
        return super.customCommandLineValidation();
    }

    @Override protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        if (HEADER != null) IOUtil.assertFileIsReadable(HEADER);
        IOUtil.assertFileIsWritable(OUTPUT);

        final VCFFileReader reader = new VCFFileReader(INPUT, false);

        final VCFHeader outHeader;
        if (HEADER != null) { // read the header from the file
            final VCFFileReader headerReader = new VCFFileReader(HEADER, false);
            if (ENFORCE_SAME_SAMPLES) {
                outHeader = headerReader.getFileHeader();
                enforceSameSamples(reader.getFileHeader(), outHeader);
            }
            else {
                outHeader = new VCFHeader(headerReader.getFileHeader().getMetaDataInInputOrder(), reader.getFileHeader().getSampleNamesInOrder());
            }
            CloserUtil.close(headerReader);
            outHeader.getInfoHeaderLines()
                    .stream()
                    .filter(infoHeaderLine -> !reader.getFileHeader().hasInfoLine(infoHeaderLine.getID()))
                    .forEach(infoHeaderLine -> log.info("INFO line found in HEADER will be added to OUTPUT: " + infoHeaderLine.getID()));
            outHeader.getFormatHeaderLines()
                    .stream()
                    .filter(formatHeaderLine -> !reader.getFileHeader().hasInfoLine(formatHeaderLine.getID()))
                    .forEach(formatHeaderLine -> log.info("FORMAT line found in HEADER will be added to OUTPUT: " + formatHeaderLine.getID()));
        }
        else { // read the input
            final ProgressLogger progress = new ProgressLogger(log, 1000000, "read");
            final VCFFileReader in = new VCFFileReader(INPUT, false);
            final VCFHeader inHeader = in.getFileHeader();

            final Map<String, VCFInfoHeaderLine> infoHeaderLines = new HashMap<>();
            final Map<String, VCFFormatHeaderLine> formatHeaderLines = new HashMap<>();
            log.info("Reading in records to re-build the header.");
            for (final VariantContext ctx : in) {
                // INFO
                for (final Map.Entry<String, Object> attribute : ctx.getAttributes().entrySet()) {
                    final String id = attribute.getKey();
                    if (!inHeader.hasInfoLine(id) && !infoHeaderLines.containsKey(id)) {
                        log.info("Will add an INFO line with id: " + id);
                        infoHeaderLines.put(id, new VCFInfoHeaderLine(id, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Missing description: this INFO line was added by Picard's FixVCFHeader"));
                    }
                }
                // FORMAT
                for (final Genotype genotype : ctx.getGenotypes()) {
                    for (final Map.Entry<String, Object> attribute : genotype.getExtendedAttributes().entrySet()) {
                        final String id = attribute.getKey();
                        if (!inHeader.hasFormatLine(id) && !formatHeaderLines.containsKey(id)) {
                            log.info("Will add FORMAT line with id: " + id);
                            formatHeaderLines.put(id, new VCFFormatHeaderLine(id, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Missing description: this FORMAT line was added by Picard's FixVCFHeader"));
                        }
                    }
                }
                progress.record(ctx.getContig(), ctx.getStart());
                if (0 < CHECK_FIRST_N_RECORDS && CHECK_FIRST_N_RECORDS <= progress.getCount()) break;
            }
            CloserUtil.close(in);

            // create the output header
            final Set<VCFHeaderLine> headerLines = new HashSet<>(inHeader.getMetaDataInInputOrder());
            VCFStandardHeaderLines.addStandardFormatLines(headerLines, false, Genotype.PRIMARY_KEYS); // This is very frustrating to have to add
            headerLines.addAll(infoHeaderLines.values());
            headerLines.addAll(formatHeaderLines.values());
            outHeader = new VCFHeader(headerLines, inHeader.getSampleNamesInOrder());
            log.info("VCF header re-built.");
        }

        final VariantContextWriter writer = new VariantContextWriterBuilder().setOption(Options.INDEX_ON_THE_FLY)
                .setOutputFile(OUTPUT).setReferenceDictionary(outHeader.getSequenceDictionary()).build();
        writer.writeHeader(outHeader);

        log.info("Writing the output VCF.");
        final ProgressLogger progress = new ProgressLogger(log, 1000000, "read");
        for (final VariantContext ctx : reader) {
            writer.add(ctx);
            progress.record(ctx.getContig(), ctx.getStart());
        }
        writer.close();
        CloserUtil.close(reader);

        return 0;
    }

    private void enforceSameSamples(final VCFHeader readerHeader, final VCFHeader inputHeader) {
        final ArrayList<String> readerSamples = readerHeader.getSampleNamesInOrder();
        final ArrayList<String> inputSamples = inputHeader.getSampleNamesInOrder();
        if (readerSamples.size() != inputSamples.size()) {
            throw new PicardException("The input VCF had a different # of samples than the input VCF header.");
        }
        for (int i = 0; i < readerSamples.size(); i++) {
            if (!readerSamples.get(i).equals(inputSamples.get(i))) {
                throw new PicardException(String.format("Mismatch in the %dth sample: '%s' != '%s'", i, readerSamples.get(i), inputSamples.get(i)));
            }
        }
    }
}
