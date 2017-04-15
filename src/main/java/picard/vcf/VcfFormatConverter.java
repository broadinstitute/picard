/*
* Copyright (c) 2013 The Broad Institute
*
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
*
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package picard.vcf;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VcfOrBcf;

import java.io.File;

/**
 * Converts an ASCII VCF file to a binary BCF or vice versa.
 *
 * @author jgentry@broadinstitute.org
 */
@CommandLineProgramProperties(
        summary = VcfFormatConverter.USAGE_SUMMARY + VcfFormatConverter.USAGE_DETAILS,
        oneLineSummary = VcfFormatConverter.USAGE_SUMMARY,
        programGroup = VcfOrBcf.class
)
public class VcfFormatConverter extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Converts VCF to BCF or BCF to VCF.  ";
    static final String USAGE_DETAILS = "This tool converts files between the plain-text VCF format and its binary compressed equivalent, " +
            "BCF. Input and output formats are determined by file extensions specified in the file names. For best results," +
            " it is recommended to ensure that an index file is present and set the REQUIRE_INDEX option to true." +
            "<br />" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar VcfFormatConverter \\<br />" +
            "      I=input.vcf \\<br />" +
            "      O=output.bcf \\<br />" +
            "      REQUIRE_INDEX=true" +
            "</pre>"+
            "<hr />";
    // The following attributes define the command-line arguments
    public static final Log LOG = Log.getInstance(VcfFormatConverter.class);

    @Argument(doc="The BCF or VCF input file.", shortName= StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc="The BCF or VCF output file name.", shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

	@Argument(doc="Fail if an index is not available for the input VCF/BCF")
	public Boolean REQUIRE_INDEX = true;

    public static void main(final String[] argv) {
        new VcfFormatConverter().instanceMainWithExit(argv);
    }

	public VcfFormatConverter() {
		this.CREATE_INDEX = true;
	}

    @Override
    protected int doWork() {
        final ProgressLogger progress = new ProgressLogger(LOG, 10000);
        
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

	    final VCFFileReader reader = new VCFFileReader(INPUT, REQUIRE_INDEX);
	    final VCFHeader header = new VCFHeader(reader.getFileHeader());
	    final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
	    if (CREATE_INDEX && sequenceDictionary == null) {
		    throw new PicardException("A sequence dictionary must be available in the input file when creating indexed output.");
	    }

        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                .setOutputFile(OUTPUT)
                .setReferenceDictionary(sequenceDictionary);
        if (CREATE_INDEX)
            builder.setOption(Options.INDEX_ON_THE_FLY);
        else
            builder.unsetOption(Options.INDEX_ON_THE_FLY);
        final VariantContextWriter writer = builder.build();
        writer.writeHeader(header);
	    final CloseableIterator<VariantContext> iterator = reader.iterator();

	    while (iterator.hasNext()) {
		    final VariantContext context = iterator.next();
            writer.add(context);
            progress.record(context.getContig(), context.getStart());
        }

	    CloserUtil.close(iterator);
	    CloserUtil.close(reader);
        writer.close();

        return 0;
    }
}
