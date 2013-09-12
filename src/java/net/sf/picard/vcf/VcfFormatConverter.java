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

package net.sf.picard.vcf;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.picard.util.ProgressLogger;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.CloserUtil;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.VCFFileReader;
import org.broadinstitute.variant.vcf.VCFHeader;

import java.io.File;
import java.util.EnumSet;

/**
 * Converts an ASCII VCF file to a binary BCF or vice versa.
 *
 * @author jgentry@broadinstitute.org
 */
public class VcfFormatConverter extends CommandLineProgram {
    // The following attributes define the command-line arguments
    public static final Log LOG = Log.getInstance(VcfFormatConverter.class);
    
    @Usage
    public String USAGE = getStandardUsagePreamble() +
		    "Convert a VCF file to a BCF file, or BCF to VCF.\n" + "" +
            "Input and output formats are determined by file extension.";

    @Option(doc="The BCF or VCF input file. The file format is determined by file extension.", shortName= StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Option(doc="The BCF or VCF output file. The file format is determined by file extension.", shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

	@Option(doc="Fail if an index is not available for the input VCF/BCF")
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
        
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);

	    final VCFFileReader reader = new VCFFileReader(INPUT, REQUIRE_INDEX);
	    final VCFHeader header = new VCFHeader(reader.getFileHeader());
	    final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
	    if (CREATE_INDEX && sequenceDictionary == null) {
		    throw new PicardException("A sequence dictionary must be available in the input file when creating indexed output.");
	    }
	    final EnumSet<Options> options = CREATE_INDEX ? EnumSet.of(Options.INDEX_ON_THE_FLY) : EnumSet.noneOf(Options.class);
        final VariantContextWriter writer = VariantContextWriterFactory.create(OUTPUT, sequenceDictionary, options);

        writer.writeHeader(header);

	    final CloseableIterator<VariantContext> iterator = reader.iterator();
	    while (iterator.hasNext()) {
		    final VariantContext context = iterator.next();
            writer.add(context);
            progress.record(context.getChr(), context.getStart());
        }

	    CloserUtil.close(iterator);
	    CloserUtil.close(reader);
        writer.close();

        return 0;
    }
}
