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

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.picard.util.ProgressLogger;
import net.sf.samtools.util.CloserUtil;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;

import java.io.File;

/** Converts an ASCII VCF file to a binary BCF or vice versa
 *
 * @author jgentry@broadinstitute.org
 */
public class VcfFormatConverter extends CommandLineProgram {
    // The following attributes define the command-line arguments
    public static final Log LOG = Log.getInstance(VcfFormatConverter.class);
    
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Convert a VCF file to a BCF file, or BCF to VCF.\n" + "" +
            "Input and output formats are determined by file extension.";

    @Option(doc="The BCF or VCF file to parse.", shortName= StandardOptionDefinitions.INPUT_SHORT_NAME) public File INPUT;
    @Option(doc="The BCF or VCF output file. ", shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME) public File OUTPUT;

    public static void main(final String[] argv) {
        new VcfFormatConverter().instanceMainWithExit(argv);
    }

    @Override
    protected int doWork() {
        final ProgressLogger progress = new ProgressLogger(LOG, 10000);
        
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);

        final VariantContextIterator readerIterator = VariantContextIteratorFactory.create(INPUT);
        final VariantContextWriter writer = VariantContextWriterFactory.create(OUTPUT, null);

        writer.writeHeader(readerIterator.getHeader());

        while (readerIterator.hasNext()) {
            final VariantContext v = readerIterator.next();
            writer.add(v);
            progress.record(v.getChr(), v.getStart());
        }

        CloserUtil.close(readerIterator);
        CloserUtil.close(writer);
        return 0;
    }
}
