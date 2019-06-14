/*
 * The MIT License
 *
 * Copyright (c) 2012 The Broad Institute
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

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.MergingIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VariantManipulationProgramGroup;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;

/**
 * Combines multiple variant files into a single variant file.
 *
 * <h3>Inputs</h3>
 *  <ul>
 *      <li>One or more input file in VCF format (can be gzipped, i.e. ending in ".vcf.gz", or binary compressed, i.e. ending in ".bcf").</li>
 *      <li>Optionally a sequence dictionary file (typically name ending in .dict) if the input VCF does not contain a
 *          complete contig list and if the output index is to be created (true by default).</li>
 *  </ul>
 *  <p>
 *  The input variant data must adhere to the following rules:
 *     <ul>
 *         <li>If there are samples, those must be the same across all input files.</li>
 *         <li>Input file headers must be contain compatible declarations for common annotations (INFO, FORMAT fields) and filters.</li>
 *         <li>Input files variant records must be sorted by their contig and position following the sequence dictionary provided
 *         or the header contig list.</li>
 *     </ul>
 * </p>
 * <p>
 *     You can either directly specify the list of files by specifying <code>INPUT</code> multiple times, or provide a list
 *     in a file with name ending in ".list" to <code>INPUT</code>.
 * </p>
 *
 * <h3>Outputs</h3>
 * A VCF sorted (i) according to the dictionary and (ii) by coordiante. 
 * <h3>Usage examples</h3>
 * <h4>Example 1:</h4>
 *     We combine several variant files in different formats, where at least one of them contains the contig list in its header.
 * <pre>
 *     java -jar picard.jar MergeVcfs \
 *          I=input_variants.01.vcf \
 *          I=input_variants.02.vcf.gz \
 *          O=output_variants.vcf.gz
 * </pre>
 * <h4>Example 2:</h4>
 *      Similar to example 1 but we use an input list file to specify the input files:
 * <pre>
 *     java -jar picard.jar MergeVcfs \
 *          I=input_variant_files.list \
 *          O=output_variants.vcf.gz
 * </pre>
 *
 * @since 1.0.1
 */
@CommandLineProgramProperties(
		oneLineSummary = MergeVcfs.SUMMARY_FIRST_SENTENCE,
        summary = MergeVcfs.SUMMARY,
        programGroup = VariantManipulationProgramGroup.class)
@DocumentedFeature
public class MergeVcfs extends CommandLineProgram {

	static final String SUMMARY_FIRST_SENTENCE = "Combines multiple variant files into a single variant file";
	static final String SUMMARY = "<p>" + SUMMARY_FIRST_SENTENCE + ".</p>" + 
			"<h3>Inputs</h3>" + 
			"<ul>" + 
			"      <li>One or more input file in VCF format (can be gzipped, i.e. ending in \".vcf.gz\", or binary compressed, i.e. ending in \".bcf\").</li>" + 
			"      <li>Optionally a sequence dictionary file (typically name ending in .dict) if the input VCF does not contain a" + 
			"          complete contig list and if the output index is to be created (true by default).</li>" + 
			"  </ul>" + 
			"  <p>" + 
			"  The input variant data must adhere to the following rules:</p>" + 
			"     <ul>" + 
			"         <li>If there are samples, those must be the same across all input files.</li>" + 
			"         <li>Input file headers must be contain compatible declarations for common annotations (INFO, FORMAT fields) and filters.</li>" + 
			"         <li>Input files variant records must be sorted by their contig and position following the sequence dictionary provided" + 
			"         or the header contig list.</li>" + 
			"     </ul>" + 
			" <p>You can either directly specify the list of files by specifying <code>INPUT</code> multiple times, or provide a list" + 
			"     in a file with name ending in \".list\" to <code>INPUT</code>.</p>" + 
			" <h3>Outputs</h3>" + 
			" <p>A VCF sorted (i) according to the dictionary and (ii) by coordinate.</p>" +
			" <h3>Usage examples</h3>" + 
			" <h4>Example 1:</h4>" + 
			" <p>We combine several variant files in different formats, where at least one of them contains the contig list in its header.</p>" + 
			" <pre>java -jar picard.jar MergeVcfs \\\n" + 
			"          I=input_variants.01.vcf \\\n" + 
			"          I=input_variants.02.vcf.gz \\\n" + 
			"          O=output_variants.vcf.gz</pre>" + 
			" <h4>Example 2:</h4>" + 
			" <p>Similar to example 1 but we use an input list file to specify the input files:</p>" + 
			" <pre>java -jar picard.jar MergeVcfs \\\n" + 
			"          I=input_variant_files.list \\\n" + 
			"          O=output_variants.vcf.gz</pre><hr/>";  
	  	
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
              doc="VCF or BCF input files (File format is determined by file extension), or a file having a '.list' suffix containing the path to the files, one per line.", minElements=1)
    public List<File> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The merged VCF or BCF file. File format is determined by file extension.")
    public File OUTPUT;

    @Argument(shortName = "D", doc = "The index sequence dictionary to use instead of the sequence dictionary in the input files", optional = true)
    public File SEQUENCE_DICTIONARY;

    private final static String SEQ_DICT_REQUIRED = "A sequence dictionary must be available (either through the input file or by setting it explicitly).";

    private final Log log = Log.getInstance(MergeVcfs.class);

    public MergeVcfs() {
        this.CREATE_INDEX = true;
    }

    @Override
    protected int doWork() {
        final ProgressLogger progress = new ProgressLogger(log, 10000);
        final List<String> sampleList = new ArrayList<String>();
        INPUT = IOUtil.unrollFiles(INPUT, IOUtil.VCF_EXTENSIONS);
        final Collection<CloseableIterator<VariantContext>> iteratorCollection = new ArrayList<CloseableIterator<VariantContext>>(INPUT.size());
        final Collection<VCFHeader> headers = new HashSet<VCFHeader>(INPUT.size());
        VariantContextComparator variantContextComparator = null;
        SAMSequenceDictionary sequenceDictionary = null;

        if (SEQUENCE_DICTIONARY != null) {
            sequenceDictionary = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(SEQUENCE_DICTIONARY).getFileHeader().getSequenceDictionary();
        }

        for (final File file : INPUT) {
            IOUtil.assertFileIsReadable(file);
            final VCFFileReader fileReader = new VCFFileReader(file, false);
            final VCFHeader fileHeader = fileReader.getFileHeader();
            if (fileHeader.getContigLines().isEmpty()) {
                if (sequenceDictionary == null) {
                    throw new IllegalArgumentException(SEQ_DICT_REQUIRED);
                } else {
                    fileHeader.setSequenceDictionary(sequenceDictionary);
                }
            }

            if (variantContextComparator == null) {
                variantContextComparator = fileHeader.getVCFRecordComparator();
            } else {
                if (!variantContextComparator.isCompatible(fileHeader.getContigLines())) {
                    throw new IllegalArgumentException(
                            "The contig entries in input file " + file.getAbsolutePath() + " are not compatible with the others.");
                }
            }

            if (sequenceDictionary == null) sequenceDictionary = fileHeader.getSequenceDictionary();

            if (sampleList.isEmpty()) {
                sampleList.addAll(fileHeader.getSampleNamesInOrder());
            } else {
                if (!sampleList.equals(fileHeader.getSampleNamesInOrder())) {
                    throw new IllegalArgumentException("Input file " + file.getAbsolutePath() + " has sample entries that don't match the other files.");
                }
            }

            headers.add(fileHeader);
            iteratorCollection.add(fileReader.iterator());
        }

        if (CREATE_INDEX && sequenceDictionary == null) {
            throw new PicardException(String.format("Index creation failed. %s", SEQ_DICT_REQUIRED));
        }

        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                .setOutputFile(OUTPUT)
                .setReferenceDictionary(sequenceDictionary);

        if (CREATE_INDEX) {
            builder.setOption(Options.INDEX_ON_THE_FLY);
        } else {
            builder.unsetOption(Options.INDEX_ON_THE_FLY);
        }
        final VariantContextWriter writer = builder.build();

        writer.writeHeader(new VCFHeader(VCFUtils.smartMergeHeaders(headers, false), sampleList));

        final MergingIterator<VariantContext> mergingIterator = new MergingIterator<VariantContext>(variantContextComparator, iteratorCollection);
        while (mergingIterator.hasNext()) {
            final VariantContext context = mergingIterator.next();
            writer.add(context);
            progress.record(context.getContig(), context.getStart());
        }

        CloserUtil.close(mergingIterator);
        writer.close();
        return 0;
    }
}
