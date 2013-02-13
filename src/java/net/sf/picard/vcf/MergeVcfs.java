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
package net.sf.picard.vcf;

import net.sf.picard.cmdline.CommandLineParser;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.picard.util.MergingIterator;
import net.sf.picard.util.ProgressLogger;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.List;

/**
 * Combines multiple VCF files into a single file. Input files must be sorted by their contigs
 * and, within contigs, by start position. Throws IllegalArgumentException if the contig lists
 * are not present in the input files, are not identical or if the sample lists are not the
 * same; this class uses the GATK to merge headers, which may throw exceptions if the headers
 * cannot be merged. See VCFUtils.smartMergeHeaders for details.
 *
 * An index file is created for the output file by default. Using an output file name with a
 * ".gz" extension will create gzip-compressed output.
 */
public class MergeVcfs extends CommandLineProgram {

	@Usage
	public final String USAGE =
			CommandLineParser.getStandardUsagePreamble(getClass()) +
			"Merges multiple SAM/BAM files into one file.\n";

	@Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="VCF input files", minElements=1)
	public List<File> INPUT;

	@Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="VCF file to write merged result to")
	public File OUTPUT;

	@Option(shortName="D", doc="The index sequence dictionary (required if CREATE_INDEX=true)", optional = true)
	public File SEQUENCE_DICTIONARY;

	private final Log log = Log.getInstance(MergeVcfs.class);

	public static void main(final String[] argv) {
		new MergeVcfs().instanceMainWithExit(argv);
	}

	public MergeVcfs() {
		this.CREATE_INDEX = true;
	}

	@Override
	protected int doWork() {
		final ProgressLogger progress = new ProgressLogger(log, 10000);
		final List<String> sampleList = new ArrayList<String>();
		final Collection<CloseableIterator<VariantContext>> iteratorCollection = new ArrayList<CloseableIterator<VariantContext>>(INPUT.size());
		final Collection<VCFHeader> headers = new HashSet<VCFHeader>(INPUT.size());

		VariantContextComparator variantContextComparator = null;

		for (final File file : INPUT) {
			IoUtil.assertFileIsReadable(file);
			final VariantContextIterator variantIterator = new VariantContextIterator(file);
			final VCFHeader header = variantIterator.getHeader();
			if (variantContextComparator == null) {
				variantContextComparator = new VariantContextComparator(header.getContigLines());
			} else {
				if ( ! variantContextComparator.isCompatible(header.getContigLines())) {
					throw new IllegalArgumentException(
							"The contig entries in input file " + file.getAbsolutePath() + " are not compatible with the others.");
				}
			}

			if (sampleList.isEmpty()) {
				sampleList.addAll(header.getSampleNamesInOrder());
			} else {
				if ( ! sampleList.equals(header.getSampleNamesInOrder())) {
					throw new IllegalArgumentException("Input file " + file.getAbsolutePath() + " has sample entries that don't match the other files.");
				}
			}

			headers.add(header);
			iteratorCollection.add(variantIterator);
		}

		final EnumSet<Options> options = CREATE_INDEX ? EnumSet.of(Options.INDEX_ON_THE_FLY) : EnumSet.noneOf(Options.class);
		final SAMSequenceDictionary sequenceDictionary =
				SEQUENCE_DICTIONARY != null ? VariantContextUtils.getSequenceDictionary(SEQUENCE_DICTIONARY) : null;
		final VariantContextWriter out = VariantContextUtils.getConditionallyCompressingWriter(OUTPUT, sequenceDictionary, options);

		out.writeHeader(new VCFHeader(VCFUtils.smartMergeHeaders(headers, false), sampleList));

		final MergingIterator<VariantContext> mergingIterator = new MergingIterator<VariantContext>(variantContextComparator, iteratorCollection);
		while (mergingIterator.hasNext()) {
			final VariantContext context = mergingIterator.next();
			out.add(context);
			progress.record(context.getChr(), context.getStart());
		}

		out.close();
		return 0;
	}

	protected String[] customCommandLineValidation() {
		if (this.CREATE_INDEX && (this.SEQUENCE_DICTIONARY == null)) {
			return new String[] { "If CREATE_INDEX is set a sequence dictionary must be specified." };
		}
		return null;
	}
}
