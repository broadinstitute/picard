/*
 * The MIT License
 *
 * Copyright (c) 2011-2016 The Broad Institute
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
package picard.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import static htsjdk.samtools.SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX;
import static htsjdk.samtools.SAMRecord.NO_ALIGNMENT_START;

/**
 * Reorders a SAM/BAM input file according to the order of contigs in a second reference file.
 *
 * <h3>Summary</h3>
 * Not to be confused with SortSam which sorts a SAM or BAM file with a valid sequence dictionary,
 * ReorderSam reorders reads in a SAM/BAM file to match the contig ordering in a provided reference file,
 * as determined by exact name matching of contigs.  Reads mapped to contigs absent in the new
 * reference are unmapped. Runs substantially faster if the input is an indexed BAM file.
 *
 * <h3>Example</h3>
 * <pre>
 *     java -jar picard.jar ReorderSam \
 *          INPUT=sample.bam \
 *          OUTPUT=reordered.bam \
 *          SEQUENCE_DICTIONARY=reference_with_different_order.dict
 * </pre>
 *
 * <h3>Caveats</h3>
 * Note that REFERENCE_SEQUENCE is used for reading the INPUT, (e.g. when reading cram files) not for determining
 * the order of the OUTPUT. For that you must specify the SEQUENCE_DICTIONARY argument.
 *
 * @author mdepristo
 */
@CommandLineProgramProperties(
        summary = "Not to be confused with SortSam which sorts a file with a valid sequence dictionary, " +
                "ReorderSam reorders reads in a SAM/BAM/CRAM file to match the contig ordering in a provided reference file, " +
                "as determined by exact name matching of contigs.  Reads mapped to contigs absent in the new " +
                "reference are unmapped. Runs substantially faster if the input is an indexed BAM file." +
                "\n" +
                "Example\n" +
                "\n" +
                " java -jar picard.jar ReorderSam \\\n" +
                "      INPUT=sample.bam \\\n" +
                "      OUTPUT=reordered.bam \\\n" +
                "      SEQUENCE_DICTIONARY=reference_with_different_order.dict\n",
        oneLineSummary = "Reorders reads in a SAM or BAM file to match ordering in a second reference file.",
        programGroup = ReadDataManipulationProgramGroup.class)

@DocumentedFeature
public class ReorderSam extends CommandLineProgram {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input file (SAM/BAM/CRAM) to extract reads from.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (SAM/BAM/CRAM) to write extracted reads to.")
    public File OUTPUT;

    @Argument(shortName = StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME,
            doc = "A Sequence Dictionary for the OUTPUT file (can be read from one of the " +
                    "following file types (SAM, BAM, CRAM, VCF, BCF, Interval List, Fasta, or Dict)")
    public File SEQUENCE_DICTIONARY;

    @Argument(shortName = "S", doc = "If true, allows only a partial overlap of the original contigs with the new reference " +
            "sequence contigs.  By default, this tool requires a corresponding contig in the new " +
            "reference for each read contig")
    public boolean ALLOW_INCOMPLETE_DICT_CONCORDANCE = false;

    @Argument(shortName = "U", doc = "If true, then permits mapping from a read contig to a new reference contig with the " +
            "same name but a different length.  Highly dangerous, only use if you know what you " +
            "are doing.")
    public boolean ALLOW_CONTIG_LENGTH_DISCORDANCE = false;

    private final Log log = Log.getInstance(ReorderSam.class);

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(SEQUENCE_DICTIONARY);
        IOUtil.assertFileIsWritable(OUTPUT);

        try (SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT)) {

            final SAMSequenceDictionary outputDictionary = SAMSequenceDictionaryExtractor.extractDictionary(SEQUENCE_DICTIONARY.toPath());

            if (outputDictionary == null) {
                log.error("No reference sequence dictionary found. Aborting.  You can create a sequence dictionary for the reference fasta using the CreateSequenceDictionary command.");
                return 1;
            }

            printDictionary("SAM/BAM/CRAM file", in.getFileHeader().getSequenceDictionary());
            printDictionary("Reference", outputDictionary);
            final Map<Integer, Integer> newOrder;
            try {
                newOrder = buildSequenceDictionaryMap(outputDictionary, in.getFileHeader().getSequenceDictionary());
            } catch (PicardException e) {
                log.error(e);
                return 1;
            }
            // has to be after we create the newOrder
            final SAMFileHeader outHeader = in.getFileHeader().clone();
            outHeader.setSequenceDictionary(outputDictionary);

            log.info("Writing reads...");
            if (in.hasIndex()) {
                try (SAMFileWriter out = new SAMFileWriterFactory().makeWriter(outHeader, false, OUTPUT, REFERENCE_SEQUENCE)) {

                    // write the reads in contig order
                    for (final SAMSequenceRecord contig : in.getFileHeader().getSequenceDictionary().getSequences()) {
                        log.info("writing the reads from " + contig.getSequenceName());
                        final SAMRecordIterator it = in.query(contig.getSequenceName(), 0, 0, false);
                        writeReads(out, it, newOrder, contig.getSequenceName());
                    }

                    // don't forget the unmapped reads
                    writeReads(out, in.queryUnmapped(), newOrder, "unmapped");
                }
            } else {
                try (SAMFileWriter out = new SAMFileWriterFactory().makeWriter(outHeader, false, OUTPUT, REFERENCE_SEQUENCE)) {
                    writeReads(out, in.iterator(), newOrder, "All reads");
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return 0;
    }

    /**
     * Low-level helper function that returns the new reference index for oldIndex according to the
     * ordering map newOrder.  Read is provided in case an error occurs, so that an informative message
     * can be made.
     */
    private int newOrderIndex(SAMRecord read, int oldIndex, Map<Integer, Integer> newOrder) {
        if (oldIndex == NO_ALIGNMENT_REFERENCE_INDEX) {
            return NO_ALIGNMENT_REFERENCE_INDEX; // unmapped read
        } else {
            final Integer n = newOrder.get(oldIndex);

            if (n == null) {
                throw new PicardException("BUG: no mapping found for read " + read.getSAMString());
            } else {
                return n;
            }
        }
    }

    /**
     * Helper function that writes reads from iterator it into writer out, updating each SAMRecord along the way
     * according to the newOrder mapping from dictionary index -> index.  Name is used for printing only.
     */
    private void writeReads(final SAMFileWriter out,
                            final SAMRecordIterator it,
                            final Map<Integer, Integer> newOrder,
                            final String name) {
        long counter = 0;
        log.info("  Processing " + name);

        while (it.hasNext()) {
            counter++;
            final SAMRecord read = it.next();
            final int oldRefIndex = read.getReferenceIndex();
            final int oldMateIndex = read.getMateReferenceIndex();
            final int newRefIndex = newOrderIndex(read, oldRefIndex, newOrder);

            read.setHeader(out.getFileHeader());
            read.setReferenceIndex(newRefIndex);

            // read becoming unmapped
            if (oldRefIndex != NO_ALIGNMENT_REFERENCE_INDEX &&
                    newRefIndex == NO_ALIGNMENT_REFERENCE_INDEX) {
                read.setAlignmentStart(NO_ALIGNMENT_START);
                read.setReadUnmappedFlag(true);
                read.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
                read.setMappingQuality(SAMRecord.NO_MAPPING_QUALITY);
            }

            final int newMateIndex = newOrderIndex(read, oldMateIndex, newOrder);
            if (oldMateIndex != NO_ALIGNMENT_REFERENCE_INDEX &&
                    newMateIndex == NO_ALIGNMENT_REFERENCE_INDEX) { // mate becoming unmapped
                read.setMateAlignmentStart(NO_ALIGNMENT_START);
                read.setMateUnmappedFlag(true);
                read.setAttribute(SAMTag.MC.name(), null);      // Set the Mate Cigar String to null
            }
            read.setMateReferenceIndex(newMateIndex);

            out.addAlignment(read);
        }

        it.close();
        log.info("Wrote " + counter + " reads");
    }

    /**
     * Constructs a mapping from read sequence records index -> new sequence dictionary index for use in
     * reordering the reference index and mate reference index in each read.  -1 (SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX)
     * means unmapped.
     */
    private Map<Integer, Integer> buildSequenceDictionaryMap(final SAMSequenceDictionary refDict,
                                                             final SAMSequenceDictionary readsDict) {
        Map<Integer, Integer> newOrder = new HashMap<>();

        log.info("Reordering SAM/BAM/CRAM file:");
        for (final SAMSequenceRecord refRec : refDict.getSequences()) {
            final SAMSequenceRecord readsRec = readsDict.getSequence(refRec.getSequenceName());

            if (readsRec != null) {
                if (refRec.getSequenceLength() != readsRec.getSequenceLength()) {
                    String msg = String.format("Discordant contig lengths: read %s LN=%d, ref %s LN=%d",
                            readsRec.getSequenceName(), readsRec.getSequenceLength(),
                            refRec.getSequenceName(), refRec.getSequenceLength());
                    if (ALLOW_CONTIG_LENGTH_DISCORDANCE) {
                        log.warn(msg);
                    } else {
                        throw new PicardException(msg);
                    }
                }

                log.info(String.format("  Reordering read contig %s [index=%d] to => ref contig %s [index=%d]%n",
                        readsRec.getSequenceName(), readsRec.getSequenceIndex(),
                        refRec.getSequenceName(), refRec.getSequenceIndex()));
                newOrder.put(readsRec.getSequenceIndex(), refRec.getSequenceIndex());
            }
        }

        for (SAMSequenceRecord readsRec : readsDict.getSequences()) {
            if (!newOrder.containsKey(readsRec.getSequenceIndex())) {
                if (ALLOW_INCOMPLETE_DICT_CONCORDANCE) {
                    newOrder.put(readsRec.getSequenceIndex(), NO_ALIGNMENT_REFERENCE_INDEX);
                } else {
                    throw new PicardException("New reference sequence does not contain a matching contig for " + readsRec.getSequenceName());
                }
            }
        }

        return newOrder;
    }

    /**
     * Helper function to print out a sequence dictionary
     */
    private void printDictionary(String name, SAMSequenceDictionary dict) {
        log.info(name);
        for (final SAMSequenceRecord contig : dict.getSequences()) {
            log.info(String.format("   SN=%s LN=%d%n", contig.getSequenceName(), contig.getSequenceLength()));
        }
    }
}
