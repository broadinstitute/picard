/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

/**
 * Reorders a SAM/BAM input file according to the order of contigs in a second reference sequence
 *
 * @author mdepristo
 */
@CommandLineProgramProperties(
        usage = "Not to be confused with SortSam which sorts a SAM or BAM file with a valid sequence dictionary, " +
                "ReorderSam reorders reads in a SAM/BAM file to match the contig ordering in a provided reference file, " +
                "as determined by exact name matching of contigs.  Reads mapped to contigs absent in the new " +
                "reference are dropped. Runs substantially faster if the input is an indexed BAM file.",
        usageShort = "Reorders reads in a SAM or BAM file to match ordering in reference",
        programGroup = SamOrBam.class
)
public class ReorderSam extends CommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input file (bam or sam) to extract reads from.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (bam or sam) to write extracted reads to.")
    public File OUTPUT;

    @Option(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence to reorder reads to match.  " +
            "A sequence dictionary corresponding to the reference fasta is required.  Create one with CreateSequenceDictionary.jar.")
    public File REFERENCE;

    @Option(shortName = "S", doc = "If true, then allows only a partial overlap of the BAM contigs with the new reference " +
            "sequence contigs.  By default, this tool requires a corresponding contig in the new " +
            "reference for each read contig")
    public boolean ALLOW_INCOMPLETE_DICT_CONCORDANCE = false;

    @Option(shortName = "U", doc = "If true, then permits mapping from a read contig to a new reference contig with the " +
            "same name but a different length.  Highly dangerous, only use if you know what you " +
            "are doing.")
    public boolean ALLOW_CONTIG_LENGTH_DISCORDANCE = false;

    private final Log log = Log.getInstance(ReorderSam.class);

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        new ReorderSam().instanceMainWithExit(argv);
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(REFERENCE);
        IOUtil.assertFileIsWritable(OUTPUT);

        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);

        ReferenceSequenceFile reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE);
        SAMSequenceDictionary refDict = reference.getSequenceDictionary();

        if (refDict == null) {
            log.error("No reference sequence dictionary found. Aborting.  You can create a sequence dictionary for the reference fasta using CreateSequenceDictionary.jar.");
            CloserUtil.close(in);
            return 1;
        }

        printDictionary("SAM/BAM file", in.getFileHeader().getSequenceDictionary());
        printDictionary("Reference", refDict);
        Map<Integer, Integer> newOrder = buildSequenceDictionaryMap(refDict, in.getFileHeader().getSequenceDictionary());

        // has to be after we create the newOrder
        SAMFileHeader outHeader = in.getFileHeader().clone();
        outHeader.setSequenceDictionary(refDict);

        log.info("Writing reads...");
        if (in.hasIndex()) {
            final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(outHeader, true, OUTPUT);

            // write the reads in contig order
            for (final SAMSequenceRecord contig : refDict.getSequences()) {
                final SAMRecordIterator it = in.query(contig.getSequenceName(), 0, 0, false);
                writeReads(out, it, newOrder, contig.getSequenceName());
            }
            // don't forget the unmapped reads
            writeReads(out, in.queryUnmapped(), newOrder, "unmapped");
            out.close();
        } else {
            SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(outHeader, false, OUTPUT);
            writeReads(out, in.iterator(), newOrder, "All reads");
            out.close();
        }

        // cleanup
        CloserUtil.close(in);
        return 0;
    }

    /**
     * Low-level helper function that returns the new reference index for oldIndex according to the
     * ordering map newOrder.  Read is provided in case an error occurs, so that an informative message
     * can be made.
     */
    private int newOrderIndex(SAMRecord read, int oldIndex, Map<Integer, Integer> newOrder) {
        if (oldIndex == -1)
            return -1; // unmapped read
        else {
            final Integer n = newOrder.get(oldIndex);

            if (n == null) throw new PicardException("BUG: no mapping found for read " + read.format());
            else return n;
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

            final int newMateIndex = newOrderIndex(read, oldMateIndex, newOrder);
            if (oldMateIndex != -1 && newMateIndex == -1) { // becoming unmapped
                read.setMateAlignmentStart(0);
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
     * reordering the reference index and mate reference index in each read.  -1 means unmapped.
     */
    private Map<Integer, Integer> buildSequenceDictionaryMap(final SAMSequenceDictionary refDict,
                                                             final SAMSequenceDictionary readsDict) {
        Map<Integer, Integer> newOrder = new HashMap<Integer, Integer>();

        log.info("Reordering SAM/BAM file:");
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
                if (ALLOW_INCOMPLETE_DICT_CONCORDANCE)
                    newOrder.put(readsRec.getSequenceIndex(), -1);
                else
                    throw new PicardException("New reference sequence does not contain a matching contig for " + readsRec.getSequenceName());
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
            log.info("  SN=%s LN=%d%n", contig.getSequenceName(), contig.getSequenceLength());
        }
    }
}
