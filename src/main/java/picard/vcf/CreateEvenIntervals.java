/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.*;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VcfOrBcf;

import java.io.File;
import java.util.*;

/**
 * Traverses one or more gVCF files in order to generate a list of intervals that have approximately
 * the same number of non-reference positions in all of them.
 */
@CommandLineProgramProperties(
        usage = "Traverses one or more gVCF files in order to generate a list of intervals that have approximately " +
                "the same number of non-reference positions in all of them.  Assumes files contain sequence dictionaries.",
        usageShort = "Creates even interval list from one or more gVCF files",
        programGroup = VcfOrBcf.class
)
public class CreateEvenIntervals extends CommandLineProgram {

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="gVCF input files (File format is determined by file extension), or a file having a '.list' suffix containing the path to the files.", minElements=1)
    public List<File> INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output Picard Interval List")
    public File OUTPUT;

    @Option(shortName = "N", doc = "The number of even intervals to create", optional = true)
    private Integer NUM_INTERVALS = 1000;

    private final Log log = Log.getInstance(CreateEvenIntervals.class);
    private ProgressLogger progress = new ProgressLogger(log, 500000, "read");

    public static void main(final String[] argv) {
        new CreateEvenIntervals().instanceMainWithExit(argv);
    }

    public CreateEvenIntervals() { }

    @Override
    protected int doWork() {
        INPUT = IOUtil.unrollFiles(INPUT, IOUtil.VCF_EXTENSIONS);
        IOUtil.assertFileIsWritable(OUTPUT);

        // mapping from contig to a bitset representing "covered" positions
        final HashMap<String, BitSet> recordSets = new HashMap<>();
        IntervalList intervals = null;
        SAMSequenceDictionary dictionary = null;

        int fileCount = 0;
        for (final File file : INPUT) {
            IOUtil.assertFileIsReadable(file);
            final VCFFileReader fileReader = new VCFFileReader(file, false);
            log.info("Searching through gvcf file #" + ++fileCount);

            // initialize if this is the first time through
            if (intervals == null) {
                dictionary = fileReader.getFileHeader().getSequenceDictionary();
                initializeRecordSets(recordSets, dictionary);
                intervals = initializeIntervalList(dictionary);
            }
            // otherwise, confirm that this VCF contains the same sequence dictionary
            else {
                dictionary.assertSameDictionary(fileReader.getFileHeader().getSequenceDictionary());
            }

            markRecords(recordSets, fileReader.iterator());
        }

        if (intervals != null) {
            makeIntervals(recordSets, intervals, dictionary);
            // sort and write the output
            intervals.sorted().write(OUTPUT);
        }
        return 0;
    }

    private void initializeRecordSets(final HashMap<String, BitSet> recordSets,
                                      final SAMSequenceDictionary dictionary) {
        for (final SAMSequenceRecord contig : dictionary.getSequences()) {
            recordSets.put(contig.getSequenceName(), new BitSet(contig.getSequenceLength() + 1));
        }
    }

    private IntervalList initializeIntervalList(final SAMSequenceDictionary dictionary) {
        SAMFileHeader samFileHeader = new SAMFileHeader();
        samFileHeader.setSequenceDictionary(dictionary);
        return new IntervalList(samFileHeader);
    }

    private void markRecords(final HashMap<String, BitSet> recordSets,
                             final CloseableIterator<VariantContext> iterator) {
        String currentContig = null;
        BitSet currentRecordSet = null;

        while (iterator.hasNext()) {
            final VariantContext record = iterator.next();
            // ignore hom ref blocks (characterized by just the single <NON-REF> alternate allele)
            if (record.getAlternateAlleles().size() == 1)
                continue;

            if (currentContig == null || !currentContig.equals(record.getContig())) {
                currentContig = record.getContig();
                currentRecordSet = recordSets.get(currentContig);
            }

            currentRecordSet.set(record.getStart());
            progress.record(record.getContig(), record.getStart());
        }
    }

    private void makeIntervals(final HashMap<String, BitSet> recordSets,
                               final IntervalList intervals,
                               final SAMSequenceDictionary dictionary) {

        final long totalUniquePositions = countBits(recordSets);
        log.info("Unique positions: " + totalUniquePositions);
        final long recordsPerInterval = (totalUniquePositions / NUM_INTERVALS) + 1;
        log.info("Records per interval: " + recordsPerInterval);

        for (final Map.Entry<String, BitSet> contig : recordSets.entrySet()) {
            // don't output intervals spanning empty contigs
            if (contig.getValue().cardinality() == 0)
                continue;

            int start = 1;
            int currentIntervalCount = 0;
            int currentPosition = start;
            final int contigLength = dictionary.getSequence(contig.getKey()).getSequenceLength();

            while (currentPosition <= contigLength) {
                if (contig.getValue().get(currentPosition)) {
                    if (++currentIntervalCount == recordsPerInterval) {
                        intervals.add(new Interval(contig.getKey(), start, currentPosition));
                        start = currentPosition + 1;
                        currentIntervalCount = 0;
                    }
                }
                currentPosition++;
            }

            if (start < contigLength) {
                intervals.add(new Interval(contig.getKey(), start, contigLength));
            }
        }
    }

    private long countBits(final HashMap<String, BitSet> recordSets) {
        long count = 0L;
        for (final BitSet contig : recordSets.values()) {
            count += contig.cardinality();
        }
        return count;
    }

}
