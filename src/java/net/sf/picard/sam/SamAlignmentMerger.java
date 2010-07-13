package net.sf.picard.sam;

import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;
import net.sf.samtools.util.Log;
import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SortingCollection;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

/**
 * Class that takes in a set of alignment information in SAM format and merges it with the set
 * of all reads for which alignment was attempted, stored in an unmapped SAM file.
 *
 * @author ktibbett@broadinstitute.org
 */
public class SamAlignmentMerger extends AbstractAlignmentMerger {

    private final Log log = Log.getInstance(SamAlignmentMerger.class);
    private final File alignedSamFile;
    private final boolean pairedRun;
    private final int maxGaps;
    private final SAMFileReader reader;

    /**
     * Constructor
     *
     * @param unmappedBamFile   The BAM file that was used as the input to the Maq aligner, which will
     *                          include info on all the reads that did not map
     * @param targetBamFile     The file to which to write the merged SAM records
     * @param referenceFasta    The reference sequence for the map files
     * @param programRecord     Program record for taget file SAMRecords created.
     * @param clipAdapters      Whether adapters marked in BAM file are clipped from the read
     * @param bisulfiteSequence Whether the reads are bisulfite sequence
     * @param pairedRun         Whether the run is a paired-end run
     * @param jumpingLibrary    Whether this is a jumping library
     * @param alignedReadsOnly  Whether to output only those reads that have alignment data
     * @param alignedSamFile    The SAM file with alignment information
     * @param maxGaps           The maximum number of insertions or deletions permitted in an
     *                          alignment.  Alignments with more than this many gaps will be ignored.
     */
    public SamAlignmentMerger (final File unmappedBamFile, final File targetBamFile, final File referenceFasta,
                 final SAMProgramRecord programRecord, final boolean clipAdapters, final boolean bisulfiteSequence,
                 final boolean pairedRun, final boolean jumpingLibrary, final boolean alignedReadsOnly,
                 final File alignedSamFile, final int maxGaps) {

        super(unmappedBamFile, targetBamFile, referenceFasta, clipAdapters, bisulfiteSequence,
              jumpingLibrary, alignedReadsOnly, programRecord);

        IoUtil.assertFileIsReadable(alignedSamFile);
        this.alignedSamFile = alignedSamFile;
        this.pairedRun = pairedRun;
        this.maxGaps = maxGaps;
        this.reader = new SAMFileReader(this.alignedSamFile);
        this.reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        if (programRecord == null) {
            if (reader.getFileHeader().getProgramRecords().size() == 1) {
                setProgramRecord(reader.getFileHeader().getProgramRecords().get(0));
            }
        }
        log.info("Processing SAM file: " + alignedSamFile.getAbsolutePath());
    }


    /**
     * Reads the aligned SAM records into a SortingCollection and returns an iterator over that collection
     */
    protected CloseableIterator<SAMRecord> getQuerynameSortedAlignedRecords() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        header.setSequenceDictionary(getSequenceDictionary());
        if (getProgramRecord() != null) {
            header.addProgramRecord(getProgramRecord());
        }

        SortingCollection<SAMRecord> alignmentSorter = SortingCollection.newInstance(SAMRecord.class,
                    new BAMRecordCodec(header), new SAMRecordQueryNameComparator(), MAX_RECORDS_IN_RAM);

        Map<String, SAMRecord> readNameToReadPending = new HashMap<String,SAMRecord>();
        int count = 0;
        for (CloseableIterator<SAMRecord> it = new CigarClippingIterator(reader.iterator()); it.hasNext();) {
            SAMRecord sam = it.next();
            sam.setReadName(cleanReadName(sam.getReadName()));
            if (pairedRun) {
                SAMRecord mate = readNameToReadPending.remove(sam.getReadName());
                if (mate != null) {
                    if ((!sam.getReadUnmappedFlag()) || (!mate.getReadUnmappedFlag())) {
                        final boolean proper = SamPairUtil.isProperPair(sam, mate, isJumpingLibrary());
                        sam.setProperPairFlag(proper);
                        mate.setProperPairFlag(proper);
                        SamPairUtil.setMateInfo(sam, mate, getHeader());

                        alignmentSorter.add(sam);
                        alignmentSorter.add(mate);
                        count += 2;
                    }
                }
                else {
                    readNameToReadPending.put(sam.getReadName(), sam);
                }
            }
            else {
                if (!sam.getReadUnmappedFlag()) {
                    alignmentSorter.add(sam);
                    count++;
                }
            }
            if (count > 0 && count % 1000000 == 0) {
                log.info("Read " + count + " records from alignment SAM/BAM.");
            }
        }
        log.info("Finished reading " + count + " total records from alignment SAM/BAM.");

        if (readNameToReadPending.size() > 0) {
            throw new PicardException("Unmatched reads left in pending map.");
        }
        
        reader.close();
        return alignmentSorter.iterator();
    }

    /**
     * For now, we only ignore those alignments that have more than <code>maxGaps</code> insertions
     * or deletions.
     */
    protected boolean ignoreAlignment(SAMRecord sam) {
        int gaps = 0;
        for (CigarElement el : sam.getCigar().getCigarElements()) {
            if (el.getOperator() == CigarOperator.I || el.getOperator() == CigarOperator.D ||
                el.getOperator() == CigarOperator.N) {
                gaps++;
            }
        }
        return gaps > maxGaps;
    }

}
