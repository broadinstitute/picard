/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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
package net.sf.picard.sam;

import net.sf.samtools.util.SequenceUtil;
import net.sf.picard.PicardException;

import java.util.*;

import net.sf.samtools.*;

/**
 * Merges SAMFileHeaders that have the same sequences into a single merged header
 * object while providing read group translation for cases where read groups
 * clash across input headers.
 *
 * @author Dave Tefft
 */
public class SamFileHeaderMerger {
    //Super Header to construct
    private final SAMFileHeader mergedHeader;
    private final Collection<SAMFileReader> readers;

    //Translation of old group ids to new group ids
    private final Map<SAMFileReader, Map<String, String>> samReadGroupIdTranslation =
            new HashMap<SAMFileReader, Map<String, String>>();

    //the read groups from different files use the same group ids
    private boolean hasReadGroupCollisions = false;

    //the program records from different files use the same program record ids
    private boolean hasProgramGroupCollisions = false;

    //Translation of old program group ids to new program group ids
    private Map<SAMFileReader, Map<String, String>> samProgramGroupIdTranslation = 
            new HashMap<SAMFileReader, Map<String, String>>();

    private boolean hasMergedSequenceDictionary = false;

    //Translation of old sequence dictionary ids to new dictionary ids
    private final Map<SAMFileReader, Map<Integer, Integer>> samSeqDictionaryIdTranslation =
            new HashMap<SAMFileReader, Map<Integer, Integer>>();

    // We don't always have access to the SAMFileReader in order to find the right mapping,
    // so also store the mapping using theSAMFileHeader
    private final Map<SAMFileHeader,  Map<Integer, Integer>> samSeqDictionaryIdTranslationViaHeader =
            new HashMap<SAMFileHeader, Map<Integer, Integer>>();

    //Letters to construct new ids from a counter
    private static final String ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";


    /**
     * Create SAMFileHeader with additional information.  Required that sequence dictionaries agree.
     *
     * @param readers sam file readers to combine
     * @param sortOrder sort order new header should have
     */
    public SamFileHeaderMerger(final Collection<SAMFileReader> readers, final SAMFileHeader.SortOrder sortOrder) {
        this(readers, sortOrder, false);
    }

    /**
     * Create SAMFileHeader with additional information.
     *
     * @param readers sam file readers to combine
     * @param sortOrder sort order new header should have
     * @param mergeDictionaries If true, merge sequence dictionaries in new header.  If false, require that
     * all input sequence dictionaries be identical.
     */
    public SamFileHeaderMerger(final Collection<SAMFileReader> readers, final SAMFileHeader.SortOrder sortOrder, final boolean mergeDictionaries) {
        this.readers = readers;
        this.mergedHeader = new SAMFileHeader();

        SAMSequenceDictionary sequenceDictionary = null;
        try {
            sequenceDictionary = getSequenceDictionary(readers);
            this.hasMergedSequenceDictionary = false;
        }
        catch (SequenceUtil.SequenceListsDifferException pe) {
            if (mergeDictionaries) {
                sequenceDictionary = mergeSequenceDictionaries(readers);
                this.hasMergedSequenceDictionary = true;
            }
            else {
                throw pe;
            }
        }

        this.mergedHeader.setSequenceDictionary(sequenceDictionary);

        // Set program that creates input alignments
        for (final SAMProgramRecord program : mergeProgramGroups(readers)) { 
            this.mergedHeader.addProgramRecord(program);
        }

        // Set read groups for merged header
        final List<SAMReadGroupRecord> readGroups = mergeReadGroups(readers);
        this.mergedHeader.setReadGroups(readGroups);
        this.mergedHeader.setGroupOrder(SAMFileHeader.GroupOrder.none);

        this.mergedHeader.setSortOrder(sortOrder);
    }

    /**
     * Checks to see if there are clashes where different readers are using the same read
     * group IDs. If they are then a new set of unique read group IDs are generated (across all
     * read groups) otherwise the original read group headers are returned.
     *
     * @param readers readers to combine
     * @return new list of readgroups constructed from all the readers
     */
    private List<SAMReadGroupRecord> mergeReadGroups(final Collection<SAMFileReader> readers) {
        final Map<String,SAMReadGroupRecord> seenGroups = new TreeMap<String,SAMReadGroupRecord>();

        // Check groups for duplicate entries 
        for (final SAMFileReader reader : readers) {
            for (final SAMReadGroupRecord group : reader.getFileHeader().getReadGroups()) {
                final String groupId = group.getReadGroupId() ;
                final SAMReadGroupRecord seenGroup = seenGroups.get(groupId);
                if (seenGroup == null) {
                    seenGroups.put(groupId, group);
                }
                else if (!group.equivalent(seenGroup)) { // same ID but different attributes
                    hasReadGroupCollisions = true;
                    break;
                }
            }
            if (hasReadGroupCollisions) break;
        }

        // If found no duplicates entries
        if (!hasReadGroupCollisions ) {
            return new ArrayList<SAMReadGroupRecord>(seenGroups.values());
        }

        // If found duplicates entries, renumber group IDs
        final List<SAMReadGroupRecord> newGroups = new ArrayList<SAMReadGroupRecord>();
        int idx=0;
        for (final SAMFileReader reader : readers) {
            final Map<String, String> idTranslation = new HashMap<String, String> ();
            samReadGroupIdTranslation.put(reader, idTranslation);
            for (final SAMReadGroupRecord group : reader.getFileHeader().getReadGroups()) {
                final String newGroupId = Integer.toString(idx);
                newGroups.add(new SAMReadGroupRecord(newGroupId, group));
                idTranslation.put(group.getReadGroupId(), newGroupId);
                idx++;
            }
        } 
        return newGroups ;
    }

    /**
     * Get the sequences off the SAMFileReader header.  Throws runtime exception if the sequence
     * are different from one another.
     *
     * @param readers readers to pull sequences from
     * @return sequences from files.  Each file should have the same sequence
     */
    private SAMSequenceDictionary getSequenceDictionary(final Collection<SAMFileReader> readers) {
        SAMSequenceDictionary sequences = null;
        for (final SAMFileReader reader : readers) {
            final SAMFileHeader header = reader.getFileHeader();

            if (sequences == null) {
                sequences = header.getSequenceDictionary();
            }
            else {
                final SAMSequenceDictionary currentSequences = header.getSequenceDictionary();
                SequenceUtil.assertSequenceDictionariesEqual(sequences, currentSequences);
            }
        }

        return sequences;
    }

    /**
     * Get the sequences from the SAMFileReader header, and merge the resulting sequence dictionaries.
     *
     * @param readers readers to pull sequences from
     * @return sequences from files.  Each file should have the same sequence
     */
    private SAMSequenceDictionary mergeSequenceDictionaries(final Collection<SAMFileReader> readers) {
        SAMSequenceDictionary sequences = new SAMSequenceDictionary();
        for (final SAMFileReader reader : readers) {
            final SAMSequenceDictionary currentSequences = reader.getFileHeader().getSequenceDictionary();
            sequences = mergeSequences(sequences, currentSequences);
        }
        // second pass, make a map of the original seqeunce id -> new sequence id
        createSequenceMapping(readers, sequences);
        return sequences;
    }

    /**
     * They've asked to merge the sequence headers.  What we support right now is finding the sequence name superset.
     *
     * @param currentDict the current dictionary, though merged entries are the superset of both dictionaries
     * @param mergingDict the sequence dictionary to merge
     * @return the superset dictionary, by sequence names
     */
    private SAMSequenceDictionary mergeSequences(SAMSequenceDictionary currentDict, SAMSequenceDictionary mergingDict) {
        LinkedList<SAMSequenceRecord> resultingDict = new LinkedList<SAMSequenceRecord>();
        LinkedList<String> resultingDictStr = new LinkedList<String>();

        // a place to hold the sequences that we haven't found a home for
        LinkedList<SAMSequenceRecord> holder = new LinkedList<SAMSequenceRecord>();

        resultingDict.addAll(currentDict.getSequences());
        for (SAMSequenceRecord r : resultingDict) {
            resultingDictStr.add(r.getSequenceName());
        }
        for (SAMSequenceRecord record : mergingDict.getSequences()) {
            if (resultingDictStr.contains(record.getSequenceName())) {
                int loc = resultingDictStr.indexOf(record.getSequenceName());
                resultingDict.addAll(loc, holder);
                holder.clear();
            } else {
                holder.add(record.clone());
            }
        }
        if (holder.size() != 0) {
            resultingDict.addAll(holder);
        }
        return new SAMSequenceDictionary(resultingDict);
    }


    /**
     * create the sequence mapping.  This map is used to convert the unmerged header sequence ID's to the merged
     * list of sequence id's.
     * @param readers the collections of readers.
     * @param masterDictionary the superset dictionary we've created.
     */
    private void createSequenceMapping(final Collection<SAMFileReader> readers, SAMSequenceDictionary masterDictionary) {
        LinkedList<String> resultingDictStr = new LinkedList<String>();
        for (SAMSequenceRecord r : masterDictionary.getSequences()) {
            resultingDictStr.add(r.getSequenceName());
        }
        for (final SAMFileReader reader : readers) {
            Map<Integer, Integer> seqMap = new HashMap<Integer, Integer>();
            SAMSequenceDictionary dict = reader.getFileHeader().getSequenceDictionary();
            for (SAMSequenceRecord rec : dict.getSequences()) {
                seqMap.put(rec.getSequenceIndex(), resultingDictStr.indexOf(rec.getSequenceName()));
            }
            this.samSeqDictionaryIdTranslation.put(reader, seqMap);
            this.samSeqDictionaryIdTranslationViaHeader.put(reader.getFileHeader(), seqMap);
        }
    }

    /**
     * Checks to see if there are clashes where different readers are using the same program
     * group IDs. If they are then a new set of unique program group IDs are generated (across all
     * program groups) otherwise the original program group headers are returned. 
     *
     * @param readers readers to combine
     * @return new list of program groups constructed from all the readers
     */

    private List<SAMProgramRecord> mergeProgramGroups(final Collection<SAMFileReader> readers) {
        final Map<String,SAMProgramRecord> seenGroups = new TreeMap<String,SAMProgramRecord>();

        // Check groups for duplicate entries 
        for (final SAMFileReader reader : readers) {
            for (final SAMProgramRecord group : reader.getFileHeader().getProgramRecords()) {
                final String groupId = group.getProgramGroupId() ;
                final SAMProgramRecord seenGroup = seenGroups.get(groupId);
                if (seenGroup == null) {
                    seenGroups.put(groupId, group);
                }
                else if (!group.equivalent(seenGroup)) { // same ID but different attributes
                    hasProgramGroupCollisions = true;
                    break;
                }
            }
            if (hasProgramGroupCollisions) break;
        }

        // If found no duplicates entries
        if (!hasProgramGroupCollisions ) {
            return new ArrayList<SAMProgramRecord>(seenGroups.values());
        }

        // If found duplicates entries, renumber group IDs
        final List<SAMProgramRecord> newGroups = new ArrayList<SAMProgramRecord>();
        int idx=0;
        for (final SAMFileReader reader : readers) {
            final Map<String, String> idTranslation = new TreeMap<String, String> ();
            samProgramGroupIdTranslation.put(reader, idTranslation);
            for (final SAMProgramRecord group : reader.getFileHeader().getProgramRecords()) {
                final String newGroupId = Integer.toString(idx);
                newGroups.add(new SAMProgramRecord(newGroupId, group));
                idTranslation.put(group.getProgramGroupId(), newGroupId);
                idx++;
            }
        } 
        return newGroups ;
    }
  
    /** Returns the read group id that should be used for the input read and RG id. */
    public String getReadGroupId(final SAMFileReader reader, final String originalReadGroupId) {
        return this.samReadGroupIdTranslation.get(reader).get(originalReadGroupId);
    }

    /**
     * @param reader one of the input files
     * @param originalProgramGroupId a program group ID from the above input file
     * @return new ID from the merged list of program groups in the output file
     */
    public String getProgramGroupId(final SAMFileReader reader, final String originalProgramGroupId) {
        return this.samProgramGroupIdTranslation.get(reader).get(originalProgramGroupId); 
        }

    /** Returns true if there are read group duplicates within the merged headers. */
    public boolean hasReadGroupCollisions() {
        return this.hasReadGroupCollisions;
    }

    /** Returns true if there are program group duplicates within the merged headers. */
    public boolean hasProgramGroupCollisions() {
        return hasProgramGroupCollisions;
    }

    /** @return if we've merged the sequence dictionaries, return true */
    public boolean hasMergedSequenceDictionary() {
        return hasMergedSequenceDictionary;
    }

    /** Returns the merged header that should be written to any output merged file. */
    public SAMFileHeader getMergedHeader() {
        return this.mergedHeader;
    }

    /** Returns the collection of readers that this header merger is working with. */
    public Collection<SAMFileReader> getReaders() {
        return this.readers;
    }

    /**
     * returns the new mapping for a specified reader, given it's old sequence index
     * @param reader the reader
     * @param oldReferenceSequenceIndex the old sequence (also called reference) index
     * @return the new index value
     */
    public Integer getMergedSequenceIndex(SAMFileReader reader, Integer oldReferenceSequenceIndex) {
        final Map<Integer, Integer> mapping = this.samSeqDictionaryIdTranslation.get(reader);
        if (mapping == null) {
            throw new PicardException("No sequence dictionary mapping available for reader: " + reader);
        }

        final Integer newIndex = mapping.get(oldReferenceSequenceIndex);
        if (newIndex == null) {
            throw new PicardException("No mapping for reference index " + oldReferenceSequenceIndex + " from reader: " + reader);
        }

        return newIndex;
    }

    /**
     * Another mechanism for getting the new sequence index, for situations in which the reader is not available.
     * Note that if the SAMRecord has already had its header replaced with the merged header, this won't work.
     * @param header The original header for the input record in question.
     * @param oldReferenceSequenceIndex The original sequence index.
     * @return the new index value that is compatible with the merged sequence index.
     */
    public Integer getMergedSequenceIndex(final SAMFileHeader header, Integer oldReferenceSequenceIndex) {
        final Map<Integer, Integer> mapping = this.samSeqDictionaryIdTranslationViaHeader.get(header);
        if (mapping == null) {
            throw new PicardException("No sequence dictionary mapping available for header: " + header);
        }

        final Integer newIndex = mapping.get(oldReferenceSequenceIndex);
        if (newIndex == null) {
            throw new PicardException("No mapping for reference index " + oldReferenceSequenceIndex + " from header: " + header);
        }

        return newIndex;
    }
}
