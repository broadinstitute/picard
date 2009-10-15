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
    private final Map<SAMFileReader, Map<String, String>> samGroupIdTranslation =
            new HashMap<SAMFileReader, Map<String, String>>();

    //the groups from different files use the same group ids
    private boolean hasGroupIdDuplicates = false;

    //Translation of old program group ids to new program group ids
    private final Map<SAMFileReader, Map<String, String>> samProgramGroupIdTranslation =
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
        for (final SAMProgramRecord program : mergeSAMProgramRecordLists(readers)) {
            this.mergedHeader.addProgramRecord(program);
        }

        // Set read groups for merged header
        final List<SAMReadGroupRecord> readGroups = getReadGroups(readers);
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
    private List<SAMReadGroupRecord> getReadGroups(final Collection<SAMFileReader> readers) {
        // Read groups as read from the readers
        final List<SAMReadGroupRecord> orginalReadGroups = new ArrayList<SAMReadGroupRecord>();

        // Read group with new ids that don't confict
        final List<SAMReadGroupRecord> modifiedReadGroups = new ArrayList<SAMReadGroupRecord>();

        //set to see if there are duplicate group ids and whether or not we need to modify them
        final Set<String> groupIdsSeenBefore = new HashSet<String>();

        int x = 0;
        this.hasGroupIdDuplicates = false;

        for (final SAMFileReader reader : readers) {
            final SAMFileHeader header = reader.getFileHeader();
            final Map<String, String> idTranslation = new HashMap<String, String>();

            // Iterate over read groups to find conflicting ids
            for (final SAMReadGroupRecord readGroup : header.getReadGroups()) {
                final String groupId = readGroup.getReadGroupId();
                final String newGroupId = createNewId(x++);

                // Check to see if same group id is used in two different readers
                if (groupIdsSeenBefore.contains(groupId)) {
                    hasGroupIdDuplicates = true;
                }
                groupIdsSeenBefore.add(groupId);

                // Creates a new read group with the new id and copies all it's attributes
                final SAMReadGroupRecord groupRecordWithNewId = copyReadGroup(readGroup, newGroupId);

                orginalReadGroups.add(readGroup);
                modifiedReadGroups.add(groupRecordWithNewId);

                idTranslation.put(groupId, newGroupId);
            }

            // Add id tranlation for updating SamRecords with new ids if neccessary
            this.samGroupIdTranslation.put(reader, idTranslation);
        }

        // return approriate readgroups whether or not the new ids have to be used
        if (this.hasGroupIdDuplicates) {
            return modifiedReadGroups;
        }
        else {
            return orginalReadGroups;
        }
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
     * Find the alignment program that produced the readers.  If there are more than one
     * generate a new program represents that
     *
     * @param readers SAMFileReaders to pull program information from
     * @return SAMProgram record that represents all the readers
     */
    // TODO: this needs to be fixed up to support multiple program records (PIC-15)
    private List<SAMProgramRecord> mergeSAMProgramRecordLists(final Collection<SAMFileReader> readers) {
        final boolean programMixed = false;
        final List<SAMProgramRecord> ret = new ArrayList<SAMProgramRecord>();
        int nextProgramGroupId = 0;
        for (final SAMFileReader reader : readers) {
            final SAMFileHeader header = reader.getFileHeader();
            final Map<String, String> idTranslation = new HashMap<String, String>();
            for (final SAMProgramRecord oldProgramRecord : header.getProgramRecords()) {
                boolean foundMatch = false;
                for (final SAMProgramRecord newProgramRecord : ret) {
                    if (newProgramRecord.equivalent(oldProgramRecord)) {
                        idTranslation.put(oldProgramRecord.getProgramGroupId(), newProgramRecord.getProgramGroupId());
                        foundMatch = true;
                        break;
                    }
                }
                if (!foundMatch) {
                    final SAMProgramRecord newProgramRecord = new SAMProgramRecord(Integer.toString(nextProgramGroupId++));
                    copyProgramGroupAttributes(oldProgramRecord, newProgramRecord);
                    ret.add(newProgramRecord);
                    idTranslation.put(oldProgramRecord.getProgramGroupId(), newProgramRecord.getProgramGroupId());
                }
            }
            samProgramGroupIdTranslation.put(reader, idTranslation);
        }
        return ret;
    }

    private void copyProgramGroupAttributes(final SAMProgramRecord oldProgramRecord, final SAMProgramRecord newProgramRecord) {
        for (final Map.Entry<String, Object> entry : oldProgramRecord.getAttributes()) {
            newProgramRecord.setAttribute(entry.getKey(), entry.getValue());
        }
    }


    /**
     * Copies all the attribute of a readgroup to a new readgroup with a new id
     *
     * @param readGroup  the group to be copied
     * @param modifiedId the id for the new readgroup
     * @return new read group
     */
    private SAMReadGroupRecord copyReadGroup(final SAMReadGroupRecord readGroup, final String modifiedId) {
        final SAMReadGroupRecord retval = new SAMReadGroupRecord(modifiedId);
        retval.setLibrary(readGroup.getLibrary());
        retval.setSample(readGroup.getSample());

        for (final Map.Entry<String, Object> attr : readGroup.getAttributes()) {
            retval.setAttribute(attr.getKey(), attr.getValue());
        }

        return retval;
    }


    /**
     * Creates a base 26 representation of an int
     *
     * @param n int to covert to letter representation
     * @return string rep for an int eg 0 = A  27 = AB
     */
    protected static String createNewId(int n) {
        final int base = ALPHABET.length();

        String s = "";
        while (true) {
            final int r = n % base;
            s = ALPHABET.charAt(r) + s;
            n = n / base;
            if (n == 0) {
                return s;
            }
            n -= 1;
        }
    }

    /** Returns the read group id that should be used for the input read and RG id. */
    public String getReadGroupId(final SAMFileReader reader, final String originalReadGroupId) {
        return this.samGroupIdTranslation.get(reader).get(originalReadGroupId);
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
    public boolean hasGroupIdDuplicates() {
        return this.hasGroupIdDuplicates;
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
