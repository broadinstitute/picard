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

import java.util.*;

import net.sf.picard.PicardException;
import net.sf.samtools.AbstractSAMHeaderRecord;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.SequenceUtil;

/**
 * Merges SAMFileHeaders that have the same sequences into a single merged header
 * object while providing read group translation for cases where read groups
 * clash across input headers.
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
    // so also store the mapping using theSAMFileHeader.
    // This is an IdentityHashMap because it can be quite expensive to compute the hashCode for
    // large SAMFileHeaders.  It is possible that two input files will have identical headers so that
    // the regular HashMap would fold them together, but the value stored in each of the two
    // Map entries will be the same, so it should not hurt anything.
    private final Map<SAMFileHeader,  Map<Integer, Integer>> samSeqDictionaryIdTranslationViaHeader =
            new IdentityHashMap<SAMFileHeader, Map<Integer, Integer>>();

    //HeaderRecordFactory that creates SAMReadGroupRecord instances.
    private static final HeaderRecordFactory<SAMReadGroupRecord> READ_GROUP_RECORD_FACTORY = new HeaderRecordFactory<SAMReadGroupRecord>() {
        public SAMReadGroupRecord createRecord(String id, SAMReadGroupRecord srcReadGroupRecord) {
            return new SAMReadGroupRecord(id, srcReadGroupRecord);
        }
    };

    //HeaderRecordFactory that creates SAMProgramRecord instances.
    private static final HeaderRecordFactory<SAMProgramRecord> PROGRAM_RECORD_FACTORY = new HeaderRecordFactory<SAMProgramRecord>() {
        public SAMProgramRecord createRecord(String id, SAMProgramRecord srcProgramRecord) {
            return new SAMProgramRecord(id, srcProgramRecord);
        }
    };

    //comparator used to sort lists of program group and read group records
    private static final Comparator<AbstractSAMHeaderRecord> RECORD_ID_COMPARATOR = new Comparator<AbstractSAMHeaderRecord>() {
        public int compare(AbstractSAMHeaderRecord o1, AbstractSAMHeaderRecord o2) {
            return o1.getId().compareTo(o2.getId());
        }
    };

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

        SAMSequenceDictionary sequenceDictionary;
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

        for (final SAMFileReader reader : readers) {
            for (final String comment : reader.getFileHeader().getComments()) {
                this.mergedHeader.addComment(comment);
            }
        }
    }

    /**
     * Checks to see if there are clashes where different readers are using the same read
     * group IDs. If yes, then those IDs that collided are remapped.
     *
     * @param readers readers to combine
     * @return new list of read groups constructed from all the readers
     */
    private List<SAMReadGroupRecord> mergeReadGroups(final Collection<SAMFileReader> readers) {
        //prepare args for mergeHeaderRecords(..) call
        final HashSet<String> idsThatAreAlreadyTaken = new HashSet<String>();

        final List<HeaderRecordAndFileReader<SAMReadGroupRecord>> readGroupsToProcess = new LinkedList<HeaderRecordAndFileReader<SAMReadGroupRecord>>();
        for (final SAMFileReader reader : readers) {
            for (final SAMReadGroupRecord readGroup : reader.getFileHeader().getReadGroups()) {
                //verify that there are no existing id collisions in this input file
                if(!idsThatAreAlreadyTaken.add(readGroup.getId()))
                    throw new PicardException("Input file: " + reader + " contains more than one RG with the same id (" + readGroup.getId() + ")");

                readGroupsToProcess.add(new HeaderRecordAndFileReader<SAMReadGroupRecord>(readGroup, reader));
            }
            idsThatAreAlreadyTaken.clear();
        }

        final List<SAMReadGroupRecord> result = new LinkedList<SAMReadGroupRecord>();

        hasReadGroupCollisions = mergeHeaderRecords(readGroupsToProcess, READ_GROUP_RECORD_FACTORY, idsThatAreAlreadyTaken, samReadGroupIdTranslation, result);

        //sort the result list by record id
        Collections.sort(result, RECORD_ID_COMPARATOR);

        return result;
    }


    /**
     * Checks to see if there are clashes where different readers are using the same program
     * group IDs. If yes, then those IDs that collided are remapped.
     *
     * @param readers readers to combine
     * @return new list of program groups constructed from all the readers
     */
    private List<SAMProgramRecord> mergeProgramGroups(final Collection<SAMFileReader> readers) {

        final List<SAMProgramRecord> overallResult = new LinkedList<SAMProgramRecord>();

        //this Set will accumulate all SAMProgramRecord ids that have been encountered so far.
        final HashSet<String> idsThatAreAlreadyTaken = new HashSet<String>();

        //need to process all program groups
        List<HeaderRecordAndFileReader<SAMProgramRecord>> programGroupsLeftToProcess = new LinkedList<HeaderRecordAndFileReader<SAMProgramRecord>>();
        for (final SAMFileReader reader : readers) {
            for (final SAMProgramRecord programGroup : reader.getFileHeader().getProgramRecords()) {
                //verify that there are no existing id collisions in this input file
                if(!idsThatAreAlreadyTaken.add(programGroup.getId()))
                    throw new PicardException("Input file: " + reader + " contains more than one PG with the same id (" + programGroup.getId() + ")");

                programGroupsLeftToProcess.add(new HeaderRecordAndFileReader<SAMProgramRecord>(programGroup, reader));
            }
            idsThatAreAlreadyTaken.clear();
        }

        //A program group header (lets say ID=2 PN=B PP=1) may have a PP (previous program) attribute which chains it to
        //another program group header (lets say ID=1 PN=A) to indicate that the given file was
        //processed by program A followed by program B. These PP attributes potentially
        //connect headers into one or more tree structures. Merging is done by
        //first merging all headers that don't have PP attributes (eg. tree roots),
        //then updating and merging all headers whose PPs point to the tree-root headers,
        //and so on until all program group headers are processed.

        //currentProgramGroups is the list of records to merge next. Start by merging the programGroups that don't have a PP attribute (eg. the tree roots).
        List< HeaderRecordAndFileReader<SAMProgramRecord> > currentProgramGroups = new LinkedList<HeaderRecordAndFileReader<SAMProgramRecord>>();
        for(final Iterator<HeaderRecordAndFileReader<SAMProgramRecord>> programGroupsLeftToProcessIterator = programGroupsLeftToProcess.iterator(); programGroupsLeftToProcessIterator.hasNext(); ) {
            final HeaderRecordAndFileReader<SAMProgramRecord> pair = programGroupsLeftToProcessIterator.next();
            if(pair.getHeaderRecord().getAttribute(SAMProgramRecord.PREVIOUS_PROGRAM_GROUP_ID_TAG) == null) {
                programGroupsLeftToProcessIterator.remove();
                currentProgramGroups.add(pair);
            }
        }

        //merge currentProgramGroups
        while(!currentProgramGroups.isEmpty())
        {
            final List<SAMProgramRecord> currentResult = new LinkedList<SAMProgramRecord>();

            hasProgramGroupCollisions |= mergeHeaderRecords(currentProgramGroups, PROGRAM_RECORD_FACTORY, idsThatAreAlreadyTaken, samProgramGroupIdTranslation, currentResult);

            //add currentResults to overallResults
            overallResult.addAll(currentResult);

            //apply the newly-computed id translations to currentProgramGroups and programGroupsLeftToProcess
            currentProgramGroups = translateIds(currentProgramGroups, samProgramGroupIdTranslation, false);
            programGroupsLeftToProcess = translateIds(programGroupsLeftToProcess, samProgramGroupIdTranslation, true);

            //find all records in programGroupsLeftToProcess whose ppId points to a record that was just processed (eg. a record that's in currentProgramGroups),
            //and move them to the list of programGroupsToProcessNext.
            LinkedList<HeaderRecordAndFileReader<SAMProgramRecord>> programGroupsToProcessNext = new LinkedList<HeaderRecordAndFileReader<SAMProgramRecord>>();
            for(final Iterator<HeaderRecordAndFileReader<SAMProgramRecord>> programGroupsLeftToProcessIterator = programGroupsLeftToProcess.iterator(); programGroupsLeftToProcessIterator.hasNext(); ) {
                final HeaderRecordAndFileReader<SAMProgramRecord> pairLeftToProcess = programGroupsLeftToProcessIterator.next();
                final Object ppIdOfRecordLeftToProcess = pairLeftToProcess.getHeaderRecord().getAttribute(SAMProgramRecord.PREVIOUS_PROGRAM_GROUP_ID_TAG);
                //find what currentProgramGroups this ppId points to (NOTE: they have to come from the same file)
                for(final HeaderRecordAndFileReader<SAMProgramRecord> justProcessedPair : currentProgramGroups) {
                    String idJustProcessed = justProcessedPair.getHeaderRecord().getId();
                    if(pairLeftToProcess.getFileReader() == justProcessedPair.getFileReader() && ppIdOfRecordLeftToProcess.equals(idJustProcessed)) {
                        programGroupsLeftToProcessIterator.remove();
                        programGroupsToProcessNext.add(pairLeftToProcess);
                        break;
                    }
                }
            }

            currentProgramGroups = programGroupsToProcessNext;
        }

        //verify that all records were processed
        if(!programGroupsLeftToProcess.isEmpty()) {
            StringBuffer errorMsg = new StringBuffer(programGroupsLeftToProcess.size() + " program groups weren't processed. Do their PP ids point to existing PGs? \n");
            for( final HeaderRecordAndFileReader<SAMProgramRecord> pair : programGroupsLeftToProcess ) {
                SAMProgramRecord record = pair.getHeaderRecord();
                errorMsg.append("@PG ID:"+record.getProgramGroupId()+" PN:"+record.getProgramName()+" PP:"+record.getPreviousProgramGroupId() +"\n");
            }
            throw new PicardException(errorMsg.toString());
        }

        //sort the result list by record id
        Collections.sort(overallResult, RECORD_ID_COMPARATOR);

        return overallResult;
    }


    /**
     * Utility method that takes a list of program groups and remaps all their
     * ids (including ppIds if requested) using the given idTranslationTable.
     *
     * NOTE: when remapping, this method creates new SAMProgramRecords and
     * doesn't mutate any records in the programGroups list.
     *
     * @param programGroups The program groups to translate.
     * @param idTranslationTable The translation table.
     * @param translatePpIds Whether ppIds should be translated as well.
     *
     * @return The list of translated records.
     */
    private List<HeaderRecordAndFileReader<SAMProgramRecord>> translateIds(
            List<HeaderRecordAndFileReader<SAMProgramRecord>> programGroups,
            Map<SAMFileReader, Map<String, String>> idTranslationTable,
            boolean translatePpIds) {

        //go through programGroups and translate any IDs and PPs based on the idTranslationTable.
        List<HeaderRecordAndFileReader<SAMProgramRecord>> result = new LinkedList<HeaderRecordAndFileReader<SAMProgramRecord>>();
        for(final HeaderRecordAndFileReader<SAMProgramRecord> pair : programGroups ) {
            final SAMProgramRecord record = pair.getHeaderRecord();
            final String id = record.getProgramGroupId();
            final String ppId = (String) record.getAttribute(SAMProgramRecord.PREVIOUS_PROGRAM_GROUP_ID_TAG);

            final SAMFileReader reader = pair.getFileReader();
            final Map<String, String> translations = idTranslationTable.get(reader);

            //see if one or both ids need to be translated
            SAMProgramRecord translatedRecord = null;
            if(translations != null)
            {
                String translatedId = translations.get( id );
                String translatedPpId = translatePpIds ? translations.get( ppId ) : null;

                boolean needToTranslateId = translatedId != null && !translatedId.equals(id);
                boolean needToTranslatePpId = translatedPpId != null && !translatedPpId.equals(ppId);

                if(needToTranslateId && needToTranslatePpId) {
                    translatedRecord = new SAMProgramRecord(translatedId, record);
                    translatedRecord.setAttribute(SAMProgramRecord.PREVIOUS_PROGRAM_GROUP_ID_TAG, translatedPpId);
                } else if(needToTranslateId) {
                    translatedRecord = new SAMProgramRecord(translatedId, record);
                } else if(needToTranslatePpId) {
                    translatedRecord = new SAMProgramRecord(id, record);
                    translatedRecord.setAttribute(SAMProgramRecord.PREVIOUS_PROGRAM_GROUP_ID_TAG, translatedPpId);
                }
            }

            if(translatedRecord != null) {
                result.add(new HeaderRecordAndFileReader<SAMProgramRecord>(translatedRecord, reader));
            } else {
                result.add(pair); //keep the original record
            }
        }

        return result;
    }


    /**
     * Utility method for merging a List of AbstractSAMHeaderRecords. If it finds
     * records that have identical ids and attributes, it will collapse them
     * into one record. If it finds records that have identical ids but
     * non-identical attributes, this is treated as a collision. When collision happens,
     * the records' ids are remapped, and an old-id to new-id mapping is added to the idTranslationTable.
     *
     * NOTE: Non-collided records also get recorded in the idTranslationTable as
     * old-id to old-id. This way, an idTranslationTable lookup should never return null.
     *
     * @param headerRecords The header records to merge.
     * @param headerRecordFactory Constructs a specific subclass of AbstractSAMHeaderRecord.
     * @param idsThatAreAlreadyTaken If the id of a headerRecord matches an id in this set, it will be treated as a collision, and the headRecord's id will be remapped.
     * @param idTranslationTable When records collide, their ids are remapped, and an old-id to new-id
     *      mapping is added to the idTranslationTable. Non-collided records also get recorded in the idTranslationTable as
     *      old-id to old-id. This way, an idTranslationTable lookup should never return null.
     *
     * @param result The list of merged header records.
     *
     * @return True if there were collisions.
     */
    private <RecordType extends AbstractSAMHeaderRecord> boolean mergeHeaderRecords(final List<HeaderRecordAndFileReader<RecordType>> headerRecords, HeaderRecordFactory<RecordType> headerRecordFactory,
            final HashSet<String> idsThatAreAlreadyTaken, Map<SAMFileReader, Map<String, String>> idTranslationTable, List<RecordType> result) {

        //The outer Map bins the header records by their ids. The nested Map further collapses
        //header records which, in addition to having the same id, also have identical attributes.
        //In other words, each key in the nested map represents one or more
        //header records which have both identical ids and identical attributes. The List of
        //SAMFileReaders keeps track of which readers these header record(s) came from.
        final Map<String, Map<RecordType, List<SAMFileReader>>> idToRecord =
            new HashMap<String, Map<RecordType, List<SAMFileReader>>>();

        //Populate the idToRecord and seenIds data structures
        for (final HeaderRecordAndFileReader<RecordType> pair : headerRecords) {
            final RecordType record = pair.getHeaderRecord();
            final SAMFileReader reader = pair.getFileReader();
            final String recordId = record.getId();
            Map<RecordType, List<SAMFileReader>> recordsWithSameId = idToRecord.get(recordId);
            if(recordsWithSameId == null) {
                recordsWithSameId = new LinkedHashMap<RecordType, List<SAMFileReader>>();
                idToRecord.put(recordId, recordsWithSameId);
            }

            List<SAMFileReader> fileReaders = recordsWithSameId.get(record);
            if(fileReaders == null) {
                fileReaders = new LinkedList<SAMFileReader>();
                recordsWithSameId.put(record, fileReaders);
            }

            fileReaders.add(reader);
        }

        //Resolve any collisions between header records by remapping their ids.
        boolean hasCollisions = false;
        for (final Map.Entry<String, Map<RecordType, List<SAMFileReader>>> entry : idToRecord.entrySet() )
        {
            final String recordId = entry.getKey();
            final Map<RecordType, List<SAMFileReader>> recordsWithSameId = entry.getValue();


            for( Map.Entry<RecordType, List<SAMFileReader>> recordWithUniqueAttr : recordsWithSameId.entrySet()) {
                final RecordType record = recordWithUniqueAttr.getKey();
                final List<SAMFileReader> fileReaders = recordWithUniqueAttr.getValue();

                String newId;
                if(!idsThatAreAlreadyTaken.contains(recordId)) {
                    //don't remap 1st record. If there are more records
                    //with this id, they will be remapped in the 'else'.
                    newId = recordId;
                    idsThatAreAlreadyTaken.add(recordId);
                } else {
                    //there is more than one record with this id.
                    hasCollisions = true;

                    //find a unique newId for this record
                    int idx=1;
                    while(idsThatAreAlreadyTaken.contains(newId = recordId + "." + Integer.toString(idx++)))
                        ;

                    idsThatAreAlreadyTaken.add( newId );
                }

                for(SAMFileReader fileReader : fileReaders) {
                    Map<String, String> readerTranslationTable = idTranslationTable.get(fileReader);
                    if(readerTranslationTable == null) {
                        readerTranslationTable = new HashMap<String, String>();
                        idTranslationTable.put(fileReader, readerTranslationTable);
                    }
                    readerTranslationTable.put(recordId, newId);
                }

                result.add( headerRecordFactory.createRecord(newId, record) );
            }
        }

        return hasCollisions;
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
     * @param mergeIntoDict the result of merging so far.  All SAMSequenceRecords in here have been cloned from the originals.
     * @param mergeFromDict A new sequence dictionary to merge into mergeIntoDict.
     * @return A new sequence dictionary that resulting from merging the two inputs.
     */
    private SAMSequenceDictionary mergeSequences(SAMSequenceDictionary mergeIntoDict, SAMSequenceDictionary mergeFromDict) {

        // a place to hold the sequences that we haven't found a home for, in the order the appear in mergeFromDict.
        LinkedList<SAMSequenceRecord> holder = new LinkedList<SAMSequenceRecord>();

        // Return value will be created from this.
        LinkedList<SAMSequenceRecord> resultingDict = new LinkedList<SAMSequenceRecord>();
        for (final SAMSequenceRecord sequenceRecord : mergeIntoDict.getSequences()) {
            resultingDict.add(sequenceRecord);
        }

        // Index into resultingDict of previous SAMSequenceRecord from mergeFromDict that already existed in mergeIntoDict.
        int prevloc = -1;
        // Previous SAMSequenceRecord from mergeFromDict that already existed in mergeIntoDict.
        SAMSequenceRecord previouslyMerged = null;

        for (SAMSequenceRecord sequenceRecord : mergeFromDict.getSequences()) {
            // Does it already exist in resultingDict?
            int loc = getIndexOfSequenceName(resultingDict, sequenceRecord.getSequenceName());
            if (loc == -1) {
                // If doesn't already exist in resultingDict, save it an decide where to insert it later.
                holder.add(sequenceRecord.clone());
            } else if (prevloc > loc) {
                // If sequenceRecord already exists in resultingDict, but prior to the previous one
                // from mergeIntoDict that already existed, cannot merge.
                throw new PicardException("Cannot merge sequence dictionaries because sequence " +
                        sequenceRecord.getSequenceName() + " and " + previouslyMerged.getSequenceName() +
                " are in different orders in two input sequence dictionaries.");
            } else {
                // Since sequenceRecord already exists in resultingDict, don't need to add it.
                // Add in all the sequences prior to it that have been held in holder.
                resultingDict.addAll(loc, holder);
                // Remember the index of sequenceRecord so can check for merge imcompatibility.
                prevloc = loc + holder.size();
                previouslyMerged = sequenceRecord;
                holder.clear();
            }
        }
        // Append anything left in holder.
        if (holder.size() != 0) {
            resultingDict.addAll(holder);
        }
        return new SAMSequenceDictionary(resultingDict);
    }

    /**
     * Find sequence in list.
     * @param list List to search for the sequence name.
     * @param sequenceName Name to search for.
     * @return Index of SAMSequenceRecord with the given name in list, or -1 if not found.
     */
    private static int getIndexOfSequenceName(final List<SAMSequenceRecord> list, final String sequenceName) {
        for (int i = 0; i < list.size(); ++i) {
            if (list.get(i).getSequenceName().equals(sequenceName)) {
                return i;
            }
        }
        return -1;
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


    /**
     * Implementations of this interface are used by mergeHeaderRecords(..) to instantiate
     * specific subclasses of AbstractSAMHeaderRecord.
     */
    private static interface HeaderRecordFactory<RecordType extends AbstractSAMHeaderRecord> {

       /**
        * Constructs a new instance of RecordType.
        * @param id The id of the new record.
        * @param srcRecord Except for the id, the new record will be a copy of this source record.
        */
        public RecordType createRecord(final String id, RecordType srcRecord);
    }

    /**
     * Struct that groups together a subclass of AbstractSAMHeaderRecord with the
     * SAMFileReader that it came from.
     */
    private static class HeaderRecordAndFileReader<RecordType extends AbstractSAMHeaderRecord> {
        private RecordType headerRecord;
        private SAMFileReader samFileReader;

        public HeaderRecordAndFileReader(RecordType headerRecord, SAMFileReader samFileReader) {
            this.headerRecord = headerRecord;
            this.samFileReader = samFileReader;
        }

        public RecordType getHeaderRecord() {
            return headerRecord;
        }
        public SAMFileReader getFileReader() {
            return samFileReader;
        }
    }
}
