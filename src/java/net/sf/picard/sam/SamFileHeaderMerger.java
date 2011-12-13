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
    private Collection<SAMFileReader> readers;
    private final Collection<SAMFileHeader> headers;

    //Translation of old group ids to new group ids
    private final Map<SAMFileHeader, Map<String, String>> samReadGroupIdTranslation =
            new IdentityHashMap<SAMFileHeader, Map<String, String>>();

    //the read groups from different files use the same group ids
    private boolean hasReadGroupCollisions = false;

    //the program records from different files use the same program record ids
    private boolean hasProgramGroupCollisions = false;

    //Translation of old program group ids to new program group ids
    private final Map<SAMFileHeader, Map<String, String>> samProgramGroupIdTranslation =
            new IdentityHashMap<SAMFileHeader, Map<String, String>>();

    private boolean hasMergedSequenceDictionary = false;

    // Translation of old sequence dictionary ids to new dictionary ids
    // This is an IdentityHashMap because it can be quite expensive to compute the hashCode for
    // large SAMFileHeaders.  It is possible that two input files will have identical headers so that
    // the regular HashMap would fold them together, but the value stored in each of the two
    // Map entries will be the same, so it should not hurt anything.
    private final Map<SAMFileHeader,  Map<Integer, Integer>> samSeqDictionaryIdTranslationViaHeader =
            new IdentityHashMap<SAMFileHeader, Map<Integer, Integer>>();

    //HeaderRecordFactory that creates SAMReadGroupRecord instances.
    private static final HeaderRecordFactory<SAMReadGroupRecord> READ_GROUP_RECORD_FACTORY = new HeaderRecordFactory<SAMReadGroupRecord>() {
        public SAMReadGroupRecord createRecord(final String id, final SAMReadGroupRecord srcReadGroupRecord) {
            return new SAMReadGroupRecord(id, srcReadGroupRecord);
        }
    };

    //HeaderRecordFactory that creates SAMProgramRecord instances.
    private static final HeaderRecordFactory<SAMProgramRecord> PROGRAM_RECORD_FACTORY = new HeaderRecordFactory<SAMProgramRecord>() {
        public SAMProgramRecord createRecord(final String id, final SAMProgramRecord srcProgramRecord) {
            return new SAMProgramRecord(id, srcProgramRecord);
        }
    };

    //comparator used to sort lists of program group and read group records
    private static final Comparator<AbstractSAMHeaderRecord> RECORD_ID_COMPARATOR = new Comparator<AbstractSAMHeaderRecord>() {
        public int compare(final AbstractSAMHeaderRecord o1, final AbstractSAMHeaderRecord o2) {
            return o1.getId().compareTo(o2.getId());
        }
    };

    /**
     * Create SAMFileHeader with additional information.  Required that sequence dictionaries agree.
     *
     * @param readers sam file readers to combine
     * @param sortOrder sort order new header should have
     * @deprecated replaced by SamFileHeaderMerger(Collection<SAMFileHeader>, SAMFileHeader.SortOrder, boolean)
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
     * @deprecated replaced by SamFileHeaderMerger(Collection<SAMFileHeader>, SAMFileHeader.SortOrder, boolean)
     */
    public SamFileHeaderMerger(final Collection<SAMFileReader> readers, final SAMFileHeader.SortOrder sortOrder, final boolean mergeDictionaries) {
        this(sortOrder, getHeadersFromReaders(readers), mergeDictionaries);
        this.readers = readers;
    }

    /**
     * Create SAMFileHeader with additional information..  This is the preferred constructor.
     *
     * @param sortOrder sort order new header should have
     * @param headers sam file headers to combine
     * @param mergeDictionaries If true, merge sequence dictionaries in new header.  If false, require that
     * all input sequence dictionaries be identical.
     */
    public SamFileHeaderMerger(final SAMFileHeader.SortOrder sortOrder, final Collection<SAMFileHeader> headers, final boolean mergeDictionaries) {
        this.headers = headers;
        this.mergedHeader = new SAMFileHeader();

        SAMSequenceDictionary sequenceDictionary;
        try {
            sequenceDictionary = getSequenceDictionary(headers);
            this.hasMergedSequenceDictionary = false;
        }
        catch (SequenceUtil.SequenceListsDifferException pe) {
            if (mergeDictionaries) {
                sequenceDictionary = mergeSequenceDictionaries(headers);
                this.hasMergedSequenceDictionary = true;
            }
            else {
                throw pe;
            }
        }

        this.mergedHeader.setSequenceDictionary(sequenceDictionary);

        // Set program that creates input alignments
        for (final SAMProgramRecord program : mergeProgramGroups(headers)) {
            this.mergedHeader.addProgramRecord(program);
        }

        // Set read groups for merged header
        final List<SAMReadGroupRecord> readGroups = mergeReadGroups(headers);
        this.mergedHeader.setReadGroups(readGroups);
        this.mergedHeader.setGroupOrder(SAMFileHeader.GroupOrder.none);

        this.mergedHeader.setSortOrder(sortOrder);

        for (final SAMFileHeader header : headers) {
            for (final String comment : header.getComments()) {
                this.mergedHeader.addComment(comment);
            }
        }
    }

    // Utilility method to make use with old constructor
    private static List<SAMFileHeader> getHeadersFromReaders(final Collection<SAMFileReader> readers) {
        final List<SAMFileHeader> headers = new ArrayList<SAMFileHeader>(readers.size());
        for (final SAMFileReader reader : readers) {
            headers.add(reader.getFileHeader());
        }
        return headers;
    }


    /**
     * Checks to see if there are clashes where different readers are using the same read
     * group IDs. If yes, then those IDs that collided are remapped.
     *
     * @param headers headers to combine
     * @return new list of read groups constructed from all the readers
     */
    private List<SAMReadGroupRecord> mergeReadGroups(final Collection<SAMFileHeader> headers) {
        //prepare args for mergeHeaderRecords(..) call
        final HashSet<String> idsThatAreAlreadyTaken = new HashSet<String>();

        final List<HeaderRecordAndFileHeader<SAMReadGroupRecord>> readGroupsToProcess = new LinkedList<HeaderRecordAndFileHeader<SAMReadGroupRecord>>();
        for (final SAMFileHeader header : headers) {
            for (final SAMReadGroupRecord readGroup : header.getReadGroups()) {
                //verify that there are no existing id collisions in this input file
                if(!idsThatAreAlreadyTaken.add(readGroup.getId()))
                    throw new PicardException("Input file: " + header + " contains more than one RG with the same id (" + readGroup.getId() + ")");

                readGroupsToProcess.add(new HeaderRecordAndFileHeader<SAMReadGroupRecord>(readGroup, header));
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
     * @param headers headers to combine
     * @return new list of program groups constructed from all the readers
     */
    private List<SAMProgramRecord> mergeProgramGroups(final Collection<SAMFileHeader> headers) {

        final List<SAMProgramRecord> overallResult = new LinkedList<SAMProgramRecord>();

        //this Set will accumulate all SAMProgramRecord ids that have been encountered so far.
        final HashSet<String> idsThatAreAlreadyTaken = new HashSet<String>();

        //need to process all program groups
        List<HeaderRecordAndFileHeader<SAMProgramRecord>> programGroupsLeftToProcess = new LinkedList<HeaderRecordAndFileHeader<SAMProgramRecord>>();
        for (final SAMFileHeader header : headers) {
            for (final SAMProgramRecord programGroup : header.getProgramRecords()) {
                //verify that there are no existing id collisions in this input file
                if(!idsThatAreAlreadyTaken.add(programGroup.getId()))
                    throw new PicardException("Input file: " + header + " contains more than one PG with the same id (" + programGroup.getId() + ")");

                programGroupsLeftToProcess.add(new HeaderRecordAndFileHeader<SAMProgramRecord>(programGroup, header));
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
        List< HeaderRecordAndFileHeader<SAMProgramRecord> > currentProgramGroups = new LinkedList<HeaderRecordAndFileHeader<SAMProgramRecord>>();
        for(final Iterator<HeaderRecordAndFileHeader<SAMProgramRecord>> programGroupsLeftToProcessIterator = programGroupsLeftToProcess.iterator(); programGroupsLeftToProcessIterator.hasNext(); ) {
            final HeaderRecordAndFileHeader<SAMProgramRecord> pair = programGroupsLeftToProcessIterator.next();
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
            final LinkedList<HeaderRecordAndFileHeader<SAMProgramRecord>> programGroupsToProcessNext = new LinkedList<HeaderRecordAndFileHeader<SAMProgramRecord>>();
            for(final Iterator<HeaderRecordAndFileHeader<SAMProgramRecord>> programGroupsLeftToProcessIterator = programGroupsLeftToProcess.iterator(); programGroupsLeftToProcessIterator.hasNext(); ) {
                final HeaderRecordAndFileHeader<SAMProgramRecord> pairLeftToProcess = programGroupsLeftToProcessIterator.next();
                final Object ppIdOfRecordLeftToProcess = pairLeftToProcess.getHeaderRecord().getAttribute(SAMProgramRecord.PREVIOUS_PROGRAM_GROUP_ID_TAG);
                //find what currentProgramGroups this ppId points to (NOTE: they have to come from the same file)
                for(final HeaderRecordAndFileHeader<SAMProgramRecord> justProcessedPair : currentProgramGroups) {
                    final String idJustProcessed = justProcessedPair.getHeaderRecord().getId();
                    if(pairLeftToProcess.getFileHeader() == justProcessedPair.getFileHeader() && ppIdOfRecordLeftToProcess.equals(idJustProcessed)) {
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
            final StringBuffer errorMsg = new StringBuffer(programGroupsLeftToProcess.size() + " program groups weren't processed. Do their PP ids point to existing PGs? \n");
            for( final HeaderRecordAndFileHeader<SAMProgramRecord> pair : programGroupsLeftToProcess ) {
                final SAMProgramRecord record = pair.getHeaderRecord();
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
    private List<HeaderRecordAndFileHeader<SAMProgramRecord>> translateIds(
            final List<HeaderRecordAndFileHeader<SAMProgramRecord>> programGroups,
            final Map<SAMFileHeader, Map<String, String>> idTranslationTable,
            final boolean translatePpIds) {

        //go through programGroups and translate any IDs and PPs based on the idTranslationTable.
        final List<HeaderRecordAndFileHeader<SAMProgramRecord>> result = new LinkedList<HeaderRecordAndFileHeader<SAMProgramRecord>>();
        for(final HeaderRecordAndFileHeader<SAMProgramRecord> pair : programGroups ) {
            final SAMProgramRecord record = pair.getHeaderRecord();
            final String id = record.getProgramGroupId();
            final String ppId = (String) record.getAttribute(SAMProgramRecord.PREVIOUS_PROGRAM_GROUP_ID_TAG);

            final SAMFileHeader header = pair.getFileHeader();
            final Map<String, String> translations = idTranslationTable.get(header);

            //see if one or both ids need to be translated
            SAMProgramRecord translatedRecord = null;
            if(translations != null)
            {
                final String translatedId = translations.get( id );
                final String translatedPpId = translatePpIds ? translations.get( ppId ) : null;

                final boolean needToTranslateId = translatedId != null && !translatedId.equals(id);
                final boolean needToTranslatePpId = translatedPpId != null && !translatedPpId.equals(ppId);

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
                result.add(new HeaderRecordAndFileHeader<SAMProgramRecord>(translatedRecord, header));
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
    private <RecordType extends AbstractSAMHeaderRecord> boolean mergeHeaderRecords(final List<HeaderRecordAndFileHeader<RecordType>> headerRecords, final HeaderRecordFactory<RecordType> headerRecordFactory,
            final HashSet<String> idsThatAreAlreadyTaken, final Map<SAMFileHeader, Map<String, String>> idTranslationTable, final List<RecordType> result) {

        //The outer Map bins the header records by their ids. The nested Map further collapses
        //header records which, in addition to having the same id, also have identical attributes.
        //In other words, each key in the nested map represents one or more
        //header records which have both identical ids and identical attributes. The List of
        //SAMFileHeaders keeps track of which readers these header record(s) came from.
        final Map<String, Map<RecordType, List<SAMFileHeader>>> idToRecord =
            new HashMap<String, Map<RecordType, List<SAMFileHeader>>>();

        //Populate the idToRecord and seenIds data structures
        for (final HeaderRecordAndFileHeader<RecordType> pair : headerRecords) {
            final RecordType record = pair.getHeaderRecord();
            final SAMFileHeader header = pair.getFileHeader();
            final String recordId = record.getId();
            Map<RecordType, List<SAMFileHeader>> recordsWithSameId = idToRecord.get(recordId);
            if(recordsWithSameId == null) {
                recordsWithSameId = new LinkedHashMap<RecordType, List<SAMFileHeader>>();
                idToRecord.put(recordId, recordsWithSameId);
            }

            List<SAMFileHeader> fileHeaders = recordsWithSameId.get(record);
            if(fileHeaders == null) {
                fileHeaders = new LinkedList<SAMFileHeader>();
                recordsWithSameId.put(record, fileHeaders);
            }

            fileHeaders.add(header);
        }

        //Resolve any collisions between header records by remapping their ids.
        boolean hasCollisions = false;
        for (final Map.Entry<String, Map<RecordType, List<SAMFileHeader>>> entry : idToRecord.entrySet() )
        {
            final String recordId = entry.getKey();
            final Map<RecordType, List<SAMFileHeader>> recordsWithSameId = entry.getValue();


            for( final Map.Entry<RecordType, List<SAMFileHeader>> recordWithUniqueAttr : recordsWithSameId.entrySet()) {
                final RecordType record = recordWithUniqueAttr.getKey();
                final List<SAMFileHeader> fileHeaders = recordWithUniqueAttr.getValue();

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

                for(final SAMFileHeader fileHeader : fileHeaders) {
                    Map<String, String> readerTranslationTable = idTranslationTable.get(fileHeader);
                    if(readerTranslationTable == null) {
                        readerTranslationTable = new HashMap<String, String>();
                        idTranslationTable.put(fileHeader, readerTranslationTable);
                    }
                    readerTranslationTable.put(recordId, newId);
                }

                result.add( headerRecordFactory.createRecord(newId, record) );
            }
        }

        return hasCollisions;
    }


    /**
     * Get the sequences off the SAMFileHeader.  Throws runtime exception if the sequence
     * are different from one another.
     *
     * @param headers headers to pull sequences from
     * @return sequences from files.  Each file should have the same sequence
     */
    private SAMSequenceDictionary getSequenceDictionary(final Collection<SAMFileHeader> headers) {
        SAMSequenceDictionary sequences = null;
        for (final SAMFileHeader header : headers) {

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
     * Get the sequences from the SAMFileHeader, and merge the resulting sequence dictionaries.
     *
     * @param headers headers to pull sequences from
     * @return sequences from files.  Each file should have the same sequence
     */
    private SAMSequenceDictionary mergeSequenceDictionaries(final Collection<SAMFileHeader> headers) {
        SAMSequenceDictionary sequences = new SAMSequenceDictionary();
        for (final SAMFileHeader header : headers) {
            final SAMSequenceDictionary currentSequences = header.getSequenceDictionary();
            sequences = mergeSequences(sequences, currentSequences);
        }
        // second pass, make a map of the original seqeunce id -> new sequence id
        createSequenceMapping(headers, sequences);
        return sequences;
    }

    /**
     * They've asked to merge the sequence headers.  What we support right now is finding the sequence name superset.
     *
     * @param mergeIntoDict the result of merging so far.  All SAMSequenceRecords in here have been cloned from the originals.
     * @param mergeFromDict A new sequence dictionary to merge into mergeIntoDict.
     * @return A new sequence dictionary that resulting from merging the two inputs.
     */
    private SAMSequenceDictionary mergeSequences(final SAMSequenceDictionary mergeIntoDict, final SAMSequenceDictionary mergeFromDict) {

        // a place to hold the sequences that we haven't found a home for, in the order the appear in mergeFromDict.
        final LinkedList<SAMSequenceRecord> holder = new LinkedList<SAMSequenceRecord>();

        // Return value will be created from this.
        final LinkedList<SAMSequenceRecord> resultingDict = new LinkedList<SAMSequenceRecord>();
        for (final SAMSequenceRecord sequenceRecord : mergeIntoDict.getSequences()) {
            resultingDict.add(sequenceRecord);
        }

        // Index into resultingDict of previous SAMSequenceRecord from mergeFromDict that already existed in mergeIntoDict.
        int prevloc = -1;
        // Previous SAMSequenceRecord from mergeFromDict that already existed in mergeIntoDict.
        SAMSequenceRecord previouslyMerged = null;

        for (final SAMSequenceRecord sequenceRecord : mergeFromDict.getSequences()) {
            // Does it already exist in resultingDict?
            final int loc = getIndexOfSequenceName(resultingDict, sequenceRecord.getSequenceName());
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
     * @param headers the collections of headers.
     * @param masterDictionary the superset dictionary we've created.
     */
    private void createSequenceMapping(final Collection<SAMFileHeader> headers, final SAMSequenceDictionary masterDictionary) {
        final LinkedList<String> resultingDictStr = new LinkedList<String>();
        for (final SAMSequenceRecord r : masterDictionary.getSequences()) {
            resultingDictStr.add(r.getSequenceName());
        }
        for (final SAMFileHeader header : headers) {
            final Map<Integer, Integer> seqMap = new HashMap<Integer, Integer>();
            final SAMSequenceDictionary dict = header.getSequenceDictionary();
            for (final SAMSequenceRecord rec : dict.getSequences()) {
                seqMap.put(rec.getSequenceIndex(), resultingDictStr.indexOf(rec.getSequenceName()));
            }
            this.samSeqDictionaryIdTranslationViaHeader.put(header, seqMap);
        }
    }



    /**
     * Returns the read group id that should be used for the input read and RG id.
     *
     * @deprecated replaced by getReadGroupId(SAMFileHeader, String)
     * */
    public String getReadGroupId(final SAMFileReader reader, final String originalReadGroupId) {
        return getReadGroupId(reader.getFileHeader(), originalReadGroupId);
    }

    /** Returns the read group id that should be used for the input read and RG id. */
    public String getReadGroupId(final SAMFileHeader header, final String originalReadGroupId) {
        return this.samReadGroupIdTranslation.get(header).get(originalReadGroupId);
    }

    /**
     * @param reader one of the input files
     * @param originalProgramGroupId a program group ID from the above input file
     * @return new ID from the merged list of program groups in the output file
     * @deprecated replaced by getProgramGroupId(SAMFileHeader, String)
     */
    public String getProgramGroupId(final SAMFileReader reader, final String originalProgramGroupId) {
        return getProgramGroupId(reader.getFileHeader(), originalProgramGroupId);
    }

    /**
     * @param header one of the input headers
     * @param originalProgramGroupId a program group ID from the above input file
     * @return new ID from the merged list of program groups in the output file
     */
    public String getProgramGroupId(final SAMFileHeader header, final String originalProgramGroupId) {
        return this.samProgramGroupIdTranslation.get(header).get(originalProgramGroupId);
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

    /** Returns the collection of readers that this header merger is working with. May return null.
     * @deprecated replaced by getHeaders()
     */
    public Collection<SAMFileReader> getReaders() {
        return this.readers;
    }

    /** Returns the collection of readers that this header merger is working with.
     */
    public Collection<SAMFileHeader> getHeaders() {
        return this.headers;
    }

    /**
     * returns the new mapping for a specified reader, given it's old sequence index
     * @param reader the reader
     * @param oldReferenceSequenceIndex the old sequence (also called reference) index
     * @return the new index value
     * @deprecated replaced by getMergedSequenceIndex(SAMFileHeader, Integer)
     */
    public Integer getMergedSequenceIndex(final SAMFileReader reader, final Integer oldReferenceSequenceIndex) {
        return this.getMergedSequenceIndex(reader.getFileHeader(), oldReferenceSequenceIndex);
    }

    /**
     * Another mechanism for getting the new sequence index, for situations in which the reader is not available.
     * Note that if the SAMRecord has already had its header replaced with the merged header, this won't work.
     * @param header The original header for the input record in question.
     * @param oldReferenceSequenceIndex The original sequence index.
     * @return the new index value that is compatible with the merged sequence index.
     */
    public Integer getMergedSequenceIndex(final SAMFileHeader header, final Integer oldReferenceSequenceIndex) {
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
     * SAMFileHeader that it came from.
     */
    private static class HeaderRecordAndFileHeader<RecordType extends AbstractSAMHeaderRecord> {
        private final RecordType headerRecord;
        private final SAMFileHeader samFileHeader;

        public HeaderRecordAndFileHeader(final RecordType headerRecord, final SAMFileHeader samFileHeader) {
            this.headerRecord = headerRecord;
            this.samFileHeader = samFileHeader;
        }

        public RecordType getHeaderRecord() {
            return headerRecord;
        }
        public SAMFileHeader getFileHeader() {
            return samFileHeader;
        }
    }
}
