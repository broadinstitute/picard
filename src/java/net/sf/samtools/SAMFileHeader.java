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
package net.sf.samtools;


import net.sf.samtools.util.StringLineReader;

import java.lang.reflect.Constructor;
import java.util.*;
import java.io.StringWriter;

/**
 * Header information from a SAM or BAM file.
 */
public class SAMFileHeader extends AbstractSAMHeaderRecord
{
    public static final String VERSION_TAG = "VN";
    public static final String SORT_ORDER_TAG = "SO";
    public static final String GROUP_ORDER_TAG = "GO";
    public static final String CURRENT_VERSION = "1.0";
    private static final String NEXT_VERSION = "1.3";
    public static final Set<String> ACCEPTABLE_VERSIONS = new HashSet<String>(Arrays.asList(CURRENT_VERSION, NEXT_VERSION));

    /**
     * These tags are of known type, so don't need a type field in the text representation.
     */
    public static final Set<String> STANDARD_TAGS =
            new HashSet<String>(Arrays.asList(VERSION_TAG, SORT_ORDER_TAG, GROUP_ORDER_TAG));

    Set<String> getStandardTags() {
        return STANDARD_TAGS;
    }

    /**
     * Ways in which a SAM or BAM may be sorted.
     */
    public enum SortOrder {

        unsorted(null),
        queryname(SAMRecordQueryNameComparator.class),
        coordinate(SAMRecordCoordinateComparator.class);

        private final Class<? extends SAMRecordComparator> comparator;

        SortOrder(final Class<? extends SAMRecordComparator> comparatorClass) {
            this.comparator = comparatorClass;
        }

        /**
         * @return Comparator class to sort in the specified order, or null if unsorted.
         */
        public Class<? extends SAMRecordComparator> getComparator() {
            return comparator;
        }

        /**
         * @return Comparator to sort in the specified order, or null if unsorted.
         */
        public SAMRecordComparator getComparatorInstance() {
            if (comparator != null) {
                try {
                    final Constructor<? extends SAMRecordComparator> ctor = comparator.getConstructor();
                    return ctor.newInstance();
                }
                catch (Exception e) {
                    throw new IllegalStateException("Could not instantiate a comparator for sort order: " +
                            this.name(), e);
                }
            }
            return null;
        }
    }

    public enum GroupOrder {
        none, query, reference
    }

    private List<SAMReadGroupRecord> mReadGroups =
        new ArrayList<SAMReadGroupRecord>();
    private List<SAMProgramRecord> mProgramRecords = new ArrayList<SAMProgramRecord>();
    private final Map<String, SAMReadGroupRecord> mReadGroupMap =
        new HashMap<String, SAMReadGroupRecord>();
    private final Map<String, SAMProgramRecord> mProgramRecordMap = new HashMap<String, SAMProgramRecord>();
    private SAMSequenceDictionary mSequenceDictionary = new SAMSequenceDictionary();
    final private List<String> mComments = new ArrayList<String>();
    private String textHeader;
    private final List<SAMValidationError> mValidationErrors = new ArrayList<SAMValidationError>();

    public SAMFileHeader() {
        setAttribute(VERSION_TAG, CURRENT_VERSION);
    }

    public String getVersion() {
        return (String) getAttribute("VN");
    }

    public String getCreator() {
        return (String) getAttribute("CR");
    }

    public SAMSequenceDictionary getSequenceDictionary() {
        return mSequenceDictionary;
    }

    public List<SAMReadGroupRecord> getReadGroups() {
        return Collections.unmodifiableList(mReadGroups);
    }

    /**
     * Look up sequence record by name.
     */
    public SAMSequenceRecord getSequence(final String name) {
        return mSequenceDictionary.getSequence(name);
    }

    /**
     * Look up read group record by name.
     */
    public SAMReadGroupRecord getReadGroup(final String name) {
        return mReadGroupMap.get(name);
    }

    /**
     * Replace entire sequence dictionary.  The given sequence dictionary is stored, not copied.
     */
    public void setSequenceDictionary(final SAMSequenceDictionary sequenceDictionary) {
        mSequenceDictionary = sequenceDictionary;
    }

    public void addSequence(final SAMSequenceRecord sequenceRecord) {
        mSequenceDictionary.addSequence(sequenceRecord);
    }

    /**
     * Look up a sequence record by index.  First sequence in the header is the 0th.
     * @return The corresponding sequence record, or null if the index is out of range.
     */
    public SAMSequenceRecord getSequence(final int sequenceIndex) {
        return mSequenceDictionary.getSequence(sequenceIndex);
    }

    /**
     *
     * @return Sequence index for the given sequence name, or -1 if the name is not found.
     */
    public int getSequenceIndex(final String sequenceName) {
        return mSequenceDictionary.getSequenceIndex(sequenceName);
    }

    /**
     * Replace entire list of read groups.  The given list is stored, not copied.
     */
    public void setReadGroups(final List<SAMReadGroupRecord> readGroups) {
        mReadGroups = readGroups;
        mReadGroupMap.clear();
        for (final SAMReadGroupRecord readGroupRecord : readGroups) {
            mReadGroupMap.put(readGroupRecord.getReadGroupId(), readGroupRecord);
        }
    }

    public void addReadGroup(final SAMReadGroupRecord readGroup) {
        mReadGroups.add(readGroup);
        mReadGroupMap.put(readGroup.getReadGroupId(), readGroup);
    }

    public List<SAMProgramRecord> getProgramRecords() {
        return Collections.unmodifiableList(mProgramRecords);
    }

    public void addProgramRecord(final SAMProgramRecord programRecord) {
        this.mProgramRecords.add(programRecord);
        this.mProgramRecordMap.put(programRecord.getProgramGroupId(), programRecord);
    }

    public SAMProgramRecord getProgramRecord(final String name) {
        return this.mProgramRecordMap.get(name);
    }

    /**
     * Replace entire list of program records
     * @param programRecords This list is used directly, not copied.
     */
    public void setProgramRecords(final List<SAMProgramRecord> programRecords) {
        this.mProgramRecords = programRecords;
        this.mProgramRecordMap.clear();
        for (final SAMProgramRecord programRecord : this.mProgramRecords) {
            this.mProgramRecordMap.put(programRecord.getProgramGroupId(), programRecord);
        }
    }

    /**
     * @return a new SAMProgramRecord with an ID guaranteed to not exist in this SAMFileHeader
     */
    public SAMProgramRecord createProgramRecord() {
        for (int i = 0; i < Integer.MAX_VALUE; ++i) {
            final String s = Integer.toString(i);
            if (!this.mProgramRecordMap.containsKey(s)) {
                final SAMProgramRecord ret = new SAMProgramRecord(s);
                addProgramRecord(ret);
                return ret;
            }
        }
        throw new IllegalStateException("Surprising number of SAMProgramRecords");
    }

    public SortOrder getSortOrder() {
        if (getAttribute("SO") == null) {
            return SortOrder.unsorted;
        }
        return SortOrder.valueOf((String)getAttribute("SO"));
    }

    public void setSortOrder(final SortOrder so) {
        setAttribute("SO", so.name());
    }

    public GroupOrder getGroupOrder() {
        if (getAttribute("GO") == null) {
            return GroupOrder.none;
        }
        return GroupOrder.valueOf((String)getAttribute("GO"));
    }

    public void setGroupOrder(final GroupOrder go) {
        setAttribute("GO", go.name());
    }

    /**
     * If this SAMHeader was read from a file, this property contains the header
     * as it appeared in the file, otherwise it is null.  Note that this is not a toString()
     * operation.  Changes to the SAMFileHeader object after reading from the file are not reflected in this value.
     *
     * In addition this value is only set if one of the following is true:
     *   - The size of the header is < 1,048,576 characters (1MB ascii, 2MB unicode)
     *   - There are either validation or parsing errors associated with the header
     *
     * Invalid header lines may appear in value but are not stored in the SAMFileHeader object.
     */
    public String getTextHeader() {
        return textHeader;
    }

    public void setTextHeader(final String textHeader) {
        this.textHeader = textHeader;
    }

    public List<String> getComments() {
        return Collections.unmodifiableList(mComments);
    }

    public void addComment(String comment) {
        if (!comment.startsWith(SAMTextHeaderCodec.COMMENT_PREFIX)) {
            comment = SAMTextHeaderCodec.COMMENT_PREFIX + comment;
        }
        mComments.add(comment);
    }


    /**
     * Replace existing comments with the contents of the given collection.
     */
    public void setComments(final Collection<String> comments) {
        mComments.clear();
        for (final String comment : comments) {
            addComment(comment);
        }
    }

    public List<SAMValidationError> getValidationErrors() {
        return Collections.unmodifiableList(mValidationErrors);
    }

    public void addValidationError(final SAMValidationError error) {
        mValidationErrors.add(error);
    }

    /**
     * Replace list of validation errors with the elements of the given list.
     */
    public void setValidationErrors(final Collection<SAMValidationError> errors) {
        mValidationErrors.clear();
        mValidationErrors.addAll(errors);
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final SAMFileHeader that = (SAMFileHeader) o;

        if (!attributesEqual(that)) return false;
        if (mProgramRecords != null ? !mProgramRecords.equals(that.mProgramRecords) : that.mProgramRecords != null)
            return false;
        if (mReadGroups != null ? !mReadGroups.equals(that.mReadGroups) : that.mReadGroups != null) return false;
        if (mSequenceDictionary != null ? !mSequenceDictionary.equals(that.mSequenceDictionary) : that.mSequenceDictionary != null)
            return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = attributesHashCode();
        result = 31 * result + (mSequenceDictionary != null ? mSequenceDictionary.hashCode() : 0);
        result = 31 * result + (mReadGroups != null ? mReadGroups.hashCode() : 0);
        result = 31 * result + (mProgramRecords != null ? mProgramRecords.hashCode() : 0);
        return result;
    }

    public final SAMFileHeader clone() {
        final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
        codec.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        final StringWriter stringWriter = new StringWriter();
        codec.encode(stringWriter, this);
        return codec.decode(new StringLineReader(stringWriter.toString()), "SAMFileHeader.clone");
    }
}
