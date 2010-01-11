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

import java.util.*;

/**
 * Collection of SAMSequenceRecords.
 */
public class SAMSequenceDictionary {
    private List<SAMSequenceRecord> mSequences = new ArrayList<SAMSequenceRecord>();
    private final Map<String, SAMSequenceRecord> mSequenceMap = new HashMap<String, SAMSequenceRecord>();

    public SAMSequenceDictionary() {
    }

    public SAMSequenceDictionary(final List<SAMSequenceRecord> list) {
        this();
        setSequences(list);
    }

    public List<SAMSequenceRecord> getSequences() {
        return Collections.unmodifiableList(mSequences);
    }

    public SAMSequenceRecord getSequence(final String name) {
        return mSequenceMap.get(name);
    }

    /**
     * Replaces the existing list of SAMSequenceRecords with the given list.
     * @param list This value is used directly, rather than being copied.
     */
    public void setSequences(final List<SAMSequenceRecord> list) {
        mSequences = list;
        mSequenceMap.clear();
        int index = 0;
        for (final SAMSequenceRecord record : list) {
            record.setSequenceIndex(index++);
            if (mSequenceMap.put(record.getSequenceName(), record) != null) {
                throw new IllegalArgumentException("Cannot add sequence that already exists in SAMSequenceDictionary: " +
                        record.getSequenceName());
            }
        }
    }

    public void addSequence(final SAMSequenceRecord sequenceRecord) {
        if (mSequenceMap.containsKey(sequenceRecord.getSequenceName())) {
            throw new IllegalArgumentException("Cannot add sequence that already exists in SAMSequenceDictionary: " +
                    sequenceRecord.getSequenceName());
        }
        sequenceRecord.setSequenceIndex(mSequences.size());
        mSequences.add(sequenceRecord);
        mSequenceMap.put(sequenceRecord.getSequenceName(), sequenceRecord);
    }

    /**
     * @return The SAMSequenceRecord with the given index, or null if index is out of range.
     */
    public SAMSequenceRecord getSequence(final int sequenceIndex) {
        if (sequenceIndex < 0 || sequenceIndex >= mSequences.size()) {
            return null;
        }
        return mSequences.get(sequenceIndex);
    }

    /**
     * @return The index for the given sequence name, or -1 if the name is not found.
     */
    public int getSequenceIndex(final String sequenceName) {
        final SAMSequenceRecord record = mSequenceMap.get(sequenceName);
        if (record == null) {
            return -1;
        }
        return record.getSequenceIndex();
    }

    public int size() {
        return mSequences.size();
    }

    public boolean isEmpty() {
        return mSequences.isEmpty();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        SAMSequenceDictionary that = (SAMSequenceDictionary) o;

        if (!mSequences.equals(that.mSequences)) return false;

        return true;
    }

    @Override
    public int hashCode() {
        return mSequences.hashCode();
    }
}