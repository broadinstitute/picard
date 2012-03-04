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


import net.sf.samtools.util.Iso8601Date;

import java.util.*;

/**
 * Header information about a read group.
 */
public class SAMReadGroupRecord extends AbstractSAMHeaderRecord
{
    private String mReadGroupId = null;
    public static final String READ_GROUP_ID_TAG = "ID";
    public static final String READ_GROUP_SAMPLE_TAG = "SM";
    public static final String PREDICTED_MEDIAN_INSERT_SIZE_TAG = "PI";
    public static final String DATE_RUN_PRODUCED_TAG = "DT";
    public static final String FLOW_ORDER_TAG = "FO";
    public static final String KEY_SEQUENCE_TAG = "KS";
    public static final String DESCRIPTION_TAG = "DS";
    public static final String PLATFORM_UNIT_TAG = "PU";
    public static final String SEQUENCING_CENTER_TAG = "CN";
    public static final String PLATFORM_TAG = "PL";
    public static final String LIBRARY_TAG = "LB";

    public static final Set<String> STANDARD_TAGS =
            new HashSet<String>(Arrays.asList(READ_GROUP_ID_TAG, READ_GROUP_SAMPLE_TAG, LIBRARY_TAG,
        DESCRIPTION_TAG, PLATFORM_UNIT_TAG, PREDICTED_MEDIAN_INSERT_SIZE_TAG, SEQUENCING_CENTER_TAG,
            DATE_RUN_PRODUCED_TAG, PLATFORM_TAG, FLOW_ORDER_TAG, KEY_SEQUENCE_TAG));

    public SAMReadGroupRecord(final String id) { mReadGroupId = id; }

    public SAMReadGroupRecord(final String id, final SAMReadGroupRecord srcProgramRecord) {
        mReadGroupId = id;
        for (final Map.Entry<String, String> entry : srcProgramRecord.getAttributes()) {
            setAttribute(entry.getKey(), entry.getValue());
        }
    }

    public String getId() { return getReadGroupId();  }
    public String getReadGroupId() { return mReadGroupId; }

    public String getSample() { return getAttribute(READ_GROUP_SAMPLE_TAG); }
    public void setSample(final String value) { setAttribute(READ_GROUP_SAMPLE_TAG, value); }

    public String getLibrary() { return getAttribute(LIBRARY_TAG); }
    public void setLibrary(final String value) { setAttribute(LIBRARY_TAG, value); }

    public String getPlatformUnit() { return getAttribute(PLATFORM_UNIT_TAG); }
    public void setPlatformUnit(final String pu) { setAttribute(PLATFORM_UNIT_TAG, pu); }

    public String getPlatform() { return getAttribute(PLATFORM_TAG); }
    public void setPlatform(final String platform) { setAttribute(PLATFORM_TAG, platform); }

    public Date getRunDate() {
        final String dt = getAttribute(DATE_RUN_PRODUCED_TAG);
        if (dt == null) return null;
        else return new Iso8601Date(dt);
    }

    public String getFlowOrder() { return getAttribute(FLOW_ORDER_TAG); }
    public void setFlowOrder(final String flowOrder) { setAttribute(FLOW_ORDER_TAG, flowOrder); }

    public String getKeySequence() { return getAttribute(KEY_SEQUENCE_TAG); }
    public void setKeySequence(final String keySequence) { setAttribute(KEY_SEQUENCE_TAG, keySequence); }

    /**
     * Converts to Iso8601Date if not already in that form.
     */
    public void setRunDate(Date runDate) {
        if (runDate != null && !(runDate instanceof Iso8601Date)) {
            runDate = new Iso8601Date(runDate);
        }
        setAttribute(DATE_RUN_PRODUCED_TAG, runDate != null ? runDate.toString() : null);
    }

    public String getSequencingCenter() { return getAttribute(SEQUENCING_CENTER_TAG); }
    public void setSequencingCenter(final String center) { setAttribute(SEQUENCING_CENTER_TAG, center); }

    public String getDescription() { return getAttribute(DESCRIPTION_TAG); }
    public void setDescription(final String description) { setAttribute(DESCRIPTION_TAG, description); }

    public Integer getPredictedMedianInsertSize() {
        final String stringRep = getAttribute(PREDICTED_MEDIAN_INSERT_SIZE_TAG);
        if (stringRep == null) return null;
        return Integer.parseInt(stringRep); 
    }
    public void setPredictedMedianInsertSize(final Integer predictedMedianInsertSize) {

        setAttribute(PREDICTED_MEDIAN_INSERT_SIZE_TAG, (predictedMedianInsertSize == null? null: predictedMedianInsertSize.toString())); 
    }

    /**
     * @return true if this == that except for the read group ID, which is arbitrary
     */
    public boolean equivalent(final SAMReadGroupRecord that) {
        return attributesEqual(that);
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final SAMReadGroupRecord that = (SAMReadGroupRecord) o;

        if (!attributesEqual(that)) return false;
        if (mReadGroupId != null ? !mReadGroupId.equals(that.mReadGroupId) : that.mReadGroupId != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        return mReadGroupId.hashCode();
    }

    Set<String> getStandardTags() {
        return STANDARD_TAGS;
    }
}

