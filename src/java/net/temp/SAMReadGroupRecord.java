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
    public static final String DESCRIPTION_TAG = "DS";
    public static final String PLATFORM_UNIT_TAG = "PU";
    public static final String SEQUENCING_CENTER_TAG = "CN";
    public static final String PLATFORM_TAG = "PL";
    public static final String LIBRARY_TAG = "LB";

    public static final Set<String> STANDARD_TAGS =
            new HashSet<String>(Arrays.asList(READ_GROUP_ID_TAG, READ_GROUP_SAMPLE_TAG, LIBRARY_TAG,
        DESCRIPTION_TAG, PLATFORM_UNIT_TAG, PREDICTED_MEDIAN_INSERT_SIZE_TAG, SEQUENCING_CENTER_TAG,
            DATE_RUN_PRODUCED_TAG, PLATFORM_TAG));

    public SAMReadGroupRecord(final String id) { mReadGroupId = id; }

    public SAMReadGroupRecord(final String id, SAMReadGroupRecord srcProgramRecord) {
        mReadGroupId = id;
        for (final Map.Entry<String, Object> entry : srcProgramRecord.getAttributes()) {
            setAttribute(entry.getKey(), entry.getValue());
        }
    }

    public String getReadGroupId() { return mReadGroupId; }

    public String getSample() { return (String) getAttribute("SM"); }
    public void setSample(final String value) { setAttribute("SM", value); }

    public String getLibrary() { return (String) getAttribute("LB"); }
    public void setLibrary(final String value) { setAttribute("LB", value); }

    public String getPlatformUnit() { return (String) getAttribute(PLATFORM_UNIT_TAG); }
    public void setPlatformUnit(final String pu) { setAttribute(PLATFORM_UNIT_TAG, pu); }

    public String getPlatform() { return (String) getAttribute(PLATFORM_TAG); }
    public void setPlatform(final String platform) { setAttribute(PLATFORM_TAG, platform); }

    public Date getRunDate() { return (Date) getAttribute(DATE_RUN_PRODUCED_TAG); }

    /**
     * Converts to Iso8601Date if not already in that form.
     */
    public void setRunDate(Date runDate) {
        if (runDate != null && !(runDate instanceof Iso8601Date)) {
            runDate = new Iso8601Date(runDate);
        }
        setAttribute(DATE_RUN_PRODUCED_TAG, runDate);
    }

    public String getSequencingCenter() { return (String) getAttribute(SEQUENCING_CENTER_TAG); }
    public void setSequencingCenter(final String center) { setAttribute(SEQUENCING_CENTER_TAG, center); }

    public String getDescription() { return (String) getAttribute(DESCRIPTION_TAG); }
    public void setDescription(final String description) { setAttribute(DESCRIPTION_TAG, description); }

    public Integer getPredictedMedianInsertSize() { return (Integer) getAttribute(PREDICTED_MEDIAN_INSERT_SIZE_TAG); }
    public void setPredictedMedianInsertSize(final Integer predictedMedianInsertSize) { setAttribute(PREDICTED_MEDIAN_INSERT_SIZE_TAG, predictedMedianInsertSize); }

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
        int result = mReadGroupId != null ? mReadGroupId.hashCode() : 0;
        result = 31 * result + attributesHashCode();
        return result;
    }

    Set<String> getStandardTags() {
        return STANDARD_TAGS;
    }
}

