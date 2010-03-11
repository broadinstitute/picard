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
 * In-memory representation of @PG SAM header record.
 */
public class SAMProgramRecord extends AbstractSAMHeaderRecord {
    public static final String PROGRAM_GROUP_ID_TAG = "ID";
    public static final String PROGRAM_NAME_TAG = "PN";
    public static final String PROGRAM_VERSION_TAG = "VN";
    public static final String COMMAND_LINE_TAG = "CL";
    public static final String PREVIOUS_PROGRAM_GROUP_ID_TAG = "PP";
    private String mProgramGroupId;
    public static final Set<String> STANDARD_TAGS = Collections.unmodifiableSet(
            new HashSet<String>(Arrays.asList(PROGRAM_GROUP_ID_TAG,
                    PROGRAM_NAME_TAG,
                    PROGRAM_VERSION_TAG,
                    COMMAND_LINE_TAG,
                    PREVIOUS_PROGRAM_GROUP_ID_TAG)) );

    public SAMProgramRecord(final String programGroupId) {
        this.mProgramGroupId = programGroupId;
    }

    public SAMProgramRecord(final String id, SAMProgramRecord srcProgramRecord) {
        mProgramGroupId = id;
        for (final Map.Entry<String, Object> entry : srcProgramRecord.getAttributes()) {
            setAttribute(entry.getKey(), entry.getValue());
        }
    }

    public String getId() {
        return getProgramGroupId();
    }

    public String getProgramGroupId() {
        return mProgramGroupId;
    }

    public String getProgramName() {
        return (String)getAttribute(PROGRAM_NAME_TAG);
    }

    public void setProgramName(final String name) {
        setAttribute(PROGRAM_NAME_TAG, name);
    }

    public String getProgramVersion() {
        return (String)getAttribute(PROGRAM_VERSION_TAG);
    }

    public void setProgramVersion(final String version) {
        setAttribute(PROGRAM_VERSION_TAG, version);
    }

    public String getCommandLine() {
        return (String)getAttribute(COMMAND_LINE_TAG);
    }

    public void setCommandLine(final String commandLine) {
        setAttribute(COMMAND_LINE_TAG, commandLine);
    }

    public String getPreviousProgramGroupId() {
        return (String)getAttribute(PREVIOUS_PROGRAM_GROUP_ID_TAG);
    }

    public void setPreviousProgramGroupId(final String id) {
        setAttribute(PREVIOUS_PROGRAM_GROUP_ID_TAG, id);
    }



    /**
     * @return true if this == that except for the program group ID, which is arbitrary
     */
    public boolean equivalent(final SAMProgramRecord that) {
        return attributesEqual(that);
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final SAMProgramRecord that = (SAMProgramRecord) o;

        if (!attributesEqual(that)) return false;
        if (mProgramGroupId != null ? !mProgramGroupId.equals(that.mProgramGroupId) : that.mProgramGroupId != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = mProgramGroupId != null ? mProgramGroupId.hashCode() : 0;
        result = 31 * result + attributesHashCode();
        return result;
    }

    Set<String> getStandardTags() {
        return STANDARD_TAGS;
    }
}
