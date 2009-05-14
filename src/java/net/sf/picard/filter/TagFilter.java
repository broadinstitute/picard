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
package net.sf.picard.filter;

import java.util.List;
import java.util.Arrays;

import net.sf.samtools.SAMRecord;

/**
 * Filter class for matching tag attributes in SAMRecords
 */
public class TagFilter implements SamRecordFilter {

    private final String tag;           // The key of the tag to match
    private final List<Object> values;  // The list of matching values

    /**
     * Constructor for a single value
     *
     * @param tag       the key of the tag to match
     * @param value     the value to match
     */
    public TagFilter(String tag, Object value) {
        this.tag = tag;
        this.values = Arrays.asList(value);
    }

    /**
     * Constructor for multiple values
     *
     * @param tag       the key of the tag to match
     * @param values    the matching values
     */
    public TagFilter(String tag, List<Object> values) {
        this.tag = tag;
        this.values = values;
    }

    /**
     * Determines whether a SAMRecord matches this filter
     *
     * @param record    the SAMRecord to evaluate
     * @return  true if the SAMRecord matches the filter, otherwise false
     */
    public boolean filterOut(SAMRecord record) {
        return values.contains(record.getAttribute(tag));
    }
 }
