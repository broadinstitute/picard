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
package net.sf.samtools.util;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * Use this type rather than java.util.Date in command-line options in order to get ISO 8601 parsing.
 * The ctors below truncate milliseconds, since our formatter does not write them.  Note that it is possible
 * to modify an Iso8601Date so that it has fractional seconds, but that is discouraged.
 *
 * @author alecw@broadinstitute.org
 */
public class Iso8601Date extends Date {
    private static final ThreadLocal<DateFormat> iso8601DateFormatter = new ThreadLocal<DateFormat>() {
        protected synchronized DateFormat initialValue() {
            return new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ssZ");
        }
    };

    public Iso8601Date(final String dateStr) {
        super(DateParser.parse(dateStr).getTime());
        truncateMilliseconds();
    }

    public Iso8601Date(final Date date) {
        super(date.getTime());
        truncateMilliseconds();
    }

    public String toString() {
        return iso8601DateFormatter.get().format(this);
    }

    private void truncateMilliseconds() {
        long time = getTime();
        long mod = time % 1000;
        if (mod != 0) {
            setTime(time - mod);
        }
    }
}
