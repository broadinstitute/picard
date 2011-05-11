/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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
package net.sf.picard.util;

import net.sf.picard.PicardException;
import net.sf.samtools.util.CloseableIterator;

import java.io.File;
import java.util.ConcurrentModificationException;
import java.util.HashMap;
import java.util.Map;

/**
 * Parse a tabbed text file in which columns are found by looking at a header line rather than by position.
 *
 * @author alecw@broadinstitute.org
 */
public class TabbedTextFileWithHeaderParser implements Iterable<TabbedTextFileWithHeaderParser.Row> {
    public class Row {
        private final String[] fields;
        private final String currentLine;

        Row(final String[] fields, final String source) {
            this.fields = fields;
            this.currentLine = source;
        }

        /**
         * @return Array of fields in the order they appear in the file.
         */
        public String[] getFields() {
            return fields;
        }

        public String getField(final String columnLabel) {
            return fields[columnLabelIndices.get(columnLabel)];
        }

        public String getCurrentLine() {
            return this.currentLine;
        }



    }

    class TheIterator implements CloseableIterator<Row> {

        @Override
        public boolean hasNext() {
            return parser.hasNext();
        }

        @Override
        public Row next() {
            final String source = parser.getCurrentLine();
            final String[] fields = parser.next();
            return new Row(fields, source);
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }

        @Override
        public void close() {
            extantIterator = null;
        }
    }

    /**
     * Map from column label to positional index.
     */
    private final Map<String, Integer> columnLabelIndices = new HashMap<String, Integer>();
    private final TabbedInputParser parser;
    private TheIterator extantIterator;

    public TabbedTextFileWithHeaderParser(final File file) {
        parser = new TabbedInputParser(false, file);
        if (!parser.hasNext()) {
            throw new PicardException("No header line found in file " + file);
        }
        final String[] columnLabels = parser.next();
        for (int i = 0; i < columnLabels.length; ++i) {
            columnLabelIndices.put(columnLabels[i], i);
        }
    }

    /**
     * @param columnLabel
     * @return True if the given column label appears in the header.
     */
    public boolean hasColumn(final String columnLabel) {
        return columnLabelIndices.containsKey(columnLabel);
    }

    /**
     * Creates the iterator object.  It is illegal to have more than one iterator extant
     * on the same parser object.
     */
    @Override
    public CloseableIterator<Row> iterator() {
        if (extantIterator != null) {
            throw new ConcurrentModificationException("Only one iterator allowed at a time.");
        }
        extantIterator = new TheIterator();
        return extantIterator;
    }

    /**
     * Release all resources associated with the parser.  Iteration will not work after this
     * has been called.
     */
    public void close() {
        parser.close();
    }
}
