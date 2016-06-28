package picard.util;

import htsjdk.samtools.util.CloseableIterator;
import picard.PicardException;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

/**
 * Iterate through a delimited text file in which columns are found by looking at a header line rather than by position.
 *
 * TODO: This effectively replaces TabbedTextFileWithHeaderParser although the latter hasn't been modified to use this
 * code instead.
 *
 * @author jgentry@broadinstitute.org
 */
public class DelimitedTextFileWithHeaderIterator implements CloseableIterator<DelimitedTextFileWithHeaderIterator.Row> {
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
            final Integer key = columnLabelIndices.get(columnLabel);
            if (key == null) throw new NoSuchElementException(String.format("column %s in %s", columnLabel, parser.getFileName()));
            return fields[key];
        }

        public Integer getIntegerField(final String columnLabel) {
            if (fields[columnLabelIndices.get(columnLabel)] == null)  return null;
            return Integer.parseInt(fields[columnLabelIndices.get(columnLabel)]);
        }

        public String getCurrentLine() {
            return this.currentLine;
        }
    }

    /**
     * Map from column label to positional index.
     */
    private final Map<String, Integer> columnLabelIndices = new HashMap<String, Integer>();
    private final BasicInputParser parser;

    public DelimitedTextFileWithHeaderIterator(final BasicInputParser parser) {
        this.parser = parser;
        if (!parser.hasNext()) {
            throw new PicardException("No header line found in file " + parser.getFileName());
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
     *
     * @return The set of column labels for this file in no particular order.
     */
    public Set<String> columnLabels() {
        return columnLabelIndices.keySet();
    }

    public int getCurrentLineNumber() {
        return parser.getCurrentLineNumber();
    }

    public Set<String> getColumnNames() {
        return Collections.unmodifiableSet(this.columnLabelIndices.keySet());
    }

    @Override
    public boolean hasNext() {
        return parser.hasNext();
    }

    @Override
    public Row next() {
        final String[] fields = parser.next();
        final String source = parser.getCurrentLine();
        return new Row(fields, source);
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void close() {
        parser.close();
    }
}
