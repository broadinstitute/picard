package org.broad.tribble.readers;

import net.sf.samtools.util.AbstractIterator;
import net.sf.samtools.util.CloserUtil;

import java.io.Closeable;
import java.io.IOException;

/** A simple iterator over the elements in LineReader. */
public class LineIteratorImpl extends AbstractIterator<String> implements LineIterator, Closeable {
    private final LineReader lineReader;

    /**
     * @param lineReader The line reader whose elements are to be iterated over.
     */
    public LineIteratorImpl(final LineReader lineReader) {
        this.lineReader = lineReader;
    }

    @Override
    protected String advance() {
        try {
            return lineReader.readLine();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void close() throws IOException {
        CloserUtil.close(lineReader);
    }
}
