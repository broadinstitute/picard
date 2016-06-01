/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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
package picard.illumina.parser.readers;

import picard.PicardException;

import java.io.File;
import java.nio.ByteBuffer;

/**
 * Annoyingly, there are two different files with extension .bci in NextSeq output.  This reader handles
 * the file that contains virtual file pointers into a .bcl.bgzf file.  After the header, there is a 64-bit record
 * per tile.
 */
public class BclIndexReader {
    private static final int BCI_HEADER_SIZE = 8;
    private static final int BCI_VERSION = 0;

    private final BinaryFileIterator<Long> bciIterator;
    private final int numTiles;
    private final File bciFile;
    private int nextRecordNumber = 0;

    public BclIndexReader(final File bclFile) {
        bciFile = new File(bclFile.getAbsolutePath() + ".bci");
        bciIterator = MMapBackedIteratorFactory.getLongIterator(BCI_HEADER_SIZE, bciFile);
        final ByteBuffer headerBytes = bciIterator.getHeaderBytes();
        final int actualVersion = headerBytes.getInt();
        if (actualVersion != BCI_VERSION) {
            throw new PicardException(String.format("Unexpected version number %d in %s", actualVersion, bciFile.getAbsolutePath()));
        }
        numTiles = headerBytes.getInt();
    }

    public int getNumTiles() {
        return numTiles;
    }

    public long get(final int recordNumber) {
        if (recordNumber < nextRecordNumber) {
            throw new IllegalArgumentException("Can only read forward");
        }
        if (recordNumber > nextRecordNumber) {
            bciIterator.skipElements(recordNumber - nextRecordNumber);
            nextRecordNumber = recordNumber;
        }
        ++nextRecordNumber;
        return bciIterator.getElement();
    }

    public File getBciFile() {
        return bciFile;
    }
}
