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
package picard.sam.markduplicates.util;

import htsjdk.samtools.util.SortingCollection;
import picard.PicardException;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/** Codec for ReadEnds that just outputs the primitive fields and reads them back. */
public class ReadEndsForMarkDuplicatesCodec implements SortingCollection.Codec<ReadEndsForMarkDuplicates> {
    protected DataInputStream in;
    protected DataOutputStream out;

    final protected double scaleFactor;

    public ReadEndsForMarkDuplicatesCodec(final double coordinateAccuracy) {
        this.scaleFactor = coordinateAccuracy;
    }

    public SortingCollection.Codec<ReadEndsForMarkDuplicates> clone() {
        return new ReadEndsForMarkDuplicatesCodec(this.scaleFactor);
    }

    public void setOutputStream(final OutputStream os) { this.out = new DataOutputStream(os); }

    public void setInputStream(final InputStream is) { this.in = new DataInputStream(is); }

    public DataInputStream getInputStream() {
        return in;
    }

    public DataOutputStream getOutputStream() {
        return out;
    }

    public void encode(final ReadEndsForMarkDuplicates read) {
        try {
            this.out.writeShort(read.score);
            this.out.writeShort(read.libraryId);
            this.out.writeByte(read.orientation);
            this.out.writeInt(read.read1ReferenceIndex);
            this.out.writeInt(read.read1Coordinate);
            this.out.writeLong(read.read1IndexInFile);
            this.out.writeInt(read.read2ReferenceIndex);

            if (read.orientation > ReadEnds.R) {
                this.out.writeInt(read.read2Coordinate);
                this.out.writeLong(read.read2IndexInFile);
            }

            this.out.writeShort(read.readGroup);
            this.out.writeShort(read.tile);

            // scaling this so that in cases where there may be overflow, but the local accuracy is not so important
            // the low-level bits are forgotten instead of overflowing

            this.out.writeShort((short) (read.x / this.scaleFactor));
            this.out.writeShort((short) (read.y / this.scaleFactor));
            this.out.writeByte(read.orientationForOpticalDuplicates);
            this.out.writeInt(read.duplicateSetSize);
        } catch (final IOException ioe) {
            throw new PicardException("Exception writing ReadEnds to file.", ioe);
        }
    }

    public ReadEndsForMarkDuplicates decode() {
        final ReadEndsForMarkDuplicates read = new ReadEndsForMarkDuplicates();
        try {
            // If the first read results in an EOF we've exhausted the stream
            try {
                read.score = this.in.readShort();
            } catch (final EOFException eof) {
                return null;
            }

            read.libraryId = this.in.readShort();
            read.orientation = this.in.readByte();
            read.read1ReferenceIndex = this.in.readInt();
            read.read1Coordinate = this.in.readInt();
            read.read1IndexInFile = this.in.readLong();
            read.read2ReferenceIndex = this.in.readInt();

            if (read.orientation > ReadEnds.R) {
                read.read2Coordinate = this.in.readInt();
                read.read2IndexInFile = this.in.readLong();
            }

            read.readGroup = this.in.readShort();
            read.tile = this.in.readShort();

            // reversing the scaling that was done during writing.
            read.x = (int) (this.scaleFactor * this.in.readShort());
            read.y = (int) (this.scaleFactor * this.in.readShort());

            read.orientationForOpticalDuplicates = this.in.readByte();
            read.duplicateSetSize = this.in.readInt();

            return read;
        } catch (final IOException ioe) {
            throw new PicardException("Exception writing ReadEnds to file.", ioe);
        }
    }
}
