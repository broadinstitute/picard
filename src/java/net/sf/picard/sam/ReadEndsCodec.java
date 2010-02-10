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
package net.sf.picard.sam;

import net.sf.samtools.util.SortingCollection;
import net.sf.picard.PicardException;

import java.io.*;

/** Coded for ReadEnds that just outputs the primitive fields and reads them back. */
class ReadEndsCodec implements SortingCollection.Codec<ReadEnds> {
    private DataInputStream in;
    private DataOutputStream out;

    public SortingCollection.Codec<ReadEnds> clone() {
        return new ReadEndsCodec();
    }

    public void setOutputStream(final OutputStream os) { this.out = new DataOutputStream(os); }
    public void setInputStream(final InputStream is) { this.in = new DataInputStream(is); }

    public DataInputStream getInputStream() {
        return in;
    }

    public DataOutputStream getOutputStream() {
        return out;
    }

    public void encode(final ReadEnds read) {
        try {
            this.out.writeShort(read.score);
            this.out.writeShort(read.libraryId);
            this.out.writeByte(read.orientation);
            this.out.writeInt(read.read1Sequence);
            this.out.writeInt(read.read1Coordinate);
            this.out.writeLong(read.read1IndexInFile);
            this.out.writeInt(read.read2Sequence);

            if (read.orientation > ReadEnds.R) {
                this.out.writeInt(read.read2Coordinate);
                this.out.writeLong(read.read2IndexInFile);
            }

            this.out.writeShort(read.readGroup);
            this.out.writeByte(read.tile);
            this.out.writeShort(read.x);
            this.out.writeShort(read.y);

            this.out.flush();
        }
        catch (IOException ioe) {
            throw new PicardException("Exception writing ReadEnds to file.", ioe);
        }
    }

    public ReadEnds decode() {
        final ReadEnds read = new ReadEnds();
        try {
            // If the first read results in an EOF we've exhausted the stream
            try { read.score = this.in.readShort(); }
            catch (EOFException eof) { return null; }

            read.libraryId        = this.in.readShort();
            read.orientation      = this.in.readByte();
            read.read1Sequence    = this.in.readInt();
            read.read1Coordinate  = this.in.readInt();
            read.read1IndexInFile = this.in.readLong();
            read.read2Sequence    = this.in.readInt();

            if (read.orientation > ReadEnds.R) {
                read.read2Coordinate  = this.in.readInt();
                read.read2IndexInFile = this.in.readLong();
            }

            read.readGroup = this.in.readShort();
            read.tile      = this.in.readByte();
            read.x         = this.in.readShort();
            read.y         = this.in.readShort();

            return read;
        }
        catch (IOException ioe) {
            throw new PicardException("Exception writing ReadEnds to file.", ioe);
        }
    }
}
