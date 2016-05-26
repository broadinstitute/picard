/*
  * The MIT License
  *
  * Copyright (c) 2015 The Broad Institute
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

import java.io.EOFException;
import java.io.IOException;

/**
 * Created by nhomer on 9/13/15.
 */
public class ReadEndsForMarkDuplicatesSetSizeTagsCodec extends ReadEndsForMarkDuplicatesCodec {

    @Override
    public SortingCollection.Codec<ReadEndsForMarkDuplicates> clone() {
        return new ReadEndsForMarkDuplicatesSetSizeTagsCodec();
    }

    @Override
    public void encode(final ReadEndsForMarkDuplicates read) {
        if (!(read instanceof ReadEndsForMarkDuplicatesSetSizeTags)) {
            throw new PicardException("Read was not a ReadEndsForMarkDuplicatesSetSizeTagscodes");
        }
        super.encode(read);

        try {
            final ReadEndsForMarkDuplicatesSetSizeTags val = (ReadEndsForMarkDuplicatesSetSizeTags)read;
            out.writeBytes(val.firstEncounteredReadName);
            out.writeBytes(val.representativeReadName);
            out.writeInt(val.duplicateSetSize);
        } catch (final IOException ioe) {
            throw new PicardException("Exception writing ReadEnds to file.", ioe);
        }
    }

    @Override
    public ReadEndsForMarkDuplicates decode() {
        final ReadEndsForMarkDuplicates parentRead = super.decode();
        if (null == parentRead) return null; // EOF
        final ReadEndsForMarkDuplicatesSetSizeTags read = new ReadEndsForMarkDuplicatesSetSizeTags(parentRead);
        try {
            read.firstEncounteredReadName = in.readUTF();
            read.representativeReadName = in.readUTF();
            read.duplicateSetSize = in.readInt();
            return read;
        } catch (final IOException ioe) {
            throw new PicardException("Exception writing ReadEnds to file.", ioe);
        }
    }

}


