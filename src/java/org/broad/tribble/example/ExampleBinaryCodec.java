/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
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
package org.broad.tribble.example;

import org.broad.tribble.*;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.readers.PositionalBufferedStream;

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * An example binary codec that encodes / decodes contig / start / stop values via DataInputStreams
 *
 * @author Mark DePristo
 */
public class ExampleBinaryCodec extends BinaryFeatureCodec<Feature> {
    public final static String HEADER_LINE = "# BinaryTestFeature";

    @Override
    public Feature decodeLoc(final PositionalBufferedStream stream) throws IOException {
        return decode(stream);
    }

    @Override
    public Feature decode(final PositionalBufferedStream stream) throws IOException {
        DataInputStream dis = new DataInputStream(stream);
        String contig = dis.readUTF();
        int start = dis.readInt();
        int stop = dis.readInt();
        return new BasicFeature(contig, start, stop);
    }

    @Override
    public FeatureCodecHeader readHeader(final PositionalBufferedStream stream) throws IOException {
        // Construct a reader that does not read ahead (because we don't want to consume data from the stream that is not the header)
        final AsciiLineReader nonReadAheadLineReader = new AsciiLineReader(stream);
        final List<String> headerLines = new ArrayList<String>();
        long headerLengthInBytes = 0;
        while (stream.peek() == ('#' & 0xff)) { // Look for header lines, which are prefixed by '#'.
            headerLines.add(nonReadAheadLineReader.readLine());
            headerLengthInBytes = stream.getPosition();
        }
        return new FeatureCodecHeader(headerLines, headerLengthInBytes);
    }

    @Override
    public Class<Feature> getFeatureType() {

        return Feature.class;
    }
    @Override
    public boolean canDecode(final String path) {
        return false;
    }

    /**
     * Convenience method that creates an ExampleBinaryCodec file from another feature file.
     *
     * For testing purposes really
     *
     * @param source file containing the features
     * @param dest the place to write the binary features
     * @param codec of the source file features
     * @throws IOException
     */
    public static <FEATURE_TYPE extends Feature> void convertToBinaryTest(final File source, final File dest, final FeatureCodec<FEATURE_TYPE, LineIterator> codec) throws IOException {
        final FeatureReader<FEATURE_TYPE> reader = AbstractFeatureReader.getFeatureReader(source.getAbsolutePath(), codec, false); // IndexFactory.loadIndex(idxFile));
        final OutputStream output = new FileOutputStream(dest);
        ExampleBinaryCodec.convertToBinaryTest(reader, output);
    }

    /**
     * Convenience method that creates an ExampleBinaryCodec file from another feature file.
     *
     * For testing purposes really
     *
     * @throws IOException
     */
    public static <FEATURE_TYPE extends Feature> void convertToBinaryTest(final FeatureReader<FEATURE_TYPE> reader, final OutputStream out) throws IOException {
        DataOutputStream dos = new DataOutputStream(out);
        dos.writeBytes(HEADER_LINE + "\n");
        Iterator<FEATURE_TYPE> it = reader.iterator();
        while ( it.hasNext() ) {
            final Feature f = it.next();
            dos.writeUTF(f.getChr());
            dos.writeInt(f.getStart());
            dos.writeInt(f.getEnd());
        }
        dos.close();
        reader.close();
    }
}
