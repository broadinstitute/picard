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
package net.sf.picard.annotation;

import net.sf.picard.util.Interval;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * Holds annotation of a gene for storage in an OverlapDetector.  May hold multiple transcripts for the same gene.
 * The transcripts must all be relative to the same strand.
 */
public class Gene extends Interval implements Iterable<Transcript>  {
    private Map<String, Transcript> transcripts = new HashMap<String, Transcript>();


    public Gene(String sequence, int start, int end, boolean negative, String name) {
        super(sequence, start, end, negative, name);
    }

    public Gene(String sequence, int start, int end, boolean negative, String name, Iterable<Transcript> transcriptIterable) {
        super(sequence, start, end, negative, name);
        for (final Transcript transcript : transcriptIterable) {
            addTranscript(transcript);
        }
    }

    public void addTranscript(final Transcript transcript) {
        if (transcripts.containsKey(transcript.name)) {
            throw new AnnotationException("Transcript " + transcript.name + " for gene " +
            this.getName() + " appears more than once -- skipping gene");
        }
        transcripts.put(transcript.name, transcript);
    }

    public Iterator<Transcript> iterator() {
        return transcripts.values().iterator();
    }

}
