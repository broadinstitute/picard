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
package org.broad.tribble.index.tabix;

import net.sf.samtools.*;
import org.broad.tribble.Feature;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexCreator;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * IndexCreator for Tabix.
 * Features are expected to be 1-based, inclusive.
 */
public class TabixIndexCreator implements IndexCreator {
    private final TabixFormat formatSpec;
    private final List<BinningIndexContent> indexContents = new ArrayList<BinningIndexContent>();
    private final List<String> sequenceNames = new ArrayList<String>();
    // Merely a faster way to ensure that features are added in a specific sequence name order
    private final Set<String> sequenceNamesSeen = new HashSet<String>();
    // A sequence dictionary is not required, but if it is provided all sequences names must be present in it.
    // It is used to determine the length of a sequence in order to optimize index memory allocation.
    private final SAMSequenceDictionary sequenceDictionary;

    private String currentSequenceName = null;
    private BinningIndexBuilder indexBuilder = null;
    // A feature can't be added to the index until the next feature is added because the next feature
    // defines the location of the end of the previous feature in the output file.
    private TabixFeature previousFeature = null;


    /**
     * @param sequenceDictionary is not required, but if present all features added must refer to sequences in the
     *                           dictionary.  It is used to optimize the memory needed to build the index.
     */
    public TabixIndexCreator(final SAMSequenceDictionary sequenceDictionary,
                             final TabixFormat formatSpec) {
        this.sequenceDictionary = sequenceDictionary;
        this.formatSpec = formatSpec.clone();
    }

    public TabixIndexCreator(final TabixFormat formatSpec) {
        this(null, formatSpec);
    }

    @Override
    public void addFeature(final Feature feature, final long filePosition) {
        final String sequenceName = feature.getChr();
        final int referenceIndex;
        if (sequenceName.equals(currentSequenceName)) {
            referenceIndex = sequenceNames.size() - 1;
        } else {
            referenceIndex = sequenceNames.size();
            if (currentSequenceName != null && sequenceNamesSeen.contains(sequenceName)) {
                throw new IllegalArgumentException("Sequence " + feature + " added out sequence of order");
            }
        }
        final TabixFeature thisFeature = new TabixFeature(referenceIndex, feature.getStart(), feature.getEnd(), filePosition);
        if (previousFeature != null) {
            if (previousFeature.compareTo(thisFeature) > 0) {
                throw new IllegalArgumentException(String.format("Features added out of order: previous (%s) > next (%s)",
                        previousFeature, thisFeature));
            }
            finalizeFeature(filePosition);
        }
        previousFeature = thisFeature;
        if (referenceIndex == sequenceNames.size()) {
            advanceToReference(sequenceName);
        }
    }

    private void finalizeFeature(final long featureEndPosition) {
        previousFeature.featureEndFilePosition = featureEndPosition;
        if (previousFeature.featureStartFilePosition >= previousFeature.featureEndFilePosition) {
            throw new IllegalArgumentException(String.format("Feature start position %d >= feature end position %d",
                    previousFeature.featureStartFilePosition, previousFeature.featureEndFilePosition));
        }
        indexBuilder.processFeature(previousFeature);
    }

    private void advanceToReference(final String sequenceName) {
        if (indexBuilder != null) {
            indexContents.add(indexBuilder.generateIndexContent());
        }
        // If sequence dictionary is provided, BinningIndexBuilder can reduce size of array it allocates.
        final int sequenceLength;
        if (sequenceDictionary != null) {
            sequenceLength = sequenceDictionary.getSequence(sequenceName).getSequenceLength();
        } else {
            sequenceLength = 0;
        }
        indexBuilder = new BinningIndexBuilder(sequenceNames.size(), sequenceLength);
        sequenceNames.add(sequenceName);
        currentSequenceName = sequenceName;
        sequenceNamesSeen.add(sequenceName);
    }

    @Override
    public Index finalizeIndex(final long finalFilePosition) {
        if (previousFeature != null) {
            finalizeFeature(finalFilePosition);
        }
        if (indexBuilder != null) {
            indexContents.add(indexBuilder.generateIndexContent());
        }
        // Make this as big as the sequence dictionary, even if there is not content for every sequence,
        // but truncate the sequence dictionary before its end if there are sequences in the sequence dictionary without
        // any features.
        final BinningIndexContent[] indices = indexContents.toArray(new BinningIndexContent[sequenceNames.size()]);
        return new TabixIndex(formatSpec, sequenceNames, indices);
    }


    private static class TabixFeature implements BinningIndexBuilder.FeatureToBeIndexed, Comparable<TabixFeature> {
        private final int referenceIndex;
        private final int start;
        private final int end;
        private final long featureStartFilePosition;
        // Position after this feature in the file.
        private long featureEndFilePosition = -1;

        private TabixFeature(final int referenceIndex, final int start, final int end, final long featureStartFilePosition) {
            this.referenceIndex = referenceIndex;
            this.start = start;
            this.end = end;
            this.featureStartFilePosition = featureStartFilePosition;
        }

        @Override
        public int getStart() {
            return start;
        }

        @Override
        public int getEnd() {
            return end;
        }

        /**
         *
         * @return null -- Let index builder compute this.
         */
        @Override
        public Integer getIndexingBin() {
            return null;
        }

        @Override
        public Chunk getChunk() {
            if (featureEndFilePosition == -1) {
                throw new IllegalStateException("End position is not set");
            }
            return new Chunk(featureStartFilePosition, featureEndFilePosition);
        }

        @Override
        public int compareTo(final TabixFeature other) {
            final int ret = this.referenceIndex - other.referenceIndex;
            if (ret != 0) return ret;
            return this.start - other.start;
        }

        @Override
        public String toString() {
            return "TabixFeature{" +
                    "referenceIndex=" + referenceIndex +
                    ", start=" + start +
                    ", end=" + end +
                    ", featureStartFilePosition=" + featureStartFilePosition +
                    ", featureEndFilePosition=" + featureEndFilePosition +
                    '}';
        }
    }
}
