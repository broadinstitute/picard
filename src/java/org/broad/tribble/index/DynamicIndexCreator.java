/*
 * Copyright (c) 2010, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broad.tribble.index;

import org.broad.tribble.Feature;
import org.broad.tribble.TribbleException;
import org.broad.tribble.index.interval.IntervalIndexCreator;
import org.broad.tribble.index.linear.LinearIndexCreator;
import org.broad.tribble.util.MathUtils;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeMap;


/**
 * A DynamicIndexCreator creates the proper index based on an {@link IndexFactory.IndexBalanceApproach} and
 * the characteristics of the file.  Ultimately this is either a LinearIndex or an IntervalTreeIndex, with index
 * parameters based on whether seek time or file size is to be minimized.
 */
public class DynamicIndexCreator extends TribbleIndexCreator {
    IndexFactory.IndexBalanceApproach iba;
    Map<IndexFactory.IndexType,TribbleIndexCreator> creators;

    /**
     * we're interested in two stats:
     * the longest feature and the density of features
      */
    int longestFeatureLength = 0;
    long featureCount = 0;

    MathUtils.RunningStat stats = new MathUtils.RunningStat();
    long basesSeen = 0;
    Feature lastFeature = null;
    File inputFile;

    public DynamicIndexCreator(final File inputFile, final IndexFactory.IndexBalanceApproach iba) {
        this.iba = iba;
        // get a list of index creators
        this.inputFile = inputFile;
        creators = getIndexCreators(inputFile,iba);
    }

    public Index finalizeIndex(final long finalFilePosition) {
        // finalize all of the indexes
        // return the score of the indexes we've generated
        final Map<Double,TribbleIndexCreator> mapping = scoreIndexes((double)featureCount/(double)basesSeen, creators, longestFeatureLength, iba);
        final TribbleIndexCreator creator = getMinIndex(mapping, this.iba);

        for (final Map.Entry<String, String> entry : properties.entrySet()) {
            creator.addProperty(entry.getKey(), entry.getValue());
        }

        // add our statistics to the file
        creator.addProperty("FEATURE_LENGTH_MEAN",String.valueOf(stats.mean()));
        creator.addProperty("FEATURE_LENGTH_STD_DEV",String.valueOf(stats.standardDeviation()));
        creator.addProperty("MEAN_FEATURE_VARIANCE",String.valueOf(stats.variance()));

        // add the feature count
        creator.addProperty("FEATURE_COUNT",String.valueOf(featureCount));

        // Now let's finalize and create the index itself
        return creator.finalizeIndex(finalFilePosition);
    }

    /**
     * create a list of index creators (initialized) representing the common index types we'd suspect they'd like to use
     * @param inputFile the input file to use to create the indexes
     * @return a map of index type to the best index for that balancing approach
     */
    private Map<IndexFactory.IndexType,TribbleIndexCreator> getIndexCreators(final File inputFile, final IndexFactory.IndexBalanceApproach iba) {
        final Map<IndexFactory.IndexType,TribbleIndexCreator> creators = new HashMap<IndexFactory.IndexType,TribbleIndexCreator>();

        if (iba == IndexFactory.IndexBalanceApproach.FOR_SIZE) {
            // add a linear index with the default bin size
            final LinearIndexCreator linearNormal = new LinearIndexCreator(inputFile, LinearIndexCreator.DEFAULT_BIN_WIDTH);
            creators.put(IndexFactory.IndexType.LINEAR,linearNormal);

            // create a tree index with the default size
            final IntervalIndexCreator treeNormal = new IntervalIndexCreator(inputFile, IntervalIndexCreator.DEFAULT_FEATURE_COUNT);
            creators.put(IndexFactory.IndexType.INTERVAL_TREE,treeNormal);
        }

        // this section is a little more arbitrary; we're creating indexes with a bin size that's a portion of the default; these
        // values were determined experimentally
        if (iba == IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME) {
            // create a linear index with a small bin size
            final LinearIndexCreator linearSmallBin =
                    new LinearIndexCreator(inputFile, Math.max(200, LinearIndexCreator.DEFAULT_BIN_WIDTH / 4));
            creators.put(IndexFactory.IndexType.LINEAR,linearSmallBin);

            // create a tree index with a small index size
            final IntervalIndexCreator treeSmallBin =
                    new IntervalIndexCreator(inputFile, Math.max(20, IntervalIndexCreator.DEFAULT_FEATURE_COUNT / 8));
            creators.put(IndexFactory.IndexType.INTERVAL_TREE,treeSmallBin);
        }

        return creators;
    }


    public void addFeature(final Feature f, final long filePosition) {
        // protected static Map<Double,Index> createIndex(FileBasedFeatureIterator<Feature> iterator, Map<IndexType,IndexCreator> creators, IndexBalanceApproach iba) {
        // feed each feature to the indexes we've created
        // first take care of the stats
        featureCount++;

        // calculate the number of bases seen - we have to watch out for the situation where the last record was on the previous chromosome
        basesSeen = (lastFeature == null) ? basesSeen + f.getStart() :
                ((f.getStart() - lastFeature.getStart() >= 0) ? basesSeen + (f.getStart() - lastFeature.getStart()) : basesSeen + f.getStart());

        longestFeatureLength = Math.max(longestFeatureLength,(f.getEnd()-f.getStart()) + 1);

        // push the longest feature to the running stats
        stats.push(longestFeatureLength);

        // now feed the feature to each of our creators
        for (final IndexCreator creator : creators.values()) {
            creator.addFeature(f,filePosition);
        }

        //Redundant check, done in IndexFactory
        // if the last feature is after the current feature, exception out
//        if (lastFeature != null && f.getStart() < lastFeature.getStart() && lastFeature.getChr().equals(f.getChr()))
//            throw new TribbleException.MalformedFeatureFile("We saw a record with a start of " + f.getChr() + ":" + f.getStart() +
//                    " after a record with a start of " + lastFeature.getChr() + ":" + lastFeature.getStart(), inputFile.getAbsolutePath());

        // save the last feature
        lastFeature = f;
    }

    /**
     * score the available indexes for the specified density and feature lengths
     *
     * The scoring method is trying to determine how many features would be returned for a sample one base query; or:
     * (features/seek).  For the interval index this is clear: it's the bin size (interval is binned by feature count).
     * for Linear indexes it's the density of features X the number of bins we need to retrieve (which is determined
     * by the bin size X the longest feature).
     *
     * @param densityOfFeatures the density of features (features/base)
     * @param indexes Map from IndexType -> IndexCreator
     * @param longestFeature the longest feature we've found
     * @param iba the index balancing approach
     * @return the best index available for the target indexes
     */
    protected static LinkedHashMap<Double,TribbleIndexCreator> scoreIndexes(final double densityOfFeatures, final Map<IndexFactory.IndexType,TribbleIndexCreator> indexes, final int longestFeature, final IndexFactory.IndexBalanceApproach iba) {
        if (indexes.size() < 1) throw new IllegalArgumentException("Please specify at least one index to evaluate");

        final LinkedHashMap<Double,TribbleIndexCreator> scores = new LinkedHashMap<Double,TribbleIndexCreator>();

        for (final Map.Entry<IndexFactory.IndexType,TribbleIndexCreator> entry : indexes.entrySet()) {
            // we have different scoring
            if (entry.getValue() instanceof LinearIndexCreator) {
                final double binSize = ((LinearIndexCreator)(entry.getValue())).getBinSize();
                scores.put(binSize * densityOfFeatures * Math.ceil((double) longestFeature / binSize), entry.getValue());
            } else if (entry.getValue() instanceof IntervalIndexCreator) {
                scores.put((double) ((IntervalIndexCreator)entry.getValue()).getFeaturesPerInterval(), entry.getValue());
            } else {
                throw new TribbleException.UnableToCreateCorrectIndexType("Unknown index type, we don't have a scoring method for " + entry.getValue().getClass());
            }
        }
        return scores;
    }

    /**
     * utility function to find the min of a list
     * @param scores the list of scaled features/bin scores for each index type
     * @return the best score <b>index value</b>
     */
    private TribbleIndexCreator getMinIndex(final Map<Double,TribbleIndexCreator> scores, final IndexFactory.IndexBalanceApproach iba) {
        final TreeMap<Double,TribbleIndexCreator> map = new TreeMap<Double,TribbleIndexCreator>();
        map.putAll(scores);
        
        // if we are optimizing for seek time, choose the lowest score (adjusted features/bin value), if for storage size, choose the opposite
        final TribbleIndexCreator idx = (iba != IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME) ? map.get(map.lastKey()) : map.get(map.firstKey());
        return idx;
    }

    @Override
    public void addProperty(final String key, final String value) {
        for (final TribbleIndexCreator creator : creators.values()) {
            creator.addProperty(key, value);
        }
    }
}
