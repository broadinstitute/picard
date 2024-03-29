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
package picard.util;

import htsjdk.samtools.ReservedTagConstants;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.samtools.util.Tuple;

import java.util.*;
import java.util.concurrent.atomic.AtomicReference;
import java.util.stream.Stream;

/**
 * Store one or more AdapterPairs to use to mark adapter sequence of SAMRecords.  This is a very compute-intensive process, so
 * this class implements two heuristics to reduce computation:
 * - Adapter sequences are truncated, and then any adapter pairs that become identical after truncation are collapsed into a single pair.
 * - After a specified number of reads with adapter sequence has been seen, prune the list of adapter pairs to include only the most
 *   frequently seen adapters.  For a flowcell, there should only be a single adapter pair found.
 *
 * Note that the AdapterPair object returned by all the adapterTrim* methods will not be one of the original AdapterPairs
 * passed to the ctor, but rather will be one of the truncated copies.
 */
public class AdapterMarker {
    public static final int DEFAULT_ADAPTER_LENGTH = 30;
    public static final int DEFAULT_PRUNE_ADAPTER_LIST_AFTER_THIS_MANY_ADAPTERS_SEEN = 100;
    public static final int DEFAULT_NUM_ADAPTERS_TO_KEEP = 1;

    // It is assumed that these are set once during execution, before the class is used to mark any adapters, but this is not enforced.
    private int thresholdForSelectingAdaptersToKeep = DEFAULT_PRUNE_ADAPTER_LIST_AFTER_THIS_MANY_ADAPTERS_SEEN;
    private int numAdaptersToKeep = DEFAULT_NUM_ADAPTERS_TO_KEEP;
    private int minSingleEndMatchBases = ClippingUtility.MIN_MATCH_BASES;
    private int minPairMatchBases = ClippingUtility.MIN_MATCH_PE_BASES;
    private double maxSingleEndErrorRate = ClippingUtility.MAX_ERROR_RATE;
    private double maxPairErrorRate = ClippingUtility.MAX_PE_ERROR_RATE;

    // This is AtomicReference because one thread could be matching adapters while the threshold has been crossed in another
    // thread and the array is being replaced.
    private final AtomicReference<AdapterPair[]> adapters = new AtomicReference<>();

    // All the members below are only accessed within a synchronized block.
    private boolean thresholdReached = false;
    private int numAdaptersSeen = 0;
    private final CollectionUtil.DefaultingMap<AdapterPair, Integer> seenCounts = new CollectionUtil.DefaultingMap<AdapterPair, Integer>(0);

    //Store all the sam records we have seen prior to choosing an adapter so that we can go back and fix the ones
    //that have clipping tags for adapters that were not chosen.
    private Map<AdapterPair, List<SAMRecord>> preAdapterPrunedRecords = new HashMap<>();

    /**
     * Truncates adapters to DEFAULT_ADAPTER_LENGTH
     *
     * @param originalAdapters These should be in order from longest & most likely to shortest & least likely.
     */
    public AdapterMarker(final AdapterPair... originalAdapters) {
        this(DEFAULT_ADAPTER_LENGTH, originalAdapters);
    }

    /**
     * @param adapterLength    Truncate adapters to this length.
     * @param originalAdapters These should be in order from longest & most likely to shortest & least likely.
     */
    public AdapterMarker(final int adapterLength, final AdapterPair... originalAdapters) {
        // Truncate each AdapterPair to the given length, and then combine any that end up the same after truncation.
        final ArrayList<TruncatedAdapterPair> truncatedAdapters = new ArrayList<TruncatedAdapterPair>();
        for (final AdapterPair adapter : originalAdapters) {
            final TruncatedAdapterPair truncatedAdapter = makeTruncatedAdapterPair(adapter, adapterLength);
            final int matchingIndex = truncatedAdapters.indexOf(truncatedAdapter);
            if (matchingIndex == -1) {
                truncatedAdapters.add(truncatedAdapter);
            } else {
                final TruncatedAdapterPair matchingAdapter = truncatedAdapters.get(matchingIndex);
                matchingAdapter.setName(matchingAdapter.getName() + "|" + adapter.getName());
            }
        }
        adapters.set(truncatedAdapters.toArray(new AdapterPair[truncatedAdapters.size()]));
    }

    public int getNumAdaptersToKeep() {
        return numAdaptersToKeep;
    }

    /**
     * After seeing the thresholdForSelectingAdapters number of adapters, keep up to this many of the original adapters.
     */
    public synchronized AdapterMarker setNumAdaptersToKeep(final int numAdaptersToKeep) {
        if (numAdaptersToKeep <= 0) {
            throw new IllegalArgumentException(String.format("numAdaptersToKeep should be positive: %d", numAdaptersToKeep));
        }
        this.numAdaptersToKeep = numAdaptersToKeep;
        return this;
    }

    public int getThresholdForSelectingAdaptersToKeep() {
        return thresholdForSelectingAdaptersToKeep;
    }

    /**
     * When this number of adapters have been matched, discard the least-frequently matching ones.
     *
     * @param thresholdForSelectingAdaptersToKeep set to -1 to never discard any adapters.
     */
    public synchronized AdapterMarker setThresholdForSelectingAdaptersToKeep(final int thresholdForSelectingAdaptersToKeep) {
        this.thresholdForSelectingAdaptersToKeep = thresholdForSelectingAdaptersToKeep;
        return this;
    }

    public int getMinSingleEndMatchBases() {
        return minSingleEndMatchBases;
    }

    /**
     * @param minSingleEndMatchBases When marking a single-end read, adapter must match at least this many bases.
     */
    public synchronized AdapterMarker setMinSingleEndMatchBases(final int minSingleEndMatchBases) {
        this.minSingleEndMatchBases = minSingleEndMatchBases;
        return this;
    }

    public int getMinPairMatchBases() {
        return minPairMatchBases;
    }

    /**
     * @param minPairMatchBases When marking a paired-end read, adapter must match at least this many bases.
     */
    public synchronized AdapterMarker setMinPairMatchBases(final int minPairMatchBases) {
        this.minPairMatchBases = minPairMatchBases;
        return this;
    }

    public double getMaxSingleEndErrorRate() {
        return maxSingleEndErrorRate;
    }

    /**
     * @param maxSingleEndErrorRate For single-end read, no more than this fraction of the bases that align with the adapter can
     *                              mismatch the adapter and still be considered an adapter match.
     */
    public synchronized AdapterMarker setMaxSingleEndErrorRate(final double maxSingleEndErrorRate) {
        this.maxSingleEndErrorRate = maxSingleEndErrorRate;
        return this;
    }

    public double getMaxPairErrorRate() {
        return maxPairErrorRate;
    }

    /**
     * @param maxPairErrorRate For paired-end read, no more than this fraction of the bases that align with the adapter can
     *                         mismatch the adapter and still be considered an adapter match.
     */
    public synchronized AdapterMarker setMaxPairErrorRate(final double maxPairErrorRate) {
        this.maxPairErrorRate = maxPairErrorRate;
        return this;
    }

    public AdapterPair adapterTrimIlluminaSingleRead(final SAMRecord read) {
        return adapterTrimIlluminaSingleRead(read, minSingleEndMatchBases, maxSingleEndErrorRate);
    }

    /**
     * Return the adapter to be trimmed from a read represented as an array of bytes[]
     * @param read The byte array of read data
     * @param templateIndex The paired index of the reads (1 or 2, 1 for single ended reads)
     * @return The adapter pair that matched the read and its index in the read.
     */
    public Tuple<AdapterPair, Integer> findAdapterPairAndIndexForSingleRead(final byte[] read, final int templateIndex) {
        return findAdapterPairAndIndexForSingleRead(read, minSingleEndMatchBases, maxSingleEndErrorRate, templateIndex);
    }

    public AdapterPair adapterTrimIlluminaPairedReads(final SAMRecord read1, final SAMRecord read2) {
        return adapterTrimIlluminaPairedReads(read1, read2, minPairMatchBases, maxPairErrorRate);
    }

    /**
     * Overrides defaults for minMatchBases and maxErrorRate
     */
    public AdapterPair adapterTrimIlluminaSingleRead(final SAMRecord read, final int minMatchBases, final double maxErrorRate) {
        final AdapterPair ret = ClippingUtility.adapterTrimIlluminaSingleRead(read, minMatchBases, maxErrorRate, adapters.get());
        tallyAndFixAdapters(ret, read);
        return ret;
    }

    /**
     * Return the adapter to be trimmed from a read represented as an array of bytes[]
     * @param read The byte array of read data
     * @param minMatchBases The minimum number of base matches required for adapter matching
     * @param maxErrorRate The maximum error rate allowed for adapter matching
     * @param templateIndex The paired index of the reads (1 or 2, 1 for single ended reads)
     * @return The adapter pair that matched the read and its index in the read or null.
     */
    public Tuple<AdapterPair, Integer> findAdapterPairAndIndexForSingleRead(final byte[] read,
                                                    final int minMatchBases,
                                                    final double maxErrorRate,
                                                    int templateIndex) {
        final Tuple<AdapterPair, Integer> ret = ClippingUtility.findAdapterPairAndIndexForSingleRead(read,
                minMatchBases,
                maxErrorRate,
                templateIndex,
                adapters.get());
        if (ret != null && !thresholdReached) {
            tallyFoundAdapter(ret.a, false);
        }
        return ret;
    }

    private void tallyAndFixAdapters(AdapterPair ret, SAMRecord... reads) {
        if (ret != null && !thresholdReached) {
            if (!preAdapterPrunedRecords.containsKey(ret)) {
                preAdapterPrunedRecords.put(ret, new ArrayList<>());
            }
            Arrays.stream(reads).forEach(read -> preAdapterPrunedRecords.get(ret).add(read));
            tallyFoundAdapter(ret, true);
        }
    }

    /**
     * Overrides defaults for minMatchBases and maxErrorRate
     */
    public AdapterPair adapterTrimIlluminaPairedReads(final SAMRecord read1, final SAMRecord read2,
                                                      final int minMatchBases, final double maxErrorRate) {
        final AdapterPair ret = ClippingUtility.adapterTrimIlluminaPairedReads(read1, read2, minMatchBases, maxErrorRate, adapters.get());
        tallyAndFixAdapters(ret, read1, read2);
        return ret;
    }

    /**
     * For unit testing only
     */
    AdapterPair[] getAdapters() {
        return adapters.get();
    }

    private TruncatedAdapterPair makeTruncatedAdapterPair(final AdapterPair adapterPair, final int adapterLength) {
        return new TruncatedAdapterPair("truncated " + adapterPair.getName(),
                substringAndRemoveTrailingNs(adapterPair.get3PrimeAdapterInReadOrder(), adapterLength),
                substringAndRemoveTrailingNs(adapterPair.get5PrimeAdapterInReadOrder(), adapterLength));
    }

    /**
     * Truncate to the given length, and in addition truncate any trailing Ns.
     */
    private String substringAndRemoveTrailingNs(final String s, int length) {
        length = Math.min(length, s.length());
        final byte[] bytes = StringUtil.stringToBytes(s);
        while (length > 0 && SequenceUtil.isNoCall(bytes[length - 1])) {
            length--;
        }
        return s.substring(0, length);
    }

    /**
     * Keep track of every time an adapter is found, until it is time to prune the list of adapters.
     */
    private void tallyFoundAdapter(final AdapterPair foundAdapter, final boolean fixPreviousRecords) {
        // If caller does not want adapter pruning, do nothing.
        if (thresholdForSelectingAdaptersToKeep < 1) return;
        synchronized (this) {

            // Tally this adapter
            seenCounts.put(foundAdapter, seenCounts.get(foundAdapter) + 1);

            // Keep track of the number of times an adapter has been seen.
            numAdaptersSeen += 1;

            // Reached the threshold for pruning the list.
            if (numAdaptersSeen >= thresholdForSelectingAdaptersToKeep) {

                // Sort adapters by number of times each has been seen.
                final TreeMap<Integer, AdapterPair> sortedAdapters = new TreeMap<Integer, AdapterPair>(new Comparator<Integer>() {
                    @Override
                    public int compare(final Integer integer, final Integer integer2) {
                        // Reverse of natural ordering
                        return integer2.compareTo(integer);
                    }
                });
                for (final Map.Entry<AdapterPair, Integer> entry : seenCounts.entrySet()) {
                    sortedAdapters.put(entry.getValue(), entry.getKey());
                }

                // Keep the #numAdaptersToKeep adapters that have been seen the most, plus any ties.
                final ArrayList<AdapterPair> bestAdapters = new ArrayList<AdapterPair>(numAdaptersToKeep);
                int countOfLastAdapter = Integer.MAX_VALUE;
                for (final Map.Entry<Integer, AdapterPair> entry : sortedAdapters.entrySet()) {
                    if (bestAdapters.size() >= numAdaptersToKeep) {
                        if (entry.getKey() == countOfLastAdapter) {
                            bestAdapters.add(entry.getValue());
                        } else {
                            break;
                        }
                    } else {
                        countOfLastAdapter = entry.getKey();
                        bestAdapters.add(entry.getValue());
                    }
                }
                // Replace the existing list with the pruned list.
                thresholdReached = true;
                adapters.set(bestAdapters.toArray(new AdapterPair[bestAdapters.size()]));
                if (fixPreviousRecords) {
                    fixAlreadySeenReads();
                }
            }
        }
    }

    private void fixAlreadySeenReads() {
        //remove all the reads for the selected adapters
        Arrays.stream(adapters.get()).forEach(adapter -> preAdapterPrunedRecords.remove(adapter));
        //anything left is marked with the incorrect adapter and needs its XT tag removed
        preAdapterPrunedRecords.values().forEach(readList -> readList.parallelStream().forEach(read -> {
            Stream<SAMRecord.SAMTagAndValue> filterAttributes = read.getAttributes().stream().filter(tag -> !tag.tag.equals(ReservedTagConstants.XT));
            read.clearAttributes();
            filterAttributes.forEach(tag -> read.setAttribute(tag.tag, tag.value));
        }));
        preAdapterPrunedRecords.clear();
    }

    private static final class TruncatedAdapterPair implements AdapterPair {
        String name;
        final String fivePrime, threePrime, fivePrimeReadOrder;
        final byte[] fivePrimeBytes, threePrimeBytes, fivePrimeReadOrderBytes;

        private TruncatedAdapterPair(final String name, final String threePrimeReadOrder, final String fivePrimeReadOrder) {
            this.name = name;
            this.threePrime = threePrimeReadOrder;
            this.threePrimeBytes = StringUtil.stringToBytes(threePrimeReadOrder);
            this.fivePrimeReadOrder = fivePrimeReadOrder;
            this.fivePrimeReadOrderBytes = StringUtil.stringToBytes(fivePrimeReadOrder);
            this.fivePrime = SequenceUtil.reverseComplement(fivePrimeReadOrder);
            this.fivePrimeBytes = StringUtil.stringToBytes(this.fivePrime);
        }

        public String get3PrimeAdapter() {
            return threePrime;
        }

        public String get5PrimeAdapter() {
            return fivePrime;
        }

        public String get3PrimeAdapterInReadOrder() {
            return threePrime;
        }

        public String get5PrimeAdapterInReadOrder() {
            return fivePrimeReadOrder;
        }

        public byte[] get3PrimeAdapterBytes() {
            return threePrimeBytes;
        }

        public byte[] get5PrimeAdapterBytes() {
            return fivePrimeBytes;
        }

        public byte[] get3PrimeAdapterBytesInReadOrder() {
            return threePrimeBytes;
        }

        public byte[] get5PrimeAdapterBytesInReadOrder() {
            return fivePrimeReadOrderBytes;
        }

        public String getName() {
            return this.name;
        }

        public void setName(final String name) {
            this.name = name;
        }

        // WARNING: These methods ignore the name member!
        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final TruncatedAdapterPair that = (TruncatedAdapterPair) o;

            if (!fivePrime.equals(that.fivePrime)) return false;
            if (!threePrime.equals(that.threePrime)) return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result = fivePrime.hashCode();
            result = 31 * result + threePrime.hashCode();
            return result;
        }

        @Override
        public String toString() {
            return "TruncatedAdapterPair{" +
                    "fivePrimeReadOrder='" + fivePrimeReadOrder + '\'' +
                    ", threePrime='" + threePrime + '\'' +
                    ", name='" + name + '\'' +
                    '}';
        }
    }
}
