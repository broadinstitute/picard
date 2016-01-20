/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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

package picard.analysis;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import picard.util.IlluminaUtil;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * A utility class for matching reads to adapters.
 * Note that this is different from ClippingUtility in that it tries to match the starts of reads
 * to any part of the adapter (as opposed to finding the start of the adapter anywhere in the read).
 */
public class AdapterUtility {

    //The number of bases to check in order to map a read to an adapter
    private static final int ADAPTER_MATCH_LENGTH = 16;

    // The maximum number of mismatches a read can have and still be considered as matching an adapter
    private static final int MAX_ADAPTER_ERRORS = 1;

    // byte arrays in both fwd and rc for the adapter sequences
    final byte [][] adapterKmers;

    public static List<String> DEFAULT_ADAPTER_SEQUENCE = CollectionUtil.makeList(
            IlluminaUtil.IlluminaAdapterPair.SINGLE_END.get5PrimeAdapter(),
            IlluminaUtil.IlluminaAdapterPair.SINGLE_END.get3PrimeAdapter(),
            IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get5PrimeAdapter(),
            IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get3PrimeAdapter(),
            IlluminaUtil.IlluminaAdapterPair.INDEXED.get5PrimeAdapter(),
            IlluminaUtil.IlluminaAdapterPair.INDEXED.get3PrimeAdapter()
    );
    // TODO -- consider adding DUAL_INDEXED to the list above

    public AdapterUtility(final List<String> adapterSequence) {
        adapterKmers = prepareAdapterSequences(adapterSequence);
    }

    /** Converts the supplied adapter sequences to byte arrays in both fwd and rc */
    private static byte [][] prepareAdapterSequences(final List<String> adapterSequence) {
        final Set<String> kmers = new HashSet<>();

        // Make a set of all kmers of adapterMatchLength
        for (final String seq : adapterSequence) {
            for (int i=0; i<=seq.length() - ADAPTER_MATCH_LENGTH; ++i) {
                final String kmer = seq.substring(i, i+ADAPTER_MATCH_LENGTH).toUpperCase();

                int ns = 0;
                for (final char ch : kmer.toCharArray()) if (ch == 'N') ++ns;
                if (ns <= MAX_ADAPTER_ERRORS) {
                    kmers.add(kmer);
                    kmers.add(SequenceUtil.reverseComplement(kmer));
                }
            }
        }

        // Make an array of byte[] for the kmers
        final byte [][] adapterKmers = new byte[kmers.size()][];
        int i=0;
        for (final String kmer : kmers) {
            adapterKmers[i++] = StringUtil.stringToBytes(kmer);
        }
        return adapterKmers;
    }

    /**
     * Checks the first ADAPTER_MATCH_LENGTH bases of the read against known adapter sequences and returns
     * true if the read matches an adapter sequence with MAX_ADAPTER_ERRORS mismsatches or fewer.
     *
     * @param read the basecalls for the read in the order and orientation the machine read them
     * @return true if the read matches an adapter and false otherwise
     */
    public boolean isAdapterSequence(final byte[] read) {
        if (read.length < ADAPTER_MATCH_LENGTH) return false;

        for (final byte[] adapter : adapterKmers) {
            int errors = 0;

            for (int i=0; i<adapter.length; ++i) {
                if (read[i] != adapter[i]) {
                    if (++errors > MAX_ADAPTER_ERRORS) break;
                }
            }

            if (errors <= MAX_ADAPTER_ERRORS) return true;
        }

        return false;
    }
}
