/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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

package picard.fingerprint;

import java.io.File;
import java.util.*;

/**
 * Small class to represent a genetic fingerprint as a set of HaplotypeProbabilities
 * objects that give the relative probabilities of each of the possible haplotypes
 * at a locus.
 *
 * @author Tim Fennell
 */
public class Fingerprint extends TreeMap<HaplotypeBlock, HaplotypeProbabilities> {
    private final String sample;
    private final File source;
    private final String info;

    public Fingerprint(final String sample, final File source, final String info) {
        this.sample = sample;
        this.source = source;
        this.info = info;
    }

    public String getSample() { return sample; }

    public File getSource() { return source; }

    public String getInfo() { return info; }

    public String getPrintableId() {
        return getSample() + "@" + (source == null ? "" : source.getName()) + (info == null ? "" : (":" + info));
    }

    public void add(final HaplotypeProbabilities h) {
        put(h.getHaplotype(), h);
    }

    /**
     * Merges the likelihoods from the supplied Fingerprint into the likelihoods for this fingerprint.
     */
    public void merge(final Fingerprint other) {
        final Set<HaplotypeBlock> haps = new HashSet<>();
        haps.addAll(keySet());
        haps.addAll(other.keySet());

        for (final HaplotypeBlock haplotype : haps) {
            HaplotypeProbabilities probabilities = get(haplotype);
            final HaplotypeProbabilities otherProbabilities = other.get(haplotype);
            if (probabilities == null) {
                probabilities = otherProbabilities;
                put(haplotype, probabilities);
            } else if (otherProbabilities != null) {
                probabilities.merge(otherProbabilities);
            }
        }
    }

    /**
     * Attempts to filter out haplotypes that may have suspect genotyping by removing haplotypes that reach
     * a minimum confidence score yet have a significant fraction of observations from a third or fourth allele.
     */
    public void filterSuspectSites() {
        final Iterator<Map.Entry<HaplotypeBlock, HaplotypeProbabilities>> iterator = entrySet().iterator();
        while (iterator.hasNext()) {
            final Map.Entry<HaplotypeBlock, HaplotypeProbabilities> entry = iterator.next();
            final HaplotypeProbabilities p = entry.getValue();
            if (p instanceof HaplotypeProbabilitiesFromSequence) {
                final HaplotypeProbabilitiesFromSequence probs = (HaplotypeProbabilitiesFromSequence) p;

                if (probs.getLodMostProbableGenotype() >= 3 && probs.getFractionUnexpectedAlleleObs() > 0.1) {
                    iterator.remove();
                }
            }
        }
    }
}
