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
import java.util.Collection;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Represents the results of a fingerprint comparison between one dataset and a specific
 * fingerprint file.  Implements Comparable so that better matches (higher positive LODs)
 * are sorted earlier.
 *
 * @author Tim Fennell
 */
public class MatchResults implements Comparable<MatchResults> {
    private final File fingerprintFile;
    private final String sample;
    private final double sampleLikelihood;
    private final double populationLikelihood;
    private final double LOD;
    //the lod score when assuming the left sample is tumor and the right is normal
    private final double lodTN;
    //the lod score when assuming the left sample is tumor and the right is normal
    private final double lodNT;

    public double getLodNT() {
        return lodNT;
    }

    public double getLodTN() {
        return lodTN;
    }

    private final SortedSet<LocusResult> locusResults = new TreeSet<LocusResult>();

    MatchResults(final File fingerprintFile, final String sample,
                 final double sampleLikelihood, final double populationLikelihood, final double lodTN, final double lodNT,
                 final Collection<LocusResult> locusResults) {
        this.fingerprintFile = fingerprintFile;
        this.sample = sample;
        this.sampleLikelihood = sampleLikelihood;
        this.populationLikelihood = populationLikelihood;
        this.LOD = sampleLikelihood - populationLikelihood;
        this.lodTN = lodTN;
        this.lodNT = lodNT;

        if (locusResults != null) {
            this.locusResults.addAll(locusResults);
        }
    }

    public void addLocusResult(final LocusResult result) {
        this.locusResults.add(result);
    }

    /** Provides a natural sort so that better matches (by LOD) sort earlier. */
    @Override public int compareTo(final MatchResults that) {
        if (this.LOD != that.LOD) {
            return this.LOD > that.LOD ? -1 : 1;
        }
        else {
            return this.sample.compareTo(that.sample);
        }
    }

    public String getSample() { return sample; }
    public double getSampleLikelihood() { return sampleLikelihood; }
    public double getPopulationLikelihood() { return populationLikelihood; }
    public double getLOD() { return LOD; }
    public SortedSet<LocusResult> getLocusResults() { return locusResults; }
    public File getFingerprintFile() { return this.fingerprintFile; }
}
