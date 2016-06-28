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
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Class that is used to represent the results of comparing a read group within a SAM file
 * against one or more set of fingerprint genotypes.
 *
 * @author Tim Fennell
 */
public class FingerprintResults {
    private final File samFile;
    private final String readGroup;
    private final SortedSet<MatchResults> matchResults = new TreeSet<>();

    public FingerprintResults(final File samFile, final String readGroup) {
        this.samFile = samFile;
        this.readGroup = readGroup;
    }

    public void addResults(final MatchResults matchResults) {
        this.matchResults.add(matchResults);
    }

    public File getSamFile() { return samFile; }
    public String getReadGroup() { return readGroup; }
    public SortedSet<MatchResults> getMatchResults() { return matchResults; }
}

