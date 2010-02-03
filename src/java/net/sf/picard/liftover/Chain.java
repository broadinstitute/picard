/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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
package net.sf.picard.liftover;

import net.sf.picard.PicardException;
import net.sf.picard.util.Interval;
import net.sf.picard.util.OverlapDetector;
import net.sf.samtools.util.AsciiLineReader;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.regex.Pattern;

/**
 * Holds a single chain from a UCSC chain file.  In a chain file, a chain consists of a header line followed by
 * alignment data lines.  Chain class embodies the header line, and the list of ContinuousBlocks embodies the
 * alignment data lines.
 *
 * In UCSC liftOver terminology, the "target" is the "from" genome build, and the "query" is the "to" genome build.
 * E.g. when mapping from HG18 to HG19, the HG18 coordinates are "target" and HG19 is "query."  The terminology
 * confuses me but I'm trying to be faithful to the C liftOver implementation.
 *
 * Chain coordinates are zero-based, half open.  However, there is also an Interval attribute of Chain that is in
 * standard Picard coordinates, i.e. one-based inclusive.
 *
 * @author alecw@broadinstitute.org
 */
class Chain {
    // For parsing chain file
    private static final Pattern SPLITTER = Pattern.compile("\\s");

    /** on-based, inclusive, so that Chain can be stored in an OverlapDetector */
    final Interval interval;
    /** Total score for chain. */
    final double score;
    /** target name. */
    final String tName;
    /** Overall size of target. */
    final int tSize;
    /* tStrand always +, so not stored */
    /** Start of range covered in target. */
    final int tStart;
    /** End of range covered in target. */
    final int tEnd;
    /** query name. */
    final String qName;
    /** Overall size of query. */
    final int qSize;
    /** Query strand. */
    final boolean qNegativeStrand;
    /** Start of range covered in query. */
    final int qStart;
    /** End of range covered in query. */
    final int qEnd;
    /** ID of chain in file. */
    final int id;
    private final List<ContinuousBlock> blockList = new ArrayList<ContinuousBlock>();

    /**
     * Construct a Chain from the parsed header fields.
     */
    private Chain(final double score, final String tName, final int tSize, final int tStart, final int tEnd,
          final String qName, final int qSize, final boolean qNegativeStrand,
          final int qStart, final int qEnd, final int id) {
        // Convert  to one-based, inclusive for Interval.
        interval = new Interval(tName, tStart + 1, tEnd);
        this.id = id;
        this.qEnd = qEnd;
        this.qName = qName;
        this.qNegativeStrand = qNegativeStrand;
        this.qSize = qSize;
        this.qStart = qStart;
        this.score = score;
        this.tEnd = tEnd;
        this.tName = tName;
        this.tSize = tSize;
        this.tStart = tStart;
    }


    /**
     * Holds a range that continuously lines up between target and query genome builds.
     * Indices are 0-based, half-open.
     */
    static class ContinuousBlock {
        final int tStart,tEnd;		/* Range covered in target. */
        final int qStart,qEnd;		/* Range covered in query. */
        //int score;	 	 	/* Score of block. */

        private ContinuousBlock(final int tStart, final int tEnd, final int qStart, final int qEnd) {
            this.tStart = tStart;
            this.tEnd = tEnd;
            this.qStart = qStart;
            this.qEnd = qEnd;
        }
    }

    private void addBlock(final int tStart, final int tEnd, final int qStart, final int qEnd) {
        blockList.add(new ContinuousBlock(tStart, tEnd, qStart, qEnd));
    }

    /**
     * @return The ith ContinuousBlock in this Chain.
     */
    ContinuousBlock getBlock(final int i) {
        return blockList.get(i);
    }

    /**
     * @return Unmodifiable list of ContinuousBlocks in this Chain.
     */
    List<ContinuousBlock> getBlocks() {
        return Collections.unmodifiableList(blockList);
    }

    /**
     * Read all the chains and load into an OverlapDetector.
     * @param reader Text representation of chains.
     * @param chainFile For error messages only.
     * @return OverlapDetector will all Chains from reader loaded into it.
     */
    static OverlapDetector<Chain> loadChains(final AsciiLineReader reader, final String chainFile) {
        final OverlapDetector<Chain> ret = new OverlapDetector<Chain>(0, 0);
        Chain chain;
        while ((chain = Chain.loadChain(reader, chainFile)) != null) {
            ret.addLhs(chain, chain.interval);
        }
        return ret;
    }

    /**
     * Read a single Chain from reader.
     * @param reader Text representation of chains.
     * @param chainFile For error messages only.
     * @return New Chain with associated ContinuousBlocks.
     */
    private static Chain loadChain(final AsciiLineReader reader, final String chainFile) {
        String line = reader.readLine();
        if (line == null) {
            return null;
        }
        String[] chainFields = SPLITTER.split(line);
        if (chainFields.length != 13) {
            throwChainFileParseException("chain line has wrong number of fields", chainFile, reader.getLineNumber());
        }
        if (!"chain".equals(chainFields[0])) {
            throwChainFileParseException("chain line does not start with 'chain'", chainFile, reader.getLineNumber());
        }
        double score = 0;
        String tName = null;
        int tSize = 0;
        int tStart = 0;
        int tEnd = 0;
        String qName = null;
        int qSize = 0;
        boolean qNegativeStrand = false;
        int qStart = 0;
        int qEnd = 0;
        int id = 0;
        try {
            score = Double.parseDouble(chainFields[1]);
            tName = chainFields[2];
            tSize = Integer.parseInt(chainFields[3]);
            // Strand ignored because it is always +
            tStart = Integer.parseInt(chainFields[5]);
            tEnd = Integer.parseInt(chainFields[6]);
            qName = chainFields[7];
            qSize = Integer.parseInt(chainFields[8]);
            qNegativeStrand = chainFields[9].equals("-");
            qStart = Integer.parseInt(chainFields[10]);
            qEnd = Integer.parseInt(chainFields[11]);
            id = Integer.parseInt(chainFields[12]);
        } catch (NumberFormatException e) {
            throwChainFileParseException("Invalid field", chainFile, reader.getLineNumber());
        }
        final Chain chain = new Chain(score, tName, tSize, tStart, tEnd, qName, qSize, qNegativeStrand, qStart,
                qEnd, id);
        int q = chain.qStart;
        int t = chain.tStart;
        boolean sawLastLine = false;
        while (true) {
            line = reader.readLine();
            if (line == null || line.equals("")) {
                if (!sawLastLine) {
                    throwChainFileParseException("Reached end of chain without seeing terminal block", chainFile, reader.getLineNumber());
                }
                break;
            }
            if (sawLastLine) {
                throwChainFileParseException("Terminal block seen before end of chain", chainFile, reader.getLineNumber());
            }
            String[] blockFields = SPLITTER.split(line);
            if (blockFields.length == 1) {
                sawLastLine = true;
            } else if (blockFields.length != 3) {
                throwChainFileParseException("Block line has unexpected number of fields", chainFile, reader.getLineNumber());
            }
            int size = Integer.parseInt(blockFields[0]);
            chain.addBlock(t, t+size, q, q+size);
            if (sawLastLine) {
                q = Integer.MAX_VALUE;
                t = Integer.MAX_VALUE;
            } else {
                t += Integer.parseInt(blockFields[1]) + size;
                q += Integer.parseInt(blockFields[2]) + size;
            }

        }
        return chain;
    }

    private static void throwChainFileParseException(final String message, final String chainFile, final int lineNumber) {
        throw new PicardException(message + " in chain file " + chainFile + " at line " + lineNumber);
    }
}
