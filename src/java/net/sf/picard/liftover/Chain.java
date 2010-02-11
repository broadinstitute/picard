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
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Interval;
import net.sf.picard.util.OverlapDetector;
import net.sf.samtools.util.AsciiLineReader;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.regex.Pattern;

/**
 * Holds a single chain from a UCSC chain file.  Chain file format is described here: http://genome.ucsc.edu/goldenPath/help/chain.html
 *
 * In a chain file, a chain consists of a header line followed by alignment data lines.  Chain class embodies the header
 * line, and the list of ContinuousBlocks embodies the alignment data lines.
 *
 * A continuous block represents a continuous range of the "from" genome that maps to a continuous range of the "to"
 * genome of the same length.  A chain is an ordered list of continuous blocks, with gaps between the continuous blocks.
 * All the continuous blocks in a chain must map from a single "from" sequence to a single "to" sequence.  All the
 * continuous blocks in a chain map from the positive strand in the "from" genome build to the same strand in the
 * "to" genome build.  The gaps in between the continuous blocks in a chain represent un-lift-overable regions.
 * A gap in a chain may be found in another chain (e.g. if a portion of a sequence is reversed in the "to" genome).
 *
 * In UCSC liftOver terminology, the "target" is the "from" genome build, and the "query" is the "to" genome build.
 * E.g. when mapping from HG18 to HG19, the HG18 coordinates are "target" and HG19 is "query."  The USCS terminology
 * is not used here because it confuses me.
 *
 * Chain coordinates are zero-based, half open.  However, there is also an Interval attribute of Chain that is in
 * standard Picard coordinates, i.e. one-based inclusive.
 *
 * @author alecw@broadinstitute.org
 */
class Chain {
    // For parsing chain file
    private static final Pattern SPLITTER = Pattern.compile("\\s");

    /** one-based, inclusive, so that Chain can be stored in an OverlapDetector */
    final Interval interval;
    /** Total score for chain is not used in basic liftover so not stored. */
    // final double score;
    final String fromSequenceName;
    /** Overall size of the "from" sequence. */
    final int fromSequenceSize;
    /* tStrand always +, so not stored */
    /** Start of range covered in "from" sequence. */
    final int fromChainStart;
    /** End of range covered in "from" sequence. */
    final int fromChainEnd;
    final String toSequenceName;
    /** Overall size of the "to" sequence. */
    final int toSequenceSize;
    /** "to" strand. If this is true, then the region covered by this chain is flipped in the "to" genome.  */
    final boolean toNegativeStrand;
    /** Start of range covered in "to" sequence. */
    final int toChainStart;
    /** End of range covered in "to" sequence. */
    final int toChainEnd;
    /** ID of chain in file.  */
    final int id;
    private final List<ContinuousBlock> blockList = new ArrayList<ContinuousBlock>();

    /**
     * Construct a Chain from the parsed header fields.
     */
    private Chain(final double score, final String fromSequenceName, final int fromSequenceSize, final int fromChainStart, final int fromChainEnd,
          final String toSequenceName, final int toSequenceSize, final boolean toNegativeStrand,
          final int toChainStart, final int toChainEnd, final int id) {
        // Convert  to one-based, inclusive for Interval.
        interval = new Interval(fromSequenceName, fromChainStart + 1, fromChainEnd);
        // Not used
        //this.id = id;
        this.toChainEnd = toChainEnd;
        this.toSequenceName = toSequenceName;
        this.toNegativeStrand = toNegativeStrand;
        this.toSequenceSize = toSequenceSize;
        this.toChainStart = toChainStart;
        // not used
        //this.score = score;
        this.fromChainEnd = fromChainEnd;
        this.fromSequenceName = fromSequenceName;
        this.fromSequenceSize = fromSequenceSize;
        this.fromChainStart = fromChainStart;
        this.id = id;
    }


    /**
     * Holds a range that continuously lines up between target and query genome builds.
     * Indices are 0-based, half-open.
     */
    static class ContinuousBlock {
        final int fromStart;	  /* Start of range covered in "from". */
        final int toStart;		  /* Range covered in "to". */
        final int blockLength;    /* length of continuous block of that maps btw from and to */
        //int score;	 	 	  /* Score of block. */

        private ContinuousBlock(final int fromStart, final int toStart, final int blockLength) {
            this.fromStart = fromStart;
            this.toStart = toStart;
            this.blockLength = blockLength;
        }

        /**
         * @return 0-based, half-open end of region in "from"
         */
        int getFromEnd() {
            return fromStart + blockLength;
        }

        /**
         * @return 0-based, half-open end of region in "to"
         */
        int getToEnd() {
            return toStart + blockLength;
        }
    }

    private void addBlock(final int tStart, final int qStart, final int blockLength) {
        blockList.add(new ContinuousBlock(tStart, qStart, blockLength));
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
     * @param chainFile File in UCSC chain format.
     * @return OverlapDetector will all Chains from reader loaded into it.
     */
    static OverlapDetector<Chain> loadChains(final File chainFile) {
        AsciiLineReader reader = new AsciiLineReader(IoUtil.openFileForReading(chainFile));
        final OverlapDetector<Chain> ret = new OverlapDetector<Chain>(0, 0);
        Chain chain;
        while ((chain = Chain.loadChain(reader, chainFile.toString())) != null) {
            ret.addLhs(chain, chain.interval);
        }
        reader.close();
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
        String fromSequenceName = null;
        int fromSequenceSize = 0;
        int fromChainStart = 0;
        int fromChainEnd = 0;
        String toSequenceName = null;
        int toSequenceSize = 0;
        boolean toNegativeStrand = false;
        int toChainStart = 0;
        int toChainEnd = 0;
        int id = 0;
        try {
            score = Double.parseDouble(chainFields[1]);
            fromSequenceName = chainFields[2];
            fromSequenceSize = Integer.parseInt(chainFields[3]);
            // Strand ignored because it is always +
            fromChainStart = Integer.parseInt(chainFields[5]);
            fromChainEnd = Integer.parseInt(chainFields[6]);
            toSequenceName = chainFields[7];
            toSequenceSize = Integer.parseInt(chainFields[8]);
            toNegativeStrand = chainFields[9].equals("-");
            toChainStart = Integer.parseInt(chainFields[10]);
            toChainEnd = Integer.parseInt(chainFields[11]);
            id = Integer.parseInt(chainFields[12]);
        } catch (NumberFormatException e) {
            throwChainFileParseException("Invalid field", chainFile, reader.getLineNumber());
        }
        final Chain chain = new Chain(score, fromSequenceName, fromSequenceSize, fromChainStart, fromChainEnd, toSequenceName, toSequenceSize, toNegativeStrand, toChainStart,
                toChainEnd, id);
        int toBlockStart = chain.toChainStart;
        int fromBlockStart = chain.fromChainStart;
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
            chain.addBlock(fromBlockStart, toBlockStart, size);
            if (!sawLastLine) {
                fromBlockStart += Integer.parseInt(blockFields[1]) + size;
                toBlockStart += Integer.parseInt(blockFields[2]) + size;
            }

        }
        return chain;
    }

    private static void throwChainFileParseException(final String message, final String chainFile, final int lineNumber) {
        throw new PicardException(message + " in chain file " + chainFile + " at line " + lineNumber);
    }
}
