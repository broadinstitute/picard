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
import java.io.PrintWriter;
import java.util.*;
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

    /** Score is not used in basic liftover implementation, but is stored so that chain can be written to disk. */
    final double score;
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
        this.score = score;
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

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final ContinuousBlock that = (ContinuousBlock) o;

            if (blockLength != that.blockLength) return false;
            if (fromStart != that.fromStart) return false;
            if (toStart != that.toStart) return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result = fromStart;
            result = 31 * result + toStart;
            result = 31 * result + blockLength;
            return result;
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

    void write(final PrintWriter writer) {
        writer.printf("chain\t%f\t%s\t%d\t+\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%d\n",
                score, fromSequenceName, fromSequenceSize, fromChainStart, fromChainEnd,
                toSequenceName, toSequenceSize, (toNegativeStrand? "-": "+"), toChainStart, toChainEnd, id);
        for (int i = 0; i < blockList.size() - 1; ++i) {
            final ContinuousBlock thisBlock = blockList.get(i);
            final ContinuousBlock nextBlock = blockList.get(i+1);

            final int fromGap = nextBlock.fromStart - thisBlock.getFromEnd();
            final int toGap = nextBlock.toStart - thisBlock.getToEnd();
            writer.printf("%d\t%d\t%d\n", thisBlock.blockLength, fromGap, toGap);
        }
        writer.printf("%d\n", blockList.get(blockList.size() - 1).blockLength);
        writer.println();
    }

    /**
     * Throw an exception if Chain looks strange.
     */
    void validate() {
        validatePositive("fromSequenceSize", fromSequenceSize);
        validateNonNegative("fromChainStart", fromChainStart);
        validateNonNegative("fromChainEnd", fromChainEnd);
        validatePositive("toSequenceSize", toSequenceSize);
        validateNonNegative("toChainStart", toChainStart);
        validateNonNegative("toChainEnd", toChainEnd);
        int fromLength = fromChainEnd - fromChainStart;
        validatePositive("from length", fromLength);
        int toLength = toChainEnd - toChainStart;
        validatePositive("to length", toLength);
        if (fromLength > fromSequenceSize) throw new PicardException("From chain length (" + fromLength +
                ") < from sequence length (" + fromSequenceSize + ") for chain " + id);
        if (toLength > toSequenceSize) throw new PicardException("To chain length (" + toLength +
                ") < to sequence length (" + toSequenceSize + ") for chain " + id);
        if (fromSequenceName.isEmpty()) throw new PicardException("Chain " + id + "has empty from sequence name.");
        if (toSequenceName.isEmpty()) throw new PicardException("Chain " + id + "has empty to sequence name.");
        if (blockList.isEmpty()) throw new PicardException("Chain " + id + " has empty block list.");
        final ContinuousBlock firstBlock = blockList.get(0);
        if (firstBlock.fromStart != fromChainStart) {
            throw new PicardException("First block from start != chain from start for chain " + id);
        }
        if (firstBlock.toStart != toChainStart) {
            throw new PicardException("First block to start != chain to start for chain " + id);
        }
        final ContinuousBlock lastBlock = blockList.get(blockList.size() - 1);
        if (lastBlock.getFromEnd() != fromChainEnd) {
            throw new PicardException("Last block from end != chain from end for chain " + id);
        }
        if (lastBlock.getToEnd() != toChainEnd) {
            throw new PicardException("Last block to end < chain to end for chain " + id);
        }
        for (int i = 1; i < blockList.size(); ++i) {
            final ContinuousBlock thisBlock = blockList.get(i);
            final ContinuousBlock prevBlock = blockList.get(i-1);
            if (thisBlock.fromStart < prevBlock.getFromEnd()) {
                throw new PicardException("Continuous block " + i + " from starts before previous block ends for chain " + id);
            }
            if (thisBlock.toStart < prevBlock.getToEnd()) {
                throw new PicardException("Continuous block " + i + " to starts before previous block ends for chain " + id);
            }
        }
    }

    private void validatePositive(final String attributeName, final int attribute) {
        if (attribute <= 0) {
            throw new PicardException(attributeName + " is not positive: " + attribute + " for chain " + id);
        }
    }

    private void validateNonNegative(final String attributeName, final int attribute) {
        if (attribute < 0) {
            throw new PicardException(attributeName + " is negative: " + attribute + " for chain " + id);
        }
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final Chain chain = (Chain) o;

        if (fromChainEnd != chain.fromChainEnd) return false;
        if (fromChainStart != chain.fromChainStart) return false;
        if (fromSequenceSize != chain.fromSequenceSize) return false;
        if (id != chain.id) return false;
        if (Double.compare(chain.score, score) != 0) return false;
        if (toChainEnd != chain.toChainEnd) return false;
        if (toChainStart != chain.toChainStart) return false;
        if (toNegativeStrand != chain.toNegativeStrand) return false;
        if (toSequenceSize != chain.toSequenceSize) return false;
        if (blockList != null ? !blockList.equals(chain.blockList) : chain.blockList != null) return false;
        if (fromSequenceName != null ? !fromSequenceName.equals(chain.fromSequenceName) : chain.fromSequenceName != null)
            return false;
        if (interval != null ? !interval.equals(chain.interval) : chain.interval != null) return false;
        if (toSequenceName != null ? !toSequenceName.equals(chain.toSequenceName) : chain.toSequenceName != null)
            return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        temp = score != +0.0d ? Double.doubleToLongBits(score) : 0L;
        result = (int) (temp ^ (temp >>> 32));
        result = 31 * result + (interval != null ? interval.hashCode() : 0);
        result = 31 * result + (fromSequenceName != null ? fromSequenceName.hashCode() : 0);
        result = 31 * result + fromSequenceSize;
        result = 31 * result + fromChainStart;
        result = 31 * result + fromChainEnd;
        result = 31 * result + (toSequenceName != null ? toSequenceName.hashCode() : 0);
        result = 31 * result + toSequenceSize;
        result = 31 * result + (toNegativeStrand ? 1 : 0);
        result = 31 * result + toChainStart;
        result = 31 * result + toChainEnd;
        result = 31 * result + id;
        result = 31 * result + (blockList != null ? blockList.hashCode() : 0);
        return result;
    }

    /**
     * Read all the chains and load into an OverlapDetector.
     * @param chainFile File in UCSC chain format.
     * @return OverlapDetector will all Chains from reader loaded into it.
     */
    static OverlapDetector<Chain> loadChains(final File chainFile) {
        final Set<Integer> ids = new HashSet<Integer>();
        AsciiLineReader reader = new AsciiLineReader(IoUtil.openFileForReading(chainFile));
        final OverlapDetector<Chain> ret = new OverlapDetector<Chain>(0, 0);
        Chain chain;
        while ((chain = Chain.loadChain(reader, chainFile.toString())) != null) {
            if (ids.contains(chain.id)) {
                throw new PicardException("Chain id " + chain.id + " appears more than once in chain file.");
            }
            ids.add(chain.id);
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
        chain.validate();
        return chain;
    }

    private static void throwChainFileParseException(final String message, final String chainFile, final int lineNumber) {
        throw new PicardException(message + " in chain file " + chainFile + " at line " + lineNumber);
    }
}
