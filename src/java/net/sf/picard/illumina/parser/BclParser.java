/*
 * The MIT License
 *
 * Copyright (c) 2012 The Broad Institute
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
package net.sf.picard.illumina.parser;


import net.sf.picard.illumina.parser.readers.BclReader;

import java.io.File;
import java.util.Collections;
import java.util.NoSuchElementException;
import java.util.Set;

import static net.sf.picard.util.CollectionUtil.makeSet;

/**
 * BclParser parses a number of BclFiles equal to the total of all the values in outputLengths and returns a BclData object
 * segmented based on these lengths.  The only client of this class should be IlluminaDataProvider and an test classes.  See BclReader for
 * more information on BclFiles.  BclParser provides support for reading BaseCalls and QualityScores.
 */
class BclParser extends PerTilePerCycleParser<BclData>{
    private static final int EAMSS_M2_GE_THRESHOLD = 30;
    private static final int EAMSS_S1_LT_THRESHOLD = 15; //was 15
    public static final byte MASKING_QUALITY = (byte) 0x02;

    private static final Set<IlluminaDataType> SUPPORTED_TYPES = Collections.unmodifiableSet(makeSet(IlluminaDataType.BaseCalls, IlluminaDataType.QualityScores));

    public BclParser(final File directory, final int lane, final CycleIlluminaFileMap tilesToCycleFiles, final OutputMapping outputMapping) {
        super(directory, lane, tilesToCycleFiles, outputMapping);
    }

    /** Create the BclData object segmented by the given outputLengths */
    @Override
    protected BclData makeData(final int[] outputLengths) {
        return new BclData(outputLengths);
    }

    /** Create a Bcl parser for an individual cycle and wrap it with the CycleFileParser interaface which populates
     * the correct cycle in BclData.
     * @param file The file to parse
     * @param cycle The cycle that file represents
     * @return A CycleFileParser that populates a BclData object with data for a single cycle
     */
    @Override
    protected CycleFileParser<BclData> makeCycleFileParser(final File file, final int cycle) {
        return new CycleFileParser<BclData>(){
            final OutputMapping.TwoDIndex cycleOutputIndex = outputMapping.getOutputIndexForCycle(cycle);
            BclReader reader = new BclReader(file);

            @Override
            public void close() {
                reader = null;
            }

            @Override
            public void next(final BclData ild) {
                if(!hasNext()) {
                    throw new NoSuchElementException();
                }

                final BclReader.BclValue value = reader.next();
                ild.getBases()    [cycleOutputIndex.majorIndex][cycleOutputIndex.minorIndex] = value.base;
                ild.getQualities()[cycleOutputIndex.majorIndex][cycleOutputIndex.minorIndex] = value.quality;
            }

            @Override
            public boolean hasNext() {
                return reader != null && reader.hasNext();
            }
        };
    }

    @Override
    public Set<IlluminaDataType> supportedTypes() {
        return SUPPORTED_TYPES;
    }


    @Override
    public BclData next() {
        final BclData bclData = super.next();

        final byte [][] bases     = bclData.bases;
        final byte [][] qualities = bclData.qualities;

        //first run EAMSS
        for(int i = 0; i < bases.length; i++) {
            runEamssForReadInPlace(bases[i], qualities[i]);
        }

        return bclData;
    }

    /**
     * EAMSS is an Illumina Developed Algorithm for detecting reads whose quality has deteriorated towards
     * their end and revising the quality to the masking quality (2) if this is the case.  This algorithm
     * works as follows (with one exception):
     *
     * Start at the end (high indices, at the right below) of the read and calculate an EAMSS tally at each
     * location as follow:
     * if(quality[i] < 15) tally += 1
     * if(quality[i] >= 15 and < 30) tally = tally
     * if(quality[i] >= 30) tally -= 2
     *
     *
     * For each location, keep track of this tally (e.g.)
     * Read Starts at <- this end
     * Cycle:       1  2  3  4  5  6  7  8  9
     * Bases:       A  C  T  G  G  G  T  C  A
     * Qualities:   32 32 16 15 8  10 32 2  2
     * Cycle Score: -2 -2 0  0  1  1  -2 1  1           //The EAMSS Score determined for this cycle alone
     * EAMSS TALLY: 0  0  2  2  2  1  0  2  1
     *                    X - Earliest instance of Max-Score
     *
     * You must keep track of the maximum EAMSS tally (in this case 2) and the earliest(lowest) cycle at which
     * it occurs.  If and only if, the max EAMSS tally >= 1 then from there until the end(highest cycle) of the
     * read reassign these qualities as 2 (the masking quality).  The output qualities would therefore be
     * transformed from:
     *
     * Original Qualities: 32 32 16 15 8  10 32 2  2    to
     * Final    Qualities: 32 32 2  2  2  2  2  2  2
     *                           X - Earliest instance of max-tally/end of masking
     *
     * IMPORTANT:
     * The one exception is: If the max EAMSS Tally is preceded by a long string of G basecalls (10 or more, with a single basecall exception per10 bases)
     * then the masking continues to the beginning of that string of G's. E.g.:
     *
     * Cycle:       1  2  3  4  5  6  7  8   9  10 11 12 13 14 15 16 17 18
     * Bases:       C  T  A  C  A  G  A  G   G  G  G  G  G  G  G  C  A  T
     * Qualities:   30 22 26 27 28 30 7  34  20 19 38 15 32 32 10 4  2  5
     * Cycle Score: -2  0  0  0  0 -2 1  -2  0  0  -2 0  -2 -2  1 1  1  1
     * EAMSS TALLY: -2 -5 -5 -5 -5 -5 -3 -4 -2 -2  -2 0   0  2  4 3  2  1
     *                                                          X- Earliest instance of Max-Tally
     *
     * Resulting Transformation:
     * Bases:                C  T  A  C  A  G  A   G   G  G  G  G  G  G  G  C  A  T
     * Original Qualities:   30 22 26 27 28 30 7  34  20 19 38 15 32 32 10  4  2  5
     * Final    Qualities:   30 22 26 27 28  2 2   2   2  2  2  2  2  2  2  2  2  2
     *                                                          X- Earliest instance of Max-Tally
     *                                       X - Start of EAMSS masking due to G-Run
     *
     * To further clarify the exception rule here are a few examples:
     * A C G A C G G G G G G G G G G G G G G G G G G G G A C T
     *                                                 X - Earliest instance of Max-Tally
     *     X - Start of EAMSS masking (with a two base call jump because we have 20 bases in the run already)
     *
     * T T G G A G G G G G G G G G G G G G G G G G G A G A C T
     *                                                 X - Earliest instance of Max-Tally
     *         X - We can skip this A as well as the earlier A because we have 20 or more bases in the run already
     *     X - Start of EAMSS masking (with a two base call jump because we have 20 bases in the run)
     *
     * T T G G G A A G G G G G G G G G G G G G G G G G G T T A T
     *                                                 X - Earliest instance of Max-Tally
     *           X X - WE can skip these bases because the first A counts as the first skip and as far as the length of the string of G's is
     *                 concerned, these are both counted like G's
     *           X - This A is the 20th base in the string of G's and therefore can be skipped
     *         X - Note that the A's previous to the G's are only included because there are G's further on that are within the number
     *             of allowable exceptions away (i.e. 2 in this instance), if there were NO G's after the A's you CANNOT count the A's
     *             as part of the G strings (even if no exceptions have previously occured) In other words, the end of the string of G's
     *             MUST end in a G not an "exception"
     *
     * However, if the max-tally occurs to the right of the run of Gs then this is still part of the string of G's but does count towards
     * the number of exceptions allowable
     * (e.g.)
     * T T G G G G G G G G G G A C G
     *                         X - Earliest instance of Max-tally
     *  The first index CAN be considered as an exception, the above would be masked to
     *  the following point:
     *  T T G G G G G G G G G G A C G
     *      X - End of EAMSS masking due to G-Run
     *
     * To sum up the final points, a string of G's CAN START with an exception but CANNOT END in an exception.
     *
     *
     * @param bases Bases for a single read in the cluster ( not the entire cluster )
     * @param qualities Qualities for a single read in the cluster ( not the entire cluster )
     */
    protected static void runEamssForReadInPlace(final byte [] bases, final byte [] qualities) {
        int eamssTally = 0;
        int maxTally = Integer.MIN_VALUE;
        int indexOfMax = -1;

        for(int i = bases.length - 1; i >= 0; i--) {
            final int quality = (0xff & qualities[i]);

            if(quality >= EAMSS_M2_GE_THRESHOLD) {
                eamssTally -= 2;
            } else if(quality < EAMSS_S1_LT_THRESHOLD) {
                eamssTally += 1;
            }

            if(eamssTally >= maxTally) {
                indexOfMax = i;
                maxTally = eamssTally;
            }
        }

        if(maxTally >= 1) {
            int numGs = 0;
            int exceptions = 0;

            for(int i = indexOfMax; i >= 0; i--) {
                if(bases[i] == 'G') {
                    ++numGs;
                } else {
                    Integer skip = skipBy(i, numGs, exceptions, bases);
                    if(skip != null) {
                        exceptions += skip;
                        numGs      += skip;
                        i          -= (skip-1);
                    } else {
                        break;
                    }
                }
            }

            if(numGs >= 10) {
                indexOfMax = (indexOfMax+1) - numGs;
            }
            
            for(int i = indexOfMax; i < qualities.length; i++) {
                qualities[i] = MASKING_QUALITY;
            }
        }
    }

    /** Determine whether or not the base at index is part of a skippable section in a run of G's, if so
     * return the number of bases that the section is composed of.
     * @param index Current index, which should be the index of a non-'G' base
     * @param numGs The number of bases in the current string of G's for this read
     * @param prevExceptions The number of exceptions previously detected in this string by this method
     * @param bases The bases of this read
     * @return If we have not reached our exception limit (1/every 10bases) and a G is within exceptionLimit(numGs/10)
     * indices before the current index then return index - (index of next g), else return null  Null indicates this is
     * NOT a skippable region, if we run into index 0 without finding a g then NULL is also returned
     */
    private static Integer skipBy(int index, int numGs, final int prevExceptions, final byte [] bases) {
        Integer skip = null;
        for(int backup = 1; backup <= index; backup++) {
            final int exceptionLimit = Math.max((numGs + backup) / 10, 1);
            if(prevExceptions + backup > exceptionLimit) {
                break;
            }
            if(bases[index - backup] == 'G') {
                skip = backup;
                break;
            }
        }

        return skip;
    }
}

/** A class that implements the IlluminaData interfaces provided by this parser
 * One BclData object is returned to IlluminaDataProvider per cluster and each
 * first level array in bases and qualities represents a single read in that
 * cluster */
class BclData implements BaseData, QualityData {
    public final byte [][] bases;
    public final byte [][] qualities;

    public BclData(final int[] outputLengths) {
        bases     = new byte[outputLengths.length][];
        qualities = new byte[outputLengths.length][];

        for(int i = 0; i < outputLengths.length; i++) {
            bases[i]     = new byte[outputLengths[i]];
            qualities[i] = new byte[outputLengths[i]];
        }
    }

    @Override
    public byte[][] getBases() {
        return bases;
    }

    @Override
    public byte[][] getQualities() {
        return qualities;
    }
}