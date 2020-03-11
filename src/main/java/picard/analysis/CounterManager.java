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

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Class for managing a list of Counters of integer,
 * provides methods to access data from Counters with respect to an offset.
 * Each Counter represents a certain region of a processed sequence starting from offset and
 * can accumulate some value for each locus of this sequence, for instance, read coverage.
 * @author Mariia_Zueva@epam.com, EPAM Systems, Inc. <www.epam.com>
 */

public class CounterManager {

    /**
     * Length of inner Counter arrays
     */
    private final int arrayLength;
    /**
     * Proposed length of processed reads
     */
    private final int readLength;
    /**
     * Current offset from the start of the processed sequence. If offset = 100, this means, that value in
     * Counter arrays at index 0, represents the 100th position of the sequence.
     */
    private int offset = 0;
    /**
     * List of Counter arrays for managing
     */
    private final ArrayList<Counter> arrays = new ArrayList<>(2);

    /**
     * Constructor creates new CounterManager without any Counters,
     * counters are added to CounterManager via newCounter() method
     *
     * @param arrayLength length of inner Counter arrays
     * @param readLength  proposed length of processed reads
     */
    public CounterManager(final int arrayLength, int readLength) {
        this.arrayLength = arrayLength;
        this.readLength = readLength;
    }

    /**
     * Method added for use in unit tests
     *
     * @return offset from the reference sequence, that is represented by index 0 of {@link Counter} arrays.
     */
    int getOffset() {
        return offset;
    }

    /**
     * Method added for use in unit tests
     *
     * @param offset from the reference sequence, that will be represented by index 0 of {@link Counter} arrays.
     */
    void setOffset(int offset) {
        this.offset = offset;
    }

    /**
     * Method checks that new locus position is not out of bounds of Counter arrays and there is enough
     * space in them to hold information on at least one more read of length {@link #readLength}.
     * If there is no free space, but there is accumulated information in Counter arrays after new locus position
     * that we may need, the arrays are rebased, so that 0 index of arrays represents new locus position.
     * In other case, we just clear the arrays.
     *
     * @param locusPosition position in the reference sequence,
     *                      that will be represented by index 0 of {@link Counter} arrays.
     */
    public void checkOutOfBounds(int locusPosition) {
        if (locusPosition - offset + readLength >= arrayLength) {
            if (locusPosition - offset < arrayLength) {
                rebase(locusPosition);
            } else {
                clear();
                offset = locusPosition;
            }
        }
    }

    /**
     * Rebases inner Counter arrays so that 0 index of array represents the new locus position
     *
     * @param locusPosition position in the reference sequence,
     *                      that will be represented by index 0 of {@link Counter} arrays.
     */
    private void rebase(int locusPosition) {
        if (locusPosition < offset) {
            throw new IllegalArgumentException("Position in the reference sequence is lesser than offset.");
        }
        for (Counter counter : arrays) {
            final int[] array = counter.array;
            final int skipLength = locusPosition - offset;
            System.arraycopy(array, skipLength, array, 0, arrayLength - skipLength);
            Arrays.fill(array, arrayLength - skipLength, arrayLength, 0);
        }
        offset = locusPosition;
    }

    /**
     * Clears all inner Counter arrays
     */
    public void clear() {
        for (Counter counter : arrays) {
            final int[] array = counter.array;
            Arrays.fill(array, 0);
        }
        offset = 0;
    }

    /**
     * Creates a new Counter object and adds it to the list of managed Counters.
     *
     * @return {@link Counter}, that will be managed by current {@link CounterManager}
     */
    public Counter newCounter() {
        final Counter counter = new Counter(arrayLength);
        arrays.add(counter);
        return counter;
    }

    /**
     * Class represents an integer array with methods to increment and get the values from it with respect
     * to offset of outer {@link CounterManager}.
     */
    public class Counter {

        /**
         * Inner array to store accumulated values
         */
        private final int[] array;

        /**
         * @param arrayLength length of the inner array
         */
        private Counter(int arrayLength) {
            array = new int[arrayLength];
        }

        /**
         * Increments value corresponding to reference sequence index
         *
         * @param index in the reference sequence
         */
        public void increment(int index) {
            checkIndex(index);
            array[index - offset]++;
        }

        /**
         * @param index in the reference sequence
         * @return value from the inner array, corresponding to input reference index
         */
        public int get(int index) {
            checkIndex(index);
            return array[index - offset];
        }

        private void checkIndex(int index) {
            if ((index - offset) < 0 || (index - offset) >= array.length) {
                throw new ArrayIndexOutOfBoundsException("The requested index " + index + " is out of counter bounds. " +
                        "Possible cause of exception can be wrong READ_LENGTH parameter (much smaller than actual read length)");
            }
        }
    }
}
