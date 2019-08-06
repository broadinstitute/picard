/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
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


package picard.sam.SamErrorMetric;

/**
 * Class holding a criterion to be used in the SamErrorReadFilter class
 *
 * @param <T>
 */
public abstract class SamErrorReadFilterCriterion<T extends Comparable<T>> {
    public final SamErrorReadFilter.Comparator comparator;
    public final T value;
    private boolean satisifed = false;

    /**
     * Creates a criterion object
     *
     * @param comparator Comparison operation to perform
     * @param value      Value to compare with to satisfy the criterion
     */
    public SamErrorReadFilterCriterion(final SamErrorReadFilter.Comparator comparator, final T value) {
        this.comparator = comparator;
        this.value = value;
    }

    /**
     * Method to be called in order to check whether this criterion is satisfied
     *
     * @param stratus Object to check whether it satisfies the criterion
     */
    public abstract void checkCriterion(final T stratus);

    /**
     * Returns whether or not the criterion is satisfied
     */
    public boolean isSatisifed() {
        return satisifed;
    }

    protected void setSatisifed() {
        satisifed = true;
    }

    /**
     * Resets this criterion. This needs to be called whenever a new RecordAndOffset is considered for filtering.
     */
    public void reset() {
        satisifed = false;
    }
}
