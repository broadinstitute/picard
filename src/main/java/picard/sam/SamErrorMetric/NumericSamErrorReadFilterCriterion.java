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

public class NumericSamErrorReadFilterCriterion extends SamErrorReadFilterCriterion<Integer> {
    public NumericSamErrorReadFilterCriterion(final SamErrorReadFilter.Comparator comparator, final Integer value) {
        super(comparator, value);
    }

    @Override
    public void checkCriterion(final Integer stratus) {
        switch (comparator) {
            case Equal: if (stratus.equals(value)) { setSatisifed(); } break;
            case NotEqual: if (!stratus.equals(value)) { setSatisifed(); } break;
            case Smaller: if (stratus < value) { setSatisifed(); } break;
            case SmallerOrEqual: if (stratus <= value) { setSatisifed(); } break;
            case Greater: if (stratus > value) { setSatisifed(); } break;
            case GreaterOrEqual: if (stratus >= value) { setSatisifed(); } break;
            default: throw new IllegalArgumentException("Invalid operator for this criterion type.");
        }
    }
}