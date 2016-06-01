/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

package picard.sam.markduplicates;

import htsjdk.samtools.DuplicateScoringStrategy;
import picard.cmdline.CommandLineProgram;

/**
 * @author nhomer
 */
public class SimpleMarkDuplicatesWithMateCigarTester extends AbstractMarkDuplicatesCommandLineProgramTester {

    public SimpleMarkDuplicatesWithMateCigarTester() {
        // NB: to be equivalent to MarkDuplicates we need to use SUM_OF_BASE_QUALITIES
        super(DuplicateScoringStrategy.ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH);

        addArg("MAX_RECORDS_IN_RAM=1000");
    }

    @Override
    protected CommandLineProgram getProgram() { return new SimpleMarkDuplicatesWithMateCigar(); }
}

