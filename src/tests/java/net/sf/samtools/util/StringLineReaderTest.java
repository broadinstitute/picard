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
package net.sf.samtools.util;

import org.testng.Assert;
import org.testng.annotations.Test;

public class StringLineReaderTest {

    private static final String[] TERMINATORS = {"\r", "\n", "\r\n"};
    private static final boolean[] TERMINATOR_INCLUSION = {false, true};
    private static final boolean[] LAST_LINE_TERMINATED = {false, true};

    enum EmptyLineState {
        FIRST_LINE, LAST_LINE, MIDDLE_LINE, COMPLETELY_EMPTY
    }

    /**
     * Test a bunch of combinations instead of writing a method for each.
     */
    @Test
    public void testBasic() {
        for (final String terminator : TERMINATORS) {
            for (final boolean includeTerminators : TERMINATOR_INCLUSION) {
                for (final boolean lastLineTerminated : LAST_LINE_TERMINATED) {
                    for (final EmptyLineState emptyLineState : EmptyLineState.values()) {
                        if (emptyLineState == EmptyLineState.COMPLETELY_EMPTY) {
                            emptyTestHelper(terminator, includeTerminators, lastLineTerminated);
                        } else {
                            testHelper(terminator, includeTerminators, lastLineTerminated, emptyLineState);
                        }
                    }
                }
            }
        }
    }

    /**
     * various test cases where there is no input, except perhaps a line terminator
     * @param terminator what the terminator should be in the input
     * @param includeTerminators should the test ask for terminators to be included
     * @param lastLineTerminated does the input have a terminator
     */
    private void emptyTestHelper(final String terminator, final boolean includeTerminators, final boolean lastLineTerminated) {
        final String input;
        if (lastLineTerminated) {
            input = terminator;
        } else {
            input = "";
        }
        final StringLineReader slr = new StringLineReader(input);
        final String output = slr.readLine(includeTerminators);
        if (lastLineTerminated && includeTerminators) {
            Assert.assertEquals(output, "\n");
        } else if (lastLineTerminated) {
            Assert.assertEquals(output, "");
        }
        Assert.assertNull(slr.readLine(includeTerminators));
    }

    /**
     * Test a variety of test cases in which there is more than one line.
     * @param terminator to use in the input
     * @param includeTerminators passed to readLine()
     * @param lastLineTerminated should the input end with a terminator
     * @param emptyLineState where in the input should an empty line be.
     */
    private void testHelper(final String terminator, final boolean includeTerminators, final boolean lastLineTerminated, final EmptyLineState emptyLineState) {
        final String[] lines = new String[3];
        if (emptyLineState == EmptyLineState.FIRST_LINE) {
            lines[0] = "";
            lines[1] = "Hi, Mom!";
            lines[2] = "Hi, Dad?";
        } else if (emptyLineState == EmptyLineState.LAST_LINE) {
            lines[0] = "Hi, Dad?";
            lines[1] = "Hi, Mom!";
            lines[2] = "";
        } else  if (emptyLineState == EmptyLineState.MIDDLE_LINE) {
            lines[0] = "Hi, Dad?";
            lines[1] = "";
            lines[2] = "Hi, Mom!";
        }
        String input = StringUtil.join(terminator, lines);
        if (lastLineTerminated) {
            input = input.concat(terminator);
        }
        final StringLineReader slr = new StringLineReader(input);
        for (int i = 0; i < lines.length - 1; ++i) {
            final String s = slr.readLine(includeTerminators);
            String expected = lines[i];
            if (includeTerminators) {
                expected = expected.concat("\n");
            }
            Assert.assertEquals(s, expected);
        }

        // Last line may need to be handled specially
        String s = slr.readLine(includeTerminators);
        if (!lastLineTerminated && emptyLineState == EmptyLineState.LAST_LINE) {
            Assert.assertNull(s);
        } else {
            String expected = lines[lines.length - 1];
            if (includeTerminators && lastLineTerminated) {
                expected = expected.concat("\n");
            }
            Assert.assertEquals(s, expected);
        }
        s = slr.readLine(includeTerminators);
        Assert.assertNull(s);
    }
}
