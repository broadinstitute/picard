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
package net.sf.samtools;

import org.testng.annotations.Test;
import org.testng.Assert;

import java.util.Arrays;

import net.sf.samtools.Cigar;
import net.sf.samtools.TextCigarCodec;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.BinaryCigarCodec;

public class CigarCodecTest {

    private final BinaryCigarCodec binaryCigarCodec = new BinaryCigarCodec();
    private final TextCigarCodec textCigarCodec = new TextCigarCodec();


    @Test
    public void testDefault() {
        final Cigar emptyCigar = new Cigar();
        Assert.assertEquals(emptyCigar, binaryCigarCodec.decode(new int[0]));
        final int[] binaryCigar = binaryCigarCodec.encode(emptyCigar);
        Assert.assertEquals(0, binaryCigar.length);
        Assert.assertEquals(emptyCigar, textCigarCodec.decode(SAMRecord.NO_ALIGNMENT_CIGAR));
        Assert.assertEquals(textCigarCodec.encode(emptyCigar), SAMRecord.NO_ALIGNMENT_CIGAR);
    }

    private static class Cigarette {
        final int length;
        final char op;

        private Cigarette(final int length, final char op) {
            this.length = length;
            this.op = op;
        }

        int getBinaryOp() {
            switch (op) {
                case 'M': return 0;
                case 'I': return 1;
                case 'D': return 2;
                case 'N': return 3;
                case 'S': return 4;
                case 'H': return 5;
                case 'P': return 6;
                case 'C': return 7;
                default: Assert.assertTrue(false);
            }
            return -1;
        }
    }

    private String makeTextCigar(final Cigarette[] cigarettes) {
        final StringBuilder sb = new StringBuilder();
        for (final Cigarette c : cigarettes) {
            sb.append(Integer.toString(c.length));
            sb.append(c.op);
        }
        return sb.toString();
    }

    private int[] makeBinaryCigar(final Cigarette[] cigarettes) {
        final int[] ret = new int[cigarettes.length];
        for (int i = 0; i < cigarettes.length; ++i) {
            ret[i] = cigarettes[i].length << 4 | cigarettes[i].getBinaryOp();
        }
        return ret;
    }

    @Test
    public void testSimple() {
        final Cigarette[] cigarettes = {
                new Cigarette(100, 'M'),
                new Cigarette(200, 'I'),
                new Cigarette(50, 'D'),
                new Cigarette(21, 'N'),
                new Cigarette(12, 'S'),
                new Cigarette(99, 'H'),
                new Cigarette(20, 'P'),
                new Cigarette(2, 'C'),
        };
        final String textCigar = makeTextCigar(cigarettes);
        final int[] binaryCigar = makeBinaryCigar(cigarettes);
        final Cigar fromText = textCigarCodec.decode(textCigar);
        final Cigar fromBinary = binaryCigarCodec.decode(binaryCigar);
        Assert.assertEquals(fromText, fromBinary);
        final String anotherTextCigar = textCigarCodec.encode(fromBinary);
        final int[] anotherBinaryCigar = binaryCigarCodec.encode(fromText);
        Assert.assertEquals(anotherTextCigar, textCigar);
        Assert.assertTrue(Arrays.equals(anotherBinaryCigar, binaryCigar));
    }
}
