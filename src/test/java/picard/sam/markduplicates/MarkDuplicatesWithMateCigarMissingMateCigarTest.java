/*
 * The MIT License
 *
 * Copyright (c) 2026 The Broad Institute
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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMTag;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.sam.DuplicationMetrics;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import static picard.sam.markduplicates.MateCigarAndUmiTestSupport.*;

/**
 * Tests {@link MarkDuplicatesWithMateCigar} on a pair of which only one end lacks the MC tag. Pairs pX1, pX2 and pX9
 * have their first ends forward at chr1:1000 and their second ends reverse at chr1:1200 (50M), with equal scores, so
 * pX1, the smallest name, is kept; only pX9 lacks an MC tag, on its first or on its second end.
 *
 * With SKIP_PAIRS_WITH_NO_MATE_CIGAR (the default) the end without MC is written unchanged and unmarked, and pX1 and
 * pX2 are marked as without pX9. Without it the program fails, naming the read.
 */
public class MarkDuplicatesWithMateCigarMissingMateCigarTest {

    private static File input(final boolean firstEndLacksMateCigar) throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        addPair(builder, "pX1", 0, 1000, 1200);
        addPair(builder, "pX2", 0, 1000, 1200);
        final List<SAMRecord> half = addPair(builder, "pX9", 0, 1000, 1200);
        half.get(firstEndLacksMateCigar ? 0 : 1).setAttribute(SAMTag.MC.name(), null);
        return writeCoordinateSorted(builder);
    }

    private static File[] run(final File input, final String... extraArgs) throws IOException {
        final File output = newTempFile("MarkDuplicatesWithMateCigarMissingMateCigarTest", ".sam");
        final File metrics = newTempFile("MarkDuplicatesWithMateCigarMissingMateCigarTest", ".metrics");
        final String[] args = new String[3 + extraArgs.length];
        args[0] = "INPUT=" + input.getAbsolutePath();
        args[1] = "OUTPUT=" + output.getAbsolutePath();
        args[2] = "METRICS_FILE=" + metrics.getAbsolutePath();
        System.arraycopy(extraArgs, 0, args, 3, extraArgs.length);
        Assert.assertEquals(new MarkDuplicatesWithMateCigar().instanceMain(args), 0);
        return new File[]{output, metrics};
    }

    @DataProvider(name = "endWithoutMateCigar")
    public Object[][] endWithoutMateCigar() {
        return new Object[][]{{true}, {false}};
    }

    @Test(dataProvider = "endWithoutMateCigar")
    public void testEndWithoutMateCigarIsSkipped(final boolean firstEndLacksMateCigar) throws IOException {
        final File input = input(firstEndLacksMateCigar);
        final List<SAMRecord> inputRecords = readRecords(input);
        final String skipped = firstEndLacksMateCigar ? "pX9/1" : "pX9/2";

        final File[] outputs = run(input);
        final List<SAMRecord> records = readRecords(outputs[0]);

        Assert.assertEquals(keys(records), keys(inputRecords));
        final Map<String, SAMRecord> byKey = byKey(records);
        for (final String key : List.of("pX1/1", "pX1/2", "pX2/1", "pX2/2", skipped)) {
            Assert.assertEquals(byKey.get(key).getDuplicateReadFlag(), key.startsWith("pX2"), key);
        }
        final SAMRecord skippedOutput = byKey.get(skipped).deepCopy();
        skippedOutput.setAttribute(SAMTag.PG.name(), null);
        Assert.assertEquals(skippedOutput.getSAMString(), byKey(inputRecords).get(skipped).getSAMString());
        Assert.assertNull(skippedOutput.getAttribute(SAMTag.MC.name()));

        final DuplicationMetrics metrics = readDuplicationMetrics(outputs[1]).get("libA");
        Assert.assertEquals(metrics.UNPAIRED_READS_EXAMINED, 0);
        Assert.assertEquals(metrics.READ_PAIRS_EXAMINED, 3);
        Assert.assertEquals(metrics.READ_PAIR_OPTICAL_DUPLICATES, 0);
    }

    @Test(dataProvider = "endWithoutMateCigar")
    public void testEndWithoutMateCigarFailsWhenNotSkipped(final boolean firstEndLacksMateCigar) throws IOException {
        final File input = input(firstEndLacksMateCigar);
        final PicardException e = Assert.expectThrows(PicardException.class,
                () -> run(input, "SKIP_PAIRS_WITH_NO_MATE_CIGAR=false"));
        Assert.assertTrue(e.getMessage().contains("Read pX9 was mapped and had a mapped mate, but no mate cigar"), e.getMessage());
    }
}
