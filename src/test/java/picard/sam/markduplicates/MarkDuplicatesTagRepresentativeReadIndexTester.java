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

package picard.sam.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.sam.DuplicationMetrics;
import picard.vcf.MakeVcfSampleNameMap;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * This class is an extension of AbstractMarkDuplicatesCommandLineProgramTester used to test MarkDuplicatesWithMateCigar with SAM files generated on the fly.
 * This performs the underlying tests defined by classes such as see AbstractMarkDuplicatesCommandLineProgramTest and MarkDuplicatesWithMateCigarTest.
 * @author hogstrom@broadinstitute.org
 */
public class MarkDuplicatesTagRepresentativeReadIndexTester extends AbstractMarkDuplicatesCommandLineProgramTester {

    final public Map<String, String> expectedRepresentativeReadMap = new HashMap<>();
    final Map<String, Integer> readIndexMap = new HashMap<>();
    final public Map<String, Integer> expectedSetSizeMap = new HashMap<>();
    public boolean testRepresentativeReads = false;

    public MarkDuplicatesTagRepresentativeReadIndexTester() {
        this(SAMFileHeader.SortOrder.coordinate);
    }

    public MarkDuplicatesTagRepresentativeReadIndexTester(final SAMFileHeader.SortOrder sortOrder) {
        super(SAMRecordSetBuilder.DEFAULT_DUPLICATE_SCORING_STRATEGY, sortOrder);

        addArg("TAGGING_POLICY=All");
        addArg("TAG_DUPLICATE_SET_MEMBERS=true");
    }

    @Override
    public void recordOpticalDuplicatesMarked() {}

    @Override
    protected CommandLineProgram getProgram() { return new MarkDuplicates(); }

    @Override
    void updateExpectationsHook() {
        try (final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(fastaFiles.get(samRecordSetBuilder.getHeader())).open(getOutput())) {
            int indexInFile = 0;
            for (final SAMRecord record : reader) {
                // putIfAbsent because we want to use the index of the first read with this readname
                readIndexMap.putIfAbsent(record.getReadName(), indexInFile);
                indexInFile++;
            }
        } catch (final IOException ex) {
            Assert.fail("Problem updating read file indices", ex);
        }
    }

    @Override
    void testRecordHook(final SAMRecord record){
        if (expectedRepresentativeReadMap.containsKey(record.getReadName())) {
            Assert.assertEquals(record.getAttribute("DI"), readIndexMap.get(expectedRepresentativeReadMap.get(record.getReadName())));
            Assert.assertEquals(record.getAttribute("DS"), expectedSetSizeMap.get(record.getReadName()));
        }
    }

}
