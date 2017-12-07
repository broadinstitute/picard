/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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

import htsjdk.samtools.BAMIndexer;
import htsjdk.samtools.Defaults;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.RuntimeIOException;

import java.io.File;
import java.io.IOException;

public class ParallelMarkDuplicatesTest extends MarkDuplicatesTest {

    @Override
    protected AbstractMarkDuplicatesCommandLineProgramTester getTester() {
        return new ParallelMarkDuplicatesTester();
    }

    @Override
    protected MarkDuplicates prepareMarkDuplicatesWithDuplicateSetTagging(File sam, File outputDir, File outputSam,
                                                                          File metricsFile) {
        final File inputBam = prepareBam(sam, outputDir);

        final ParallelMarkDuplicates markDuplicates = new ParallelMarkDuplicates() {
            @Override
            public String[] getOriginalArgs() {
                return new String[] {
                        "INPUT=" + inputBam.getAbsolutePath(),
                        "OUTPUT=" + outputSam.getAbsolutePath(),
                        "METRICS_FILE=" + metricsFile.getAbsolutePath(),
                        "TMP_DIR=" + outputDir.getAbsolutePath(),
                        "TAG_DUPLICATE_SET_MEMBERS=true"
                };
            }
        };

        markDuplicates.setupOpticalDuplicateFinder();
        markDuplicates.INPUT =  CollectionUtil.makeList(inputBam.getAbsolutePath());
        markDuplicates.OUTPUT = outputSam;
        markDuplicates.METRICS_FILE = metricsFile;
        markDuplicates.TAG_DUPLICATE_SET_MEMBERS = true;
        markDuplicates.PROGRAM_RECORD_ID = null;
        markDuplicates.TMP_DIR = CollectionUtil.makeList(outputDir);

        return markDuplicates;
    }

    @Override
    protected MarkDuplicates prepareMarkDuplicatesWithTaggingPolicy(File sam, File outputDir, File outputSam,
                                                                    File metricsFile, String taggingPolicy) {
        final File inputBam = prepareBam(sam, outputDir);

        final ParallelMarkDuplicates markDuplicates = new ParallelMarkDuplicates() {
            @Override
            public String[] getOriginalArgs() {
                return new String[] {
                        "INPUT=" + inputBam.getAbsolutePath(),
                        "OUTPUT=" + outputSam.getAbsolutePath(),
                        "METRICS_FILE=" + metricsFile.getAbsolutePath(),
                        "TMP_DIR=" + outputDir.getAbsolutePath(),
                        "TAGGING_POLICY=" + taggingPolicy
                };
            }
        };

        markDuplicates.setupOpticalDuplicateFinder();
        markDuplicates.INPUT =  CollectionUtil.makeList(inputBam.getAbsolutePath());
        markDuplicates.OUTPUT = outputSam;
        markDuplicates.METRICS_FILE = metricsFile;
        markDuplicates.TAGGING_POLICY = MarkDuplicates.DuplicateTaggingPolicy.valueOf(taggingPolicy);
        markDuplicates.PROGRAM_RECORD_ID = null;
        markDuplicates.TMP_DIR = CollectionUtil.makeList(outputDir);

        return markDuplicates;
    }

    private File prepareBam(File sam, File outputDir) {
        final File inputBam = new File(outputDir, sam.getName() + ".bam");
        try (
                SamReader samReader = SamReaderFactory.makeDefault().open(sam);
                SAMFileWriter bamWriter = new SAMFileWriterFactory()
                        .makeSAMOrBAMWriter(samReader.getFileHeader(), true, inputBam);
        ) {
            for (SAMRecord samRecord : samReader) {
                bamWriter.addAlignment(samRecord);
            }
            bamWriter.close();

            SamReader bamReader = SamReaderFactory.makeDefault().referenceSequence(Defaults.REFERENCE_FASTA)
                    .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
                    .open(inputBam);

            BAMIndexer.createIndex(
                    bamReader,
                    new File(inputBam.getAbsolutePath() + ".bai"));

            bamReader.close();
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }

        return inputBam;
    }
}
