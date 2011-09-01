/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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
package net.sf.picard.analysis;

import net.sf.picard.annotation.RefFlatReader.RefFlatColumns;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.samtools.*;
import net.sf.samtools.util.StringUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;

public class CollectRnaSeqMetricsTest {
    @Test
    public void basic() throws Exception {
        final String sequence = "chr1";
        final String ignoredSequence = "chrM";

        // Create a refFlat file with a single gene containing two exons, one of which is overlapped by the
        // ribosomal interval.
        final String[] refFlatFields = new String[RefFlatColumns.values().length];
        refFlatFields[RefFlatColumns.GENE_NAME.ordinal()] = "myGene";
        refFlatFields[RefFlatColumns.TRANSCRIPT_NAME.ordinal()] = "myTranscript";
        refFlatFields[RefFlatColumns.CHROMOSOME.ordinal()] = sequence;
        refFlatFields[RefFlatColumns.STRAND.ordinal()] = "+";
        refFlatFields[RefFlatColumns.TX_START.ordinal()] = "49";
        refFlatFields[RefFlatColumns.TX_END.ordinal()] = "500";
        refFlatFields[RefFlatColumns.CDS_START.ordinal()] = "74";
        refFlatFields[RefFlatColumns.CDS_END.ordinal()] = "400";
        refFlatFields[RefFlatColumns.EXON_COUNT.ordinal()] = "2";
        refFlatFields[RefFlatColumns.EXON_STARTS.ordinal()] = "49,249";
        refFlatFields[RefFlatColumns.EXON_ENDS.ordinal()] = "200,500";

        final File refFlatFile = File.createTempFile("tmp.", ".refFlat");
        refFlatFile.deleteOnExit();
        final PrintStream refFlatStream = new PrintStream(refFlatFile);
        refFlatStream.println(StringUtil.join("\t", refFlatFields));
        refFlatStream.close();

        // Create some alignments that hit the ribosomal sequence, various parts of the gene, and intergenic.
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
        // Set seed so that strandedness is consistent among runs.
        builder.setRandomSeed(0);
        final int sequenceIndex = builder.getHeader().getSequenceIndex(sequence);
        builder.addPair("pair1", sequenceIndex, 45, 475);
        builder.addPair("pair2", sequenceIndex, 90, 225);
        builder.addPair("pair3", sequenceIndex, 120, 600);
        builder.addFrag("frag1", sequenceIndex, 150, true);
        builder.addFrag("frag2", sequenceIndex, 450, true);
        builder.addFrag("frag3", sequenceIndex, 225, false);
        builder.addPair("rrnaPair", sequenceIndex, 400, 500);

        builder.addFrag("ignoredFrag", builder.getHeader().getSequenceIndex(ignoredSequence), 1, false);

        final File samFile = File.createTempFile("tmp.collectRnaSeqMetrics.", ".sam");
        final SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(builder.getHeader(), false, samFile);
        for (final SAMRecord rec: builder.getRecords()) samWriter.addAlignment(rec);
        samWriter.close();

        // Create an interval list with one ribosomal interval.
        final Interval rRnaInterval = new Interval(sequence, 300, 520, true, "rRNA");
        final IntervalList rRnaIntervalList = new IntervalList(builder.getHeader());
        rRnaIntervalList.add(rRnaInterval);
        final File rRnaIntervalsFile = File.createTempFile("tmp.rRna.", ".interval_list");
        rRnaIntervalsFile.deleteOnExit();
        rRnaIntervalList.write(rRnaIntervalsFile);

        // Generate the metrics.
        final File metricsFile = File.createTempFile("tmp.", ".rna_metrics");

        final int ret = new CollectRnaSeqMetrics().instanceMain(new String[]{
                "INPUT=" +               samFile.getAbsolutePath(),
                "OUTPUT=" +              metricsFile.getAbsolutePath(),
                "REF_FLAT=" +            refFlatFile.getAbsolutePath(),
                "RIBOSOMAL_INTERVALS=" + rRnaIntervalsFile.getAbsolutePath(),
                "STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND",
                "IGNORE_SEQUENCE=" +ignoredSequence
        });
        Assert.assertEquals(ret, 0);

        final MetricsFile<RnaSeqMetrics, Comparable<?>> output = new MetricsFile<RnaSeqMetrics, Comparable<?>>();
        output.read(new FileReader(metricsFile));

        final RnaSeqMetrics metrics = output.getMetrics().get(0);
        Assert.assertEquals(metrics.PF_ALIGNED_BASES, 396);
        Assert.assertEquals(metrics.PF_BASES, 432);
        Assert.assertEquals(metrics.RIBOSOMAL_BASES.longValue(), 108L);
        Assert.assertEquals(metrics.CODING_BASES, 136);
        Assert.assertEquals(metrics.UTR_BASES, 51);
        Assert.assertEquals(metrics.INTRONIC_BASES, 50);
        Assert.assertEquals(metrics.INTERGENIC_BASES, 51);
        Assert.assertEquals(metrics.CORRECT_STRAND_READS, 3);
        Assert.assertEquals(metrics.INCORRECT_STRAND_READS, 4);
        Assert.assertEquals(metrics.IGNORED_READS, 1);
    }
}
