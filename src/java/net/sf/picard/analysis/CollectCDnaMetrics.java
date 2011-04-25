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

import net.sf.picard.analysis.directed.CDnaMetrics;
import net.sf.picard.annotation.Gene;
import net.sf.picard.annotation.GeneAnnotationReader;
import net.sf.picard.annotation.LocusFunction;
import net.sf.picard.annotation.Transcript;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.Interval;
import net.sf.picard.util.Log;
import net.sf.picard.util.OverlapDetector;
import net.sf.samtools.*;
//import net.sf.picard.analysis.LocusFunction.*;

import java.io.File;
import java.util.Collection;
import java.util.List;

/**
 * Program to collect metrics about the alignment of cDNA to the various functional classes of loci in the genome:
 * coding, intronic, UTR, intragenic.
 */
public class CollectCDnaMetrics  extends SinglePassSamProgram {
    private static final Log LOG = Log.getInstance(CollectCDnaMetrics.class);
    @Option(doc="Gene annotations in refFlat form")
    public File REF_FLAT;

    @Option(shortName = StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME)
    public File SEQUENCE_DICTIONARY;

    @Option
    public boolean STRIP_LEADING_CHR_IN_REF_FLAT = true;

    private OverlapDetector<Gene> overlapDetector;
    private SAMSequenceDictionary sequenceDictionary;

    private CDnaMetrics metrics;

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        new CollectCDnaMetrics().instanceMainWithExit(argv);
    }

    @Override
    protected void setup(SAMFileHeader header, File samFile) {
        sequenceDictionary = new SAMFileReader(SEQUENCE_DICTIONARY).getFileHeader().getSequenceDictionary();
        overlapDetector = GeneAnnotationReader.loadRefFlat(REF_FLAT, sequenceDictionary, STRIP_LEADING_CHR_IN_REF_FLAT);
        LOG.info("Loaded " + overlapDetector.getAll().size() + " genes.");
        metrics = new CDnaMetrics();
    }

    @Override
    protected void acceptRead(SAMRecord rec, ReferenceSequence ref) {
        if (rec.getReadUnmappedFlag() || rec.getReadFailsVendorQualityCheckFlag()) return;
        final Interval readInterval = new Interval(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd());
        final Collection<Gene> overlappingGenes = overlapDetector.getOverlaps(readInterval);
        List<AlignmentBlock> alignmentBlocks = rec.getAlignmentBlocks();
        for (final AlignmentBlock alignmentBlock : alignmentBlocks) {
            final LocusFunction[] locusFunctions = new LocusFunction[alignmentBlock.getLength()];
            for (int i = 0; i < locusFunctions.length; ++i) locusFunctions[i] = LocusFunction.INTRAGENIC;
            for (final Gene gene : overlappingGenes) {
                for (final Transcript transcript : gene) {
                    transcript.getLocusFunctionForRange(alignmentBlock.getReferenceStart(), locusFunctions);
                }
            }
            for (int i = 0; i < locusFunctions.length; ++i) {
                ++metrics.ALIGNED_PF_BASES;
                switch (locusFunctions[i]) {
                    case INTRAGENIC: ++metrics.INTRAGENIC_BASES; break;
                    case INTRONIC:   ++metrics.INTRONIC_BASES;   break;
                    case UTR:        ++metrics.UTR_BASES;        break;
                    case CODING:     ++metrics.CODING_BASES;     break;
                }
            }
        }
    }

    @Override
    protected void finish() {
        if (metrics.ALIGNED_PF_BASES > 0) {
            metrics.PCT_CODING_BASES =     metrics.CODING_BASES     / (double) metrics.ALIGNED_PF_BASES;
            metrics.PCT_UTR_BASES =        metrics.UTR_BASES        / (double) metrics.ALIGNED_PF_BASES;
            metrics.PCT_INTRONIC_BASES =   metrics.INTRONIC_BASES   / (double) metrics.ALIGNED_PF_BASES;
            metrics.PCT_INTRAGENIC_BASES = metrics.INTRAGENIC_BASES / (double) metrics.ALIGNED_PF_BASES;
        }
        final MetricsFile<CDnaMetrics, Integer> file = getMetricsFile();
        file.addMetric(metrics);
        file.write(OUTPUT);
    }
}
