/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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

package picard.fingerprint;

import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import static java.lang.Math.max;

/**
 * Class that implements  SexInferencer for an input type of BAM. The idea here is to extract coverage information from the index directly
 */
public class SexInferenceEngineFromBAM extends SexInferenceEngine {

    private List<File> bamList;

    public SexInferenceEngineFromBAM(final Set<String> MALE_CHROMS, final Set<String> FEMALE_CHROMS, final List<File> bamList){
        super(MALE_CHROMS, FEMALE_CHROMS);
        this.bamList = bamList;
    }

    @Override
    protected final List<SampleXy> getSexChromCoverageDensity() {
        final List<File> bams = IOUtil.unrollFiles(bamList, BamFileIoUtils.BAM_FILE_EXTENSION);
        final List<SampleXy> samples = new ArrayList<>();
        double maxX = 0, maxY = 0;

        for (final File f : bams) {
            final SamReader sam = SamReaderFactory.make().open(f);
            final SAMSequenceDictionary dict = sam.getFileHeader().getSequenceDictionary();
            if (!sam.hasIndex())
                throw new IllegalArgumentException("Input Sequence files must have an index.");

            final BAMIndex index = sam.indexing().getIndex();
            final String sample = getSampleName(sam);

            double maleChromSize = 0, femaleChromSize = 0, maleMappedReads = 0, femaleMappedReads = 0;

            // Compute the major values needed from the BAM index, including:
            //   1) The total number of mapped reads in the BAM across all chromosomes
            //   2) The number of mapped reads to female chromosomes and to male chromosomes
            //   3) The size (in bases) of the female and male chromosomes
            long totalMapped = 0;
            for (final SAMSequenceRecord seq : dict.getSequences()) {
                final String name = seq.getSequenceName();
                final BAMIndexMetaData metadata = index.getMetaData(seq.getSequenceIndex());
                final int mappedReads = metadata.getAlignedRecordCount();
                totalMapped += mappedReads;

                if (FEMALE_CHROMS.contains(name)) {
                    femaleChromSize += seq.getSequenceLength();
                    femaleMappedReads += metadata.getAlignedRecordCount();
                } else if (MALE_CHROMS.contains(name)) {
                    maleChromSize += seq.getSequenceLength();
                    maleMappedReads += metadata.getAlignedRecordCount();
                }
            }

            try {
                sam.close();
            } catch (final IOException e) {
                throw new RuntimeIOException(e);
            }

            // Calculate the rpkm values and cache the highest ones we see for X and Y
            final double xRpkm = femaleMappedReads / (femaleChromSize / 1000) / (totalMapped / 1000000);
            final double yRpkm =   maleMappedReads / (  maleChromSize / 1000) / (totalMapped / 1000000);
            maxX = max(maxX, xRpkm);
            maxY = max(maxY, yRpkm);

            final SampleXy yx = new SampleXy(sample, xRpkm, yRpkm);
            samples.add(yx);
        }

        for (final SampleXy sample : samples) {
            //The normalization of coverage is done by looking at the x and y chromes only.
            //Since we expect the maximal coverage of X to be a 2x coverage while that on y to be a 1x coverage,
            //we divide the normalization by 2 so that the y coverage is normalized to 1/2 rather than to 1. This means that
            //on average coverage on y and coverage on x will have the same effect.
            sample.setxDensity(sample.getxDensity() / maxX);
            sample.setyDensity(sample.getyDensity() / maxY / 2);
        }
        return samples;
    }

    @Override
    protected double[][] getCentroids() {return new double[][]{{1.0, 0}, {0.5, 0.5}};}

    /**
     * Gets the sample name from the read group headers.  Throws an exception if there is not one and only one unique
     * name in the read group headers.
     */
    static String getSampleName(final SamReader in) {
        String name = null;
        for (final SAMReadGroupRecord rg : in.getFileHeader().getReadGroups()) {
            if (name == null) {
                name = rg.getSample() ;
            } else if (!name.equals(rg.getSample())) throw new IllegalStateException("BAM file with multiple samples not supported.");
        }

        if (name == null) throw new IllegalStateException("BAM files without sample names in RG headers not supported.");

        return name;
    }
}
