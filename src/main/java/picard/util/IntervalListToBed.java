/*
 * The MIT License
 *
 * Copyright (c) 2016 Tim Fennell
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
package picard.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalCodec;
import htsjdk.samtools.util.IntervalCoordinateComparator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.IntervalList.IntervalListCodec;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.IntervalsManipulationProgramGroup;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Trivially simple command line program to convert an IntervalList file to a BED file.
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = "Converts an Picard IntervalList file to a BED file.",
        oneLineSummary = "Converts an Picard IntervalList file to a BED file.",
        programGroup = IntervalsManipulationProgramGroup.class
)
@DocumentedFeature
public class IntervalListToBed extends CommandLineProgram {
    @Argument(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input IntervalList file.")
    public File INPUT;

    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output BED file.")
    public File OUTPUT;

    @Argument(doc="The score, between 0-1000, to output for each interval in the BED file.")
    public int SCORE = 500;

    @Argument(doc="If true, sort the interval list prior to outputting as BED file.")
    public boolean SORT = true;

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        try (final AbstractFeatureReader<Interval, LineIterator> intervalsReader = AbstractFeatureReader.getFeatureReader(INPUT.getPath(), new IntervalListCodec(), false);
             CloseableIterator<Interval> intervalsToUse = getSortedIntervals(intervalsReader);
             final BufferedWriter out = IOUtil.openFileForBufferedWriting(OUTPUT)) {
            intervalsToUse.stream().forEach(i -> generateOutput(i, out));
            return 0;
        } catch (Exception e) {
            throw new RuntimeIOException(e);
        }
    }

    private CloseableIterator<Interval> getSortedIntervals(AbstractFeatureReader<Interval, LineIterator> intervalsReader)  throws IOException {
        if (SORT) {
            final SAMFileHeader header = (SAMFileHeader) intervalsReader.getHeader();
            final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
            final IntervalCoordinateComparator comparator = new IntervalCoordinateComparator(header);
            SortingCollection<Interval> sortedIntervals =
                    SortingCollection.newInstance(Interval.class,
                            new IntervalCodec(sequenceDictionary),
                            comparator, 500000, TMP_DIR.stream().map(File::toPath).toArray(java.nio.file.Path[]::new));

            for (final Interval i: intervalsReader.iterator()) {
                sortedIntervals.add(i);
            }
            return sortedIntervals.iterator();
        } else {
            return intervalsReader.iterator();
        }
    }

    private void generateOutput(Interval interval, BufferedWriter outputWriter) {
        try {
            final String strand = interval.isNegativeStrand() ? "-" : "+";
            final List<?> fields = CollectionUtil.makeList(interval.getContig(), interval.getStart() - 1, interval.getEnd(), interval.getName(), SCORE, strand);
            outputWriter.append(fields.stream().map(String::valueOf).collect(Collectors.joining("\t")));
            outputWriter.newLine();
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }
}
