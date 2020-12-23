/*
 * The MIT License
 *
 * Copyright (c) 2012 The Broad Institute
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
package picard.sam;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.AsciiWriter;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.LineReader;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.List;

public class ViewSamTest extends CommandLineProgramTest {
    @Override
    public String getCommandLineProgramName() {
        return ViewSam.class.getSimpleName();
    }

    /**
     * Confirm that ViewSam retains whatever version number was in the input header.
     */
    @Test
    public void testHeaderVersion() throws Exception {
        final String oldVersionHeader = "@HD\tVN:1.3\tSO:unsorted";
        final File inputSam = File.createTempFile("ViewSamTest.input.", ".sam");
        inputSam.deleteOnExit();
        final AsciiWriter writer = new AsciiWriter(new FileOutputStream(inputSam));
        writer.write(oldVersionHeader);
        writer.write("\n");
        writer.close();
        final File viewSamOutputFile = File.createTempFile("ViewSamTest.output.", ".sam");
        viewSamOutputFile.deleteOnExit();

        final ViewSam viewSam = new ViewSam();
        viewSam.INPUT = inputSam.getAbsolutePath();

        // create a print stream to this file
        final PrintStream viewSamPrintStream = new PrintStream(viewSamOutputFile);
        // make sure the command line call exited successfully
        Assert.assertEquals(viewSam.writeSamText(viewSamPrintStream), 0);
        viewSamPrintStream.close();

        final LineReader viewSamInputReader = new BufferedLineReader(new FileInputStream(viewSamOutputFile));
        Assert.assertEquals(viewSamInputReader.readLine(), oldVersionHeader);
    }

    /**
     * Confirm that ViewSam only outputs records that overlap intervals in a provided interval file.
     */
    @Test
    public void testIntervals() throws Exception {
        // a SAM file designed to test intervals against
        final File inputSam = new File("testdata/picard/sam/viewsam_intervals_test.sam");
        // an interval file containing the intervals to run against the SAM
        final File inputIntervalsFile = new File("testdata/picard/sam/viewsam_intervals_test.interval_list");

        // create temp output file that ViewSam call get written to
        final File viewSamOutputFile = File.createTempFile("ViewSamTest.output.", ".sam");
        viewSamOutputFile.deleteOnExit();

        final ViewSam viewSam = new ViewSam();
        viewSam.INPUT = inputSam.getAbsolutePath();
        viewSam.INTERVAL_LIST = inputIntervalsFile;

        // create a print stream to this file
        final PrintStream viewSamPrintStream = new PrintStream(viewSamOutputFile);
        // make sure the command line call exited successfully
        Assert.assertEquals(viewSam.writeSamText(viewSamPrintStream), 0);
        viewSamPrintStream.close();

        // load the interval file
        final IntervalList inputIntervalsList = IntervalList.fromFile(inputIntervalsFile);
        // ViewSam internally utilizes uniqued intervals, so we will compare to the same
        final List<Interval> intervals = inputIntervalsList.uniqued().getIntervals();

        // make a reader that is not using intervals to load the output file we wrote that
        // was written by the call to ViewSam with the given interval file.  This will give us
        // the "filtered" file that we can compare to the intervals and ensure that only
        // overlapped records were written
        final SamReader samReader = SamReaderFactory.makeDefault().open(viewSamOutputFile);

        // make sure the intervals file caused at least one match to be found
        boolean foundMatches = false;

        for (final SAMRecord samRecord : samReader) {
            // make an interval representing this SAM record
            final Interval samRecordInterval = new Interval(samRecord.getContig(), samRecord.getStart(), samRecord.getEnd());
            // go through and look to see whether this SAM interval overlaps a filtering interval
            boolean samRecordIntervalOverlaps = false;
            for (final Interval interval : intervals) {
                if (interval.intersects(samRecordInterval)) {
                    samRecordIntervalOverlaps = true;
                    // mark that we have found at least one SAM record that overlaps an interval
                    foundMatches = true;
                    break;
                }
            }
            // if this SAM record does not overlap an interval, it should not have been written
            Assert.assertTrue(samRecordIntervalOverlaps, "SAM record written out was not overlapped by an interval.");
        }

        // we should have at least one SAM record written to ensure interval filtering worked correctly
        Assert.assertTrue(foundMatches, "No SAM records overlapped the given intervals.");
    }
}
