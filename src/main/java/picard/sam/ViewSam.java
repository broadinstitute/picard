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

package picard.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.AsciiWriter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamRecordIntervalIteratorFactory;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;

/**
 * Prints a SAM or BAM file to the screen.
 *
 * <p>Very simple command that just reads a SAM or BAM file and writes out the header
 * and each record to standard out. When an (optional) intervals file is specified,
 * only records overlapping those intervals will be output. </p>
 *
 * <p>All reads, just the aligned reads, or just the unaligned reads can be printed out by
 * setting AlignmentStatus accordingly. The SAM or BAM header can be printed out separately
 * using HEADER_ONLY. Only the alignment records can be printed using RECORDS_ONLY.
 * However, HEADER_ONLY and RECORDS_ONLY cannot both be specified at one time.</p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li> A SAM or BAM file to be viewed </li>
 *     <li> Optional arguments specifying which reads or records need to be viewed </li>
 * </ul>
 *
 * <p>
 * <h4>Usage example: </h4>
 * <pre>
 *     java -jar picard.jar ViewSam \
 *          I=input_reads.bam \
 *          HEADER_ONLY=true
 * </pre>
 * <p>
 * <br/>
 *
 * @author tfennell@broad.mit.edu
 */
@CommandLineProgramProperties(
        summary = ViewSam.USAGE_DETAILS,
        oneLineSummary = ViewSam.USAGE_SUMMARY,
        programGroup = DiagnosticsAndQCProgramGroup.class)
@DocumentedFeature
public class ViewSam extends CommandLineProgram {

    static final String USAGE_SUMMARY = "Prints a SAM or BAM file to the screen";
    static final String USAGE_DETAILS = "Very simple command that just reads a SAM or BAM file and" +
            "writes out the header and each record to standard out. When an (optional) intervals" +
            "file is specified, only records overlapping those intervals will be output.\n"+
            "All reads, just the aligned reads, or just the unaligned reads can be printed out by"+
            "setting AlignmentStatus accordingly. The SAM or BAM header can be printed out separately"+
            "using HEADER_ONLY. Only the alignment records can be printed using RECORDS_ONLY."+
            "However, HEADER_ONLY and RECORDS_ONLY cannot both be specified at one time."+
            "<p>"+
            "<h4>Usage example: </h4>" +
            "<pre>"+
            "java -jar picard.jar ViewSam  <br />" +
            "      I=sample.bam  <br />" +
            "      HEADER_ONLY=true" +
            "</pre>";

    public static enum AlignmentStatus {Aligned, Unaligned, All}

    public static enum PfStatus {PF, NonPF, All}

    public final String USAGE = getStandardUsagePreamble() + "Prints a SAM or BAM file to the screen.";
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The SAM or BAM file or GA4GH url to view.")
    public String INPUT;

    @Argument(doc = "Print out all reads, just the aligned reads or just the unaligned reads.")
    public AlignmentStatus ALIGNMENT_STATUS = AlignmentStatus.All;

    @Argument(doc = "Print out all reads, just the PF reads or just the non-PF reads.")
    public PfStatus PF_STATUS = PfStatus.All;

    @Argument(doc="Print the SAM header only.", optional = true)
    public boolean HEADER_ONLY = false;

    @Argument(doc="Print the alignment records only.", optional = true)
    public boolean RECORDS_ONLY = false;

    @Argument(doc = "An intervals file used to restrict what records are output.", optional = true)
    public File INTERVAL_LIST;

    @Override
    protected int doWork() {
        return writeSamText(System.out);
    }

    @Override
    protected String[] customCommandLineValidation() {
        if (HEADER_ONLY && RECORDS_ONLY) {
            return new String[]{"Cannot specify both HEADER_ONLY=true and RECORDS_ONLY=true."};
        }
        return super.customCommandLineValidation();
    }


    /**
     * This is factored out of doWork only for unit testing.
     */
    int writeSamText(PrintStream printStream) {
        try {
            final CloseableIterator<SAMRecord> samRecordsIterator;
            final SamReader samReader = SamReaderFactory.makeDefault()
                    .referenceSequence(REFERENCE_SEQUENCE)
                    .open(SamInputResource.of(INPUT));

            // if we are only using the header or we aren't using intervals, then use the reader as the iterator.
            // otherwise use the SamRecordIntervalIteratorFactory to make an interval-ed iterator
            if (HEADER_ONLY || INTERVAL_LIST == null) {
                samRecordsIterator = samReader.iterator();
            } else {
                IOUtil.assertFileIsReadable(INTERVAL_LIST);

                final List<Interval> intervals = IntervalList.fromFile(INTERVAL_LIST).uniqued().getIntervals();
                samRecordsIterator = new SamRecordIntervalIteratorFactory().makeSamRecordIntervalIterator(samReader, intervals, samReader.hasIndex());
            }
            final AsciiWriter writer = new AsciiWriter(printStream);
            final SAMFileHeader header = samReader.getFileHeader();
            if (!RECORDS_ONLY) {
                if (header.getTextHeader() != null) {
                    writer.write(header.getTextHeader());
                } else {
                    // Headers that are too large are not retained as text, so need to regenerate text
                    new SAMTextHeaderCodec().encode(writer, header, true);
                }
            }
            if (!HEADER_ONLY) {
                while (samRecordsIterator.hasNext()) {
                    final SAMRecord rec = samRecordsIterator.next();

                    if (printStream.checkError()) {
                        return 1;
                    }

                    if (this.ALIGNMENT_STATUS == AlignmentStatus.Aligned && rec.getReadUnmappedFlag()) continue;
                    if (this.ALIGNMENT_STATUS == AlignmentStatus.Unaligned && !rec.getReadUnmappedFlag()) continue;

                    if (this.PF_STATUS == PfStatus.PF && rec.getReadFailsVendorQualityCheckFlag()) continue;
                    if (this.PF_STATUS == PfStatus.NonPF && !rec.getReadFailsVendorQualityCheckFlag()) continue;
                    writer.write(rec.getSAMString());
                }
            }
            writer.flush();
            if (printStream.checkError()) {
                return 1;
            }
            CloserUtil.close(writer);
            CloserUtil.close(samRecordsIterator);
            return 0;
        } catch (IOException e) {
            throw new PicardException("Exception writing SAM text", e);
        }
    }
}
