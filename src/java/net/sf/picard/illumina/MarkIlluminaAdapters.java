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

package net.sf.picard.illumina;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.sam.ReservedTagConstants;
import net.sf.picard.util.*;
import net.sf.samtools.*;
import net.sf.samtools.util.SequenceUtil;
import net.sf.samtools.util.StringUtil;

import java.io.File;
import java.util.Iterator;

/**
 * Command line program to mark the location of adapter sequences.
 * This also outputs a histogram of metrics describing the clipped bases
 *
 * @author Tim Fennell (adapted by mborkan@broadinstitute.org)
 */
public class MarkIlluminaAdapters extends CommandLineProgram {

    // The following attributes define the command-line arguments
    @Usage
    public String USAGE =
            getStandardUsagePreamble() +  "Reads a SAM or BAM file and rewrites it with new adapter-trimming tags.\n" +
                    "Clear any existing adapter-trimming tags (XT:i:).\n" +
                    "Only works for unaligned files in query-name order.\n"+
                    "Note: This is a utility program and will not be run in the pipeline.\n";

    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;
    @Option(doc="If output is not specified, just the metrics are generated",
            shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional=true)
    public File OUTPUT;
    @Option(doc="Histogram showing counts of bases_clipped in how many reads", shortName="M")
    public File METRICS;
    @Option(doc="The minimum number of bases that must match the adapter that will be clipped. Defaults to " +
            ClippingUtility.MIN_MATCH_PE_BASES + " if paired-end, otherwise" + ClippingUtility.MIN_MATCH_BASES +
            "/nThe stricter match used when matching 2 reads will be twice this.",
            optional=true)
    public Integer MIN_MATCH_BASES;
    @Option(doc="The percentage of errors allowed when matching the adapter sequence. Defaults to " +
            ClippingUtility.MAX_PE_ERROR_RATE + " if paired-end, otherwise " + ClippingUtility.MAX_ERROR_RATE,
            optional=true)
    public Double MAX_ERROR_RATE;
    @Option(doc="Whether this is a paired-end run. ", shortName="PE")
    public Boolean PAIRED_RUN;
    @Option(doc="Which adapters to use, PAIRED_END, INDEXED, or SINGLE_END",
            mutex={"FIVE_PRIME_ADAPTER", "THREE_PRIME_ADAPTER"})
    // this probably only makes sense for paired_run where you need to specify either PAIRED_END or INDEXED?
    //                         or for non-paired_run where you need to specify either SINGLE_END or INDEXED?
    // but we won't enforce this.
    public IlluminaUtil.IlluminaAdapterPair ADAPTERS;

    @Option(doc="For specifying adapters other than standard Illumina", mutex = {"ADAPTERS"})
    public String FIVE_PRIME_ADAPTER;
    @Option(doc="For specifying adapters other than standard Illumina", mutex = {"ADAPTERS"})
    public String THREE_PRIME_ADAPTER;

    private static final Log log = Log.getInstance(MarkIlluminaAdapters.class);

    @Override
    protected String[] customCommandLineValidation() {
        // set default thresholds based on what kind of run
        if (PAIRED_RUN){
            if (MIN_MATCH_BASES == null) MIN_MATCH_BASES = ClippingUtility.MIN_MATCH_PE_BASES;
            if (MAX_ERROR_RATE == null) MAX_ERROR_RATE = ClippingUtility.MAX_PE_ERROR_RATE;
            // For paired runs, you may actually want to specify all 4 thresholds
            // so the stricter test when mismatch can be controlled.
            // We'll assume that the stricter test will be twice the min_match_bases
        } else {
            if (MIN_MATCH_BASES == null) MIN_MATCH_BASES = ClippingUtility.MIN_MATCH_BASES;
            if (MAX_ERROR_RATE == null) MAX_ERROR_RATE = ClippingUtility.MAX_ERROR_RATE;
        }
        return null;
    }

    public static void main(String[] args) {
        System.exit(new MarkIlluminaAdapters().instanceMain(args));
    }

    @Override
    protected int doWork() {
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(METRICS);

        SAMFileReader in = new SAMFileReader(INPUT);
        SAMFileWriter out = null;
        if (OUTPUT != null) {
            IoUtil.assertFileIsWritable(OUTPUT);
            out = new SAMFileWriterFactory().makeSAMOrBAMWriter(in.getFileHeader(), true, OUTPUT);
        }

        Histogram<Integer> histo = new Histogram<Integer>("clipped_bases", "read_count");

        // check sort order in the header - must be queryName for paired end runs
        if (PAIRED_RUN && !in.getFileHeader().getSortOrder().equals(SAMFileHeader.SortOrder.queryname)) {
            throw new PicardException("Input BAM file must be sorted by queryname");
        }

        final AdapterPair adapters;
        if (ADAPTERS != null) {
            adapters = ADAPTERS;
        } else {
            adapters = new CustomAdapterPair(FIVE_PRIME_ADAPTER, THREE_PRIME_ADAPTER);
        }
        // The following loop is roughly the same as "for (SAMRecord rec : in){"
        final ProgressLogger progress = new ProgressLogger(log, 1000000, "Read");
        for (Iterator<SAMRecord> iter = in.iterator(); iter.hasNext();) {
            SAMRecord rec = iter.next();

            //  clear any existing trim on rec
            rec.setAttribute(ReservedTagConstants.XT, null);

            SAMRecord rec2 = null;
            if (PAIRED_RUN) {
                if (rec.getFirstOfPairFlag() || rec.getSecondOfPairFlag()) {
                    // the secondOfPair should be the next record
                    rec2 = iter.hasNext() ? iter.next() : null;
                    if (rec2 == null) {
                        throw new PicardException("Missing second read for " + rec);
                    }

                    // clear any existing trim on rec2
                    rec2.setAttribute(ReservedTagConstants.XT, null);
                    if (!rec.getReadName().equals(rec2.getReadName())){
                        throw new PicardException("read names of two paired reads differs : " +
                                rec.getReadName() + ", " + rec2.getReadName());
                    }

                    // establish which of pair is first and which second
                    SAMRecord firstRead;
                    SAMRecord secondRead;
                    if (rec.getFirstOfPairFlag()){
                        firstRead = rec;
                        secondRead = rec2;
                    } else {
                        firstRead = rec2;
                        secondRead = rec;
                    }
                    if (!firstRead.getFirstOfPairFlag()){
                        throw new PicardException("first of two reads doesn't have getFirstOfPairFlag()");
                    }
                    if (!secondRead.getSecondOfPairFlag()){
                        throw new PicardException("second of two reads doesn't have getSecondOfPairFlag()");
                    }

                    String warnString = ClippingUtility.adapterTrimIlluminaPairedReads(firstRead, secondRead,
                            adapters, MIN_MATCH_BASES, MAX_ERROR_RATE);
                    if (warnString != null) {
                        log.info("Adapter trimming " + warnString);
                    }
                } else {
                    throw new PicardException("Non-paired reads in a paired run " + rec);
                }
            } else { // not a paired run
                ClippingUtility.adapterTrimIlluminaSingleRead(rec,
                        adapters, MIN_MATCH_BASES, MAX_ERROR_RATE);
            }

            if (out != null) out.addAlignment(rec);
            if (out != null && rec2 != null) out.addAlignment(rec2);

            Integer trimPoint = rec.getIntegerAttribute(ReservedTagConstants.XT);
            if (trimPoint != null) {
                histo.increment(rec.getReadLength() - trimPoint + 1);
            }

            progress.record(rec);
        }

        if (out != null) out.close();

        MetricsFile<?,Integer> metricsFile = getMetricsFile();
        metricsFile.setHistogram(histo);
        metricsFile.write(METRICS);

        return 0;
    }

    private class CustomAdapterPair implements AdapterPair {

        final String fivePrime, threePrime, fivePrimeReadOrder;
        final byte[]  fivePrimeBytes, threePrimeBytes, fivePrimeReadOrderBytes;

        private CustomAdapterPair(final String fivePrime, final String threePrime) {
            this.threePrime = threePrime;
            this.threePrimeBytes = StringUtil.stringToBytes(threePrime);

            this.fivePrime = fivePrime;
            this.fivePrimeReadOrder = SequenceUtil.reverseComplement(fivePrime);
            this.fivePrimeBytes = StringUtil.stringToBytes(fivePrime);
            this.fivePrimeReadOrderBytes = StringUtil.stringToBytes(fivePrimeReadOrder);
        }

        public String get3PrimeAdapter(){ return threePrime; }
        public String get5PrimeAdapter(){ return fivePrime; }
        public String get3PrimeAdapterInReadOrder(){ return threePrime; }
        public String get5PrimeAdapterInReadOrder() { return fivePrimeReadOrder; }
        public byte[] get3PrimeAdapterBytes() { return threePrimeBytes; }
        public byte[] get5PrimeAdapterBytes() { return fivePrimeBytes; }
        public byte[] get3PrimeAdapterBytesInReadOrder() { return threePrimeBytes; }
        public byte[] get5PrimeAdapterBytesInReadOrder()  { return fivePrimeReadOrderBytes; }
        public String getName() { return "Custom adapter pair"; }
    }
}
