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

import htsjdk.samtools.ReservedTagConstants;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.FastqQualityFormat;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Iso8601Date;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.QualityEncodingDetector;
import htsjdk.samtools.util.SolexaQualityConverter;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Converts a fastq file to an unaligned BAM/SAM format.
 * See <a href="http://maq.sourceforge.net/fastq.shtml">MAQ FastQ specification</a> for details.
 * Three fastq versions are supported: FastqSanger, FastqSolexa and FastqIllumina.
 * Input files can be in GZip format (end in .gz).
 */
@CommandLineProgramProperties(
        usage = "Extracts read sequences and qualities from the input fastq file and writes them into the output file in unaligned BAM format."
                + " Input files can be in GZip format (end in .gz).\n",
        usageShort = "Converts a fastq file to an unaligned BAM or SAM file",
        programGroup = SamOrBam.class
)
public class FastqToSam extends CommandLineProgram {
    private static final Log LOG = Log.getInstance(FastqToSam.class);

    @Option(shortName="F1", doc="Input fastq file (optionally gzipped) for single end data, or first read in paired end data.")
    public File FASTQ;

    @Option(shortName="F2", doc="Input fastq file (optionally gzipped) for the second read of paired end data.", optional=true)
    public File FASTQ2;

    @Option(shortName="V", doc="A value describing how the quality values are encoded in the fastq.  Either Solexa for pre-pipeline 1.3 " +
            "style scores (solexa scaling + 66), Illumina for pipeline 1.3 and above (phred scaling + 64) or Standard for phred scaled " +
            "scores with a character shift of 33.  If this value is not specified, the quality format will be detected automatically.", optional = true)
    public FastqQualityFormat QUALITY_FORMAT;

    @Option(doc="Output SAM/BAM file. ", shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME) 
    public File OUTPUT ;

    @Option(shortName="RG", doc="Read group name")
    public String READ_GROUP_NAME = "A";

    @Option(shortName="SM", doc="Sample name to insert into the read group header")
    public String SAMPLE_NAME;

    @Option(shortName="LB", doc="The library name to place into the LB attribute in the read group header", optional=true)
    public String LIBRARY_NAME;

    @Option(shortName="PU", doc="The platform unit (often run_barcode.lane) to insert into the read group header", optional=true)
    public String PLATFORM_UNIT;

    @Option(shortName="PL", doc="The platform type (e.g. illumina, solid) to insert into the read group header", optional=true)
    public String PLATFORM;

    @Option(shortName="CN", doc="The sequencing center from which the data originated", optional=true)
    public String SEQUENCING_CENTER;

    @Option(shortName = "PI", doc = "Predicted median insert size, to insert into the read group header", optional = true)
    public Integer PREDICTED_INSERT_SIZE;

    @Option(doc="Comment(s) to include in the merged output file's header.", optional=true, shortName="CO")
    public List<String> COMMENT = new ArrayList<String>();

    @Option(shortName = "DS", doc = "Inserted into the read group header", optional = true)
    public String DESCRIPTION;

    @Option(shortName = "DT", doc = "Date the run was produced, to insert into the read group header", optional = true)
    public Iso8601Date RUN_DATE;

    @Option(shortName="SO", doc="The sort order for the output sam/bam file.")
    public SortOrder SORT_ORDER = SortOrder.queryname;

    @Option(doc="Minimum quality allowed in the input fastq.  An exception will be thrown if a quality is less than this value.")
    public int MIN_Q = 0;

    @Option(doc="Maximum quality allowed in the input fastq.  An exception will be thrown if a quality is greater than this value.")
    public int MAX_Q = SAMUtils.MAX_PHRED_SCORE;

    @Option(doc="If true and this is an unpaired fastq any occurance of '/1' will be removed from the end of a read name.")
    public Boolean STRIP_UNPAIRED_MATE_NUMBER = false;

    @Option(doc="Allow (and ignore) empty lines")
    public Boolean ALLOW_AND_IGNORE_EMPTY_LINES = false;

    private static final SolexaQualityConverter solexaQualityConverter = SolexaQualityConverter.getSingleton();

    /** Stock main method. */
    public static void main(final String[] argv) {
        System.exit(new FastqToSam().instanceMain(argv));
    }

    /* Simply invokes the right method for unpaired or paired data. */
    protected int doWork() {
            final QualityEncodingDetector detector = new QualityEncodingDetector();
            final FastqReader reader = new FastqReader(FASTQ,ALLOW_AND_IGNORE_EMPTY_LINES);
            if (FASTQ2 == null) {
                detector.add(QualityEncodingDetector.DEFAULT_MAX_RECORDS_TO_ITERATE, reader);
            } else {
                final FastqReader reader2 = new FastqReader(FASTQ2,ALLOW_AND_IGNORE_EMPTY_LINES);
                detector.add(QualityEncodingDetector.DEFAULT_MAX_RECORDS_TO_ITERATE, reader, reader2);
                reader2.close();
            }
            reader.close();
            
            QUALITY_FORMAT = detector.generateBestGuess(QualityEncodingDetector.FileContext.FASTQ, QUALITY_FORMAT);
            if (detector.isDeterminationAmbiguous())
                LOG.warn("Making ambiguous determination about fastq's quality encoding; more than one format possible based on observed qualities.");
            LOG.info(String.format("Auto-detected quality format as: %s.", QUALITY_FORMAT));

        final int readCount = (FASTQ2 == null) ?  doUnpaired() : doPaired();
        LOG.info("Processed " + readCount + " fastq reads");
        return 0;
    }

    /** Creates a simple SAM file from a single fastq file. */
    protected int doUnpaired() {
        IOUtil.assertFileIsReadable(FASTQ);
        IOUtil.assertFileIsWritable(OUTPUT);
        
        final FastqReader freader = new FastqReader(FASTQ,ALLOW_AND_IGNORE_EMPTY_LINES);
        final SAMFileHeader header = createFileHeader();
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, OUTPUT);

        int readCount = 0;
        final ProgressLogger progress = new ProgressLogger(LOG);
        for ( ; freader.hasNext()  ; readCount++) {
            final FastqRecord frec = freader.next();
            final SAMRecord srec = createSamRecord(header, getReadName(frec.getReadHeader(), false) , frec, false) ;
            srec.setReadPairedFlag(false);
            writer.addAlignment(srec);
            progress.record(srec);
        }

        writer.close();
        return readCount;
    }

    /** More complicated method that takes two fastq files and builds pairing information in the SAM. */
    protected int doPaired() {
        IOUtil.assertFileIsReadable(FASTQ);
        IOUtil.assertFileIsReadable(FASTQ2);
        IOUtil.assertFileIsWritable(OUTPUT);
        
        final FastqReader freader1 = new FastqReader(FASTQ,ALLOW_AND_IGNORE_EMPTY_LINES);
        final FastqReader freader2 = new FastqReader(FASTQ2,ALLOW_AND_IGNORE_EMPTY_LINES);
        final SAMFileHeader header = createFileHeader() ;
        final SAMFileWriter writer = (new SAMFileWriterFactory()).makeSAMOrBAMWriter(header, false, OUTPUT);

        int readCount = 0;
        final ProgressLogger progress = new ProgressLogger(LOG);
        for ( ; freader1.hasNext() && freader2.hasNext() ; readCount++) {
            final FastqRecord frec1 = freader1.next();
            final FastqRecord frec2 = freader2.next();

            final String frec1Name = getReadName(frec1.getReadHeader(), true);
            final String frec2Name = getReadName(frec2.getReadHeader(), true);
            final String baseName = getBaseName(frec1Name, frec2Name, freader1, freader2);

            final SAMRecord srec1 = createSamRecord(header, baseName, frec1, true) ;
            srec1.setFirstOfPairFlag(true);
            srec1.setSecondOfPairFlag(false);
            writer.addAlignment(srec1);
            progress.record(srec1);

            final SAMRecord srec2 = createSamRecord(header, baseName, frec2, true) ;
            srec2.setFirstOfPairFlag(false);
            srec2.setSecondOfPairFlag(true);
            writer.addAlignment(srec2);
            progress.record(srec2);
        }

        writer.close();

        if (freader1.hasNext() || freader2.hasNext()) {
            throw new PicardException("Input paired fastq files must be the same length");
        }

        return readCount;
    }

    private SAMRecord createSamRecord(final SAMFileHeader header, final String baseName, final FastqRecord frec, final boolean paired) {
        final SAMRecord srec = new SAMRecord(header);
        srec.setReadName(baseName);
        srec.setReadString(frec.getReadString());
        srec.setReadUnmappedFlag(true);
        srec.setAttribute(ReservedTagConstants.READ_GROUP_ID, READ_GROUP_NAME);
        final byte[] quals = StringUtil.stringToBytes(frec.getBaseQualityString());
        convertQuality(quals, QUALITY_FORMAT);
        for (final byte qual : quals) {
            final int uQual = qual & 0xff;
            if (uQual < MIN_Q || uQual > MAX_Q) {
                throw new PicardException("Base quality " + uQual + " is not in the range " + MIN_Q + ".." +
                MAX_Q + " for read " + frec.getReadHeader());
            }
        }
        srec.setBaseQualities(quals);

        if (paired) {
            srec.setReadPairedFlag(true);
            srec.setMateUnmappedFlag(true);
        }
        return srec ;
    }

    /** Creates a simple header with the values provided on the command line. */
    private SAMFileHeader createFileHeader() {
        final SAMReadGroupRecord rgroup = new SAMReadGroupRecord(this.READ_GROUP_NAME);
        rgroup.setSample(this.SAMPLE_NAME);
        if (this.LIBRARY_NAME != null) rgroup.setLibrary(this.LIBRARY_NAME);
        if (this.PLATFORM != null) rgroup.setPlatform(this.PLATFORM);
        if (this.PLATFORM_UNIT != null) rgroup.setPlatformUnit(this.PLATFORM_UNIT);
        if (this.SEQUENCING_CENTER != null) rgroup.setSequencingCenter(SEQUENCING_CENTER);
        if (this.PREDICTED_INSERT_SIZE != null) rgroup.setPredictedMedianInsertSize(PREDICTED_INSERT_SIZE);
        if (this.DESCRIPTION != null) rgroup.setDescription(this.DESCRIPTION);
        if (this.RUN_DATE != null) rgroup.setRunDate(this.RUN_DATE);

        final SAMFileHeader header = new SAMFileHeader();
        header.addReadGroup(rgroup);

        for (final String comment : COMMENT) {
            header.addComment(comment);
        }

        header.setSortOrder(this.SORT_ORDER);
        return header ;
    }

    /** Based on the type of quality scores coming in, converts them to a numeric byte[] in phred scale. */
    void convertQuality(final byte[] quals, final FastqQualityFormat version) {
        switch (version)  {
            case Standard:
                SAMUtils.fastqToPhred(quals);
                break ;
            case Solexa:
                solexaQualityConverter.convertSolexaQualityCharsToPhredBinary(quals);
                break ;
            case Illumina:
                solexaQualityConverter.convertSolexa_1_3_QualityCharsToPhredBinary(quals);
                break ;
            }
    }

    /** Returns read baseName and asserts correct pair read name format:
     * <ul>
     * <li> Paired reads must either have the exact same read names or they must contain at least one "/"
     * <li> and the First pair read name must end with "/1" and second pair read name ends with "/2"
     * <li> The baseName (read name part before the /) must be the same for both read names
     * <li> If the read names are exactly the same but end in "/2" or "/1" then an exception will be thrown 
     * </ul>
     */
    String getBaseName(final String readName1, final String readName2, final FastqReader freader1, final FastqReader freader2) {
        String [] toks = getReadNameTokens(readName1, 1, freader1);
        final String baseName1 = toks[0] ;
        final String num1 = toks[1] ;

        toks = getReadNameTokens(readName2, 2, freader2);
        final String baseName2 = toks[0] ;
        final String num2 = toks[1];

        if (!baseName1.equals(baseName2)) {
            throw new PicardException(String.format("In paired mode, read name 1 (%s) does not match read name 2 (%s)", baseName1,baseName2));
        }

        final boolean num1Blank = StringUtil.isBlank(num1);
        final boolean num2Blank = StringUtil.isBlank(num2);
        if (num1Blank || num2Blank) {
            if(!num1Blank) throw new PicardException(error(freader1,"Pair 1 number is missing (" +readName1+ "). Both pair numbers must be present or neither."));       //num1 != blank and num2   == blank
            else if(!num2Blank) throw new PicardException(error(freader2, "Pair 2 number is missing (" +readName2+ "). Both pair numbers must be present or neither.")); //num1 == blank and num =2 != blank 
        } else {
            if (!num1.equals("1")) throw new PicardException(error(freader1,"Pair 1 number must be 1 ("+readName1+")"));
            if (!num2.equals("2")) throw new PicardException(error(freader2,"Pair 2 number must be 2 ("+readName2+")"));
        }

        return baseName1 ;
    }

    /** Breaks up read name into baseName and number separated by the last / */
    private String [] getReadNameTokens(final String readName, final int pairNum, final FastqReader freader) {
        if(readName.equals("")) throw new PicardException(error(freader,"Pair read name "+pairNum+" cannot be empty: "+readName));

        final int idx = readName.lastIndexOf("/");
        final String result[] = new String[2];

        if (idx == -1) {
            result[0] = readName;
            result[1] = null;
        } else {
            result[1] = readName.substring(idx+1, readName.length()); // should be a 1 or 2
            
            if(!result[1].equals("1") && !result[1].equals("2")) {    //if not a 1 or 2 then names must be identical
                result[0] = readName;
                result[1] = null;
            }
            else {
                result[0] = readName.substring(0,idx); // baseName
            }
        }

        return result ;
    }

    /** Little utility to give error messages corresponding to line numbers in the input files. */
    private String error(final FastqReader freader, final String str) {
        return str +" at line "+freader.getLineNumber() +" in file "+freader.getFile().getAbsolutePath();
    }

    // Read names cannot contain blanks
    private String getReadName(final String fastqHeader, final boolean paired) {
        final int idx = fastqHeader.indexOf(" ");
        String readName = (idx == -1) ? fastqHeader : fastqHeader.substring(0,idx);

        // NOTE: the while loop isn't necessarily the most efficient way to handle this but we don't
        // expect this to ever happen more than once, just trapping pathological cases
        while (STRIP_UNPAIRED_MATE_NUMBER && !paired && (readName.endsWith("/1") || readName.endsWith("/2"))) {
            // If this is an unpaired run we want to make sure that "/1" isn't tacked on the end of the read name,
            // as this can cause problems down the road in MergeBamAlignment
            readName = readName.substring(0, readName.length() - 2);
        }

        return readName;
    }

    @Override
    protected String[] customCommandLineValidation() {
        if (MIN_Q < 0) return new String[]{"MIN_Q must be >= 0"};
        if (MAX_Q > SAMUtils.MAX_PHRED_SCORE) return new String[]{"MAX_Q must be <= " + SAMUtils.MAX_PHRED_SCORE};
        return null;
    }
}
