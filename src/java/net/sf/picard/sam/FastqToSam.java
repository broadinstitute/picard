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
package net.sf.picard.sam;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.FastqQualityFormat;
import net.sf.picard.util.Log;
import net.sf.picard.util.SolexaQualityConverter;
import net.sf.samtools.*;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.util.StringUtil;

import java.io.File;

/**
 * Converts a fastq file to an unaligned BAM/SAM format.
 * See <a href="http://maq.sourceforge.net/fastq.shtml">MAQ FastQ specification</a> for details.
 * Three fastq versions are supported: FastqSanger, FastqSolexa and FastqIllumina.
 * Input files can be in GZip format (end in .gz).
 */
public class FastqToSam extends CommandLineProgram {
    private static Log LOG = Log.getInstance(FastqToSam.class);

    @Usage 
    public String USAGE = "Extracts read sequences and qualities from the input fastq file and writes them into the output file in unaligned BAM format."
        + " Input files can be in GZip format (end in .gz).\n" 
        ;

    @Option(shortName="F1", doc="Input fastq file (optionally gzipped) for single end data, or first read in paired end data.")
    public File FASTQ;

    @Option(shortName="F2", doc="Input fastq file (optionally gzipped) for the second read of paired end data.", optional=true)
    public File FASTQ2;

    @Option(shortName="V", doc="A value describing how the quality values are encoded in the fastq.  Either Solexa for pre-pipeline 1.3 " +
            "style scores (solexa scaling + 66), Illumina for pipeline 1.3 and above (phred scaling + 64) or Standard for phred scaled " +
            "scores with a character shift of 33.")
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

    @Option(shortName="SO", doc="The sort order for the output sam/bam file.")
    public SortOrder SORT_ORDER = SortOrder.queryname;

    private static final SolexaQualityConverter solexaQualityConverter = SolexaQualityConverter.getSingleton();

    /** Stock main method. */
    public static void main(final String[] argv) {
        System.exit(new FastqToSam().instanceMain(argv));
    }

    /* Simply invokes the right method for unpaired or paired data. */
    protected int doWork() {
        final int readCount = (FASTQ2 == null) ?  doUnpaired() : doPaired();
        LOG.info("Processed " + readCount + " fastq reads");
        return 0;
    }

    /** Creates a simple SAM file from a single fastq file. */
    protected int doUnpaired() {
        IoUtil.assertFileIsReadable(FASTQ);
        IoUtil.assertFileIsWritable(OUTPUT);
        
        final FastqReader freader = new FastqReader(FASTQ);
        final SAMFileHeader header = createFileHeader();
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, OUTPUT);

        int readCount = 0;
        for ( ; freader.hasNext()  ; readCount++) {
            final FastqRecord frec = freader.next();
            final SAMRecord srec = createSamRecord(header, getReadName(frec.getReadHeader()) , frec, false) ;
            srec.setReadPairedFlag(false);
            writer.addAlignment(srec);
        }

        writer.close();
        return readCount;
    }

    /** More complicated method that takes two fastq files and builds pairing information in the SAM. */
    protected int doPaired() {
        IoUtil.assertFileIsReadable(FASTQ);
        IoUtil.assertFileIsReadable(FASTQ2);
        IoUtil.assertFileIsWritable(OUTPUT);
        
        final FastqReader freader1 = new FastqReader(FASTQ);
        final FastqReader freader2 = new FastqReader(FASTQ2);
        final SAMFileHeader header = createFileHeader() ;
        final SAMFileWriter writer = (new SAMFileWriterFactory()).makeSAMOrBAMWriter(header, false, OUTPUT);

        int readCount = 0;
        for ( ; freader1.hasNext() && freader2.hasNext() ; readCount++) {
            final FastqRecord frec1 = freader1.next();
            final FastqRecord frec2 = freader2.next();

            final String frec1Name = getReadName(frec1.getReadHeader());
            final String frec2Name = getReadName(frec2.getReadHeader());
            final String baseName = getBaseName(frec1Name, frec2Name, freader1, freader2);

            final SAMRecord srec1 = createSamRecord(header, baseName, frec1, true) ;
            srec1.setFirstOfPairFlag(true);
            srec1.setSecondOfPairFlag(false);
            writer.addAlignment(srec1);

            final SAMRecord srec2 = createSamRecord(header, baseName, frec2, true) ;
            srec2.setFirstOfPairFlag(false);
            srec2.setSecondOfPairFlag(true);
            writer.addAlignment(srec2);
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

        final SAMFileHeader header = new SAMFileHeader();
        header.addReadGroup(rgroup);
        header.setSortOrder(this.SORT_ORDER);
        return header ;
    }

    /** Based on the type of quality scores coming in, converts them to a numeric byte[] in prhred scale. */
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
     * <li> A pair read name must contain at least one /
     * <li> First pair read name must end with "/1 and second pair read name ends with "/2" 
     * <li> The baseName (read name part before the /) must be the same for both read names
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

        if (StringUtil.isBlank(num1)) throw new PicardException(error(freader1,"Pair 1 number is required ("+readName1+")"));
        if (StringUtil.isBlank(num2)) throw new PicardException(error(freader2,"Pair 2 number is required ("+readName2+")"));
        if (!num1.equals("1")) throw new PicardException(error(freader1,"Pair 1 number must be 1 ("+readName1+")"));
        if (!num2.equals("2")) throw new PicardException(error(freader2,"Pair 2 number must be 2 ("+readName2+")"));

        return baseName1 ;
    }

    /** Breaks up read name into baseName and number separated by the last / */
    private String [] getReadNameTokens(final String readName, final int pairNum, final FastqReader freader) {
        final int idx = readName.lastIndexOf("/");
        if (idx == -1) throw new PicardException(error(freader,"Pair read name "+pairNum+" must have a slash: "+readName));
        final String result[] = new String[2];
        result[0] = readName.substring(0,idx); // baseName
        result[1] = readName.substring(idx+1, readName.length()); // number
        return result ;
    }

    /** Little utility to give error messages corresponding to line numbers in the input files. */
    private String error(final FastqReader freader, final String str) {
        return str +" at line "+freader.getLineNumber() +" in file "+freader.getFile().getAbsolutePath();
    }

    // Read names cannot contain blanks
    private String getReadName(final String fastaqHeader) {
        final int idx = fastaqHeader.indexOf(" ");
        return (idx == -1) ? fastaqHeader : fastaqHeader.substring(0,idx); 
    }
}
