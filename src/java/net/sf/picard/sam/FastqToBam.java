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
import net.sf.picard.util.ReadableQualityFormatType;
import net.sf.picard.util.SolexaQualityConverter;
import net.sf.samtools.*;
import net.sf.samtools.util.StringUtil;

import java.io.File;

/**
 * Converts a fastq file to an unaligned BAM/SAM format.
 * See <a href="http://maq.sourceforge.net/fastq.shtml">MAQ FastQ specification</a> for details.
 * Three fastq versions are supported: FastqSanger, FastqSolexa and FastqIllumina.
 * Input files can be in GZip format (end in .gz).
 */
public class FastqToBam extends CommandLineProgram {

    @Usage(programVersion="1.0") 
    public String USAGE = "Extracts read sequences and qualities from the input fastq file and writes them into the output file in unaligned BAM format."
        + " Input files can be in GZip format (end in .gz).\n" 
        + "Following quality formats are supported: \n"
        + "  FastqSanger - refers to Sanger style FASTQ files which encode PHRED qualities using an ASCII offset of 33 \n"
        + "  FastqSolexa - refers to (early) Solexa/Illumina style FASTQ files (from Illumina pipeline version prior to 1.3) which encode Solexa qualities using an ASCII offset of 64 \n"
        + "  FastqIllumina  - refers to recent Solexa/Illumina style FASTQ files (from Illumina pipeline version 1.3+) which encode PHRED qualities using an ASCII offset of 64"
        ;

    @Option(shortName="F", doc="Input file (fastq.gz) to extract reads from (single-end fastq or, if paired, first end of the pair fastq).")
    public File FASTQ ;

    @Option(shortName="F2", doc="Input file (fastq.gz) to extract reads from (if paired, second end of the pair fastq).", optional=true) 
    public File SECOND_END_FASTQ ;

    @Option(doc="Output SAM/BAM file. ", shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME) 
    public File OUTPUT ;

    @Option(shortName="RB", doc="Run barcode")
    public String RUN_BARCODE;

    @Option(shortName="RG", doc="Read group name")
    public String READ_GROUP_NAME;

    @Option(shortName="SM", doc="Sample name")
    public String SAMPLE_NAME;

    @Option(shortName="V", doc="FASTQ version")
    public ReadableQualityFormatType FASTQ_VERSION ;

    private static final SolexaQualityConverter solexaQualityConverter = SolexaQualityConverter.getSingleton();

    public static void main(final String[] argv) {
        System.exit(new FastqToBam().instanceMain(argv));
    }

    protected int doWork() {
        final int readCount = (SECOND_END_FASTQ == null) ?  doUnpaired() : doPaired();
        System.out.println("Processed "+readCount+" fastq reads");
        return 0;
    }


    protected int doUnpaired() {
        IoUtil.assertFileIsReadable(FASTQ);
        IoUtil.assertFileIsWritable(OUTPUT);
        
        final FastqReader freader = new FastqReader(FASTQ);
        final SAMFileHeader header = createFileHeader() ;
        final SAMFileWriter writer = (new SAMFileWriterFactory()).makeSAMOrBAMWriter(header, false, OUTPUT);

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

    protected int doPaired() {
        IoUtil.assertFileIsReadable(FASTQ);
        IoUtil.assertFileIsReadable(SECOND_END_FASTQ);
        IoUtil.assertFileIsWritable(OUTPUT);
        
        final FastqReader freader1 = new FastqReader(FASTQ);
        final FastqReader freader2 = new FastqReader(SECOND_END_FASTQ);
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

        if (freader1.hasNext() || freader2.hasNext()) {
            throw new PicardException("Input paired fastq files must be the same length");
        }

        writer.close();
        return readCount;
    }

    private SAMRecord createSamRecord(final SAMFileHeader header, final String baseName, final FastqRecord frec,
                                      final boolean paired) {
        final SAMRecord srec = new SAMRecord(header);
        srec.setReadName(RUN_BARCODE + ":" + baseName);
        srec.setReadString(frec.getReadString());
        srec.setReadUmappedFlag(true);
        srec.setAttribute(ReservedTagConstants.READ_GROUP_ID, READ_GROUP_NAME);
        final byte[] quals = StringUtil.stringToBytes(frec.getBaseQualityString());
        convertQuality(quals, FASTQ_VERSION);
        srec.setBaseQualities(quals);
        if (paired) {
            srec.setReadPairedFlag(true);
            srec.setMateUnmappedFlag(true);
        }
        return srec ;
    }

    private SAMFileHeader createFileHeader() {
        final SAMReadGroupRecord rgroup = new SAMReadGroupRecord(READ_GROUP_NAME);
        rgroup.setSample(SAMPLE_NAME);
        final SAMFileHeader header = new SAMFileHeader();
        header.addReadGroup(rgroup);
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        return header ;
    }

    void convertQuality(final byte[] quals, final ReadableQualityFormatType version) {

        switch (version)  {
            case FastqSanger : 
                SAMUtils.fastqToPhred(quals);
                break ;
            case FastqSolexa : {
                solexaQualityConverter.convertSolexaQualityCharsToPhredBinary(quals);
                break ;
                }
            case FastqIllumina  : {
                solexaQualityConverter.convertSolexa_1_3_QualityCharsToPhredBinary(quals);
                break ;
                }
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

    private String error(final FastqReader freader, final String str) {
        return str +" at line "+freader.getLineNumber() +" in file "+freader.getFile().getAbsolutePath();
    }

    // Read names cannot contain blanks
    private String getReadName(final String fastaqHeader) {
        final int idx = fastaqHeader.indexOf(" ");
        return (idx == -1) ? fastaqHeader : fastaqHeader.substring(0,idx); 
    }
}
