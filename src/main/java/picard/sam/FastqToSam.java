/*
 * The MIT License
 *
 * Copyright (c) 2009-2016 The Broad Institute
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
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqConstants.FastqExtensions;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.FastqQualityFormat;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Iso8601Date;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.QualityEncodingDetector;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SolexaQualityConverter;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import picard.nio.PicardHtsPath;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Converts a FASTQ file to an unaligned BAM or SAM file.
 * <p>
 *     Output read records will contain the original base calls and quality scores will be
 *     translated depending on the base quality score encoding: FastqSanger, FastqSolexa and FastqIllumina.
 * </p>
 * <p>
 *     There are also arguments to provide values for SAM header and read attributes that are not present in FASTQ
 *     (e.g see <code>RG</code> or <code>SM</code> below).
 * </p>
 * <h3>Inputs</h3>
 * <p>
 *     One FASTQ file name for single-end or two for pair-end sequencing input data.
 *     These files might be in gzip compressed format (when file name is ending with ".gz").
 * </p>
 * <p>
 *     Alternatively, for larger inputs you can provide a collection of FASTQ files indexed by their name (see <code>USE_SEQUENCIAL_FASTQ</code> for details below).
 * </p>
 * <p>
 *     By default, this tool will try to guess the base quality score encoding. However you can indicate it explicitly
 *     using the <code>QUALITY_FORMAT</code> argument.
 * </p>
 * <h3>Output</h3>
 * A single unaligned BAM or SAM file. By default, the records are sorted by query (read) name.
 * <h3>Usage examples</h3>
 *
 * <h4>Example 1:</h4>
 * <p>
 *     Single-end sequencing FASTQ file conversion. All reads are annotated
 *     as belonging to the "rg0013" read group that in turn is part of the sample "sample001".
 * </p>
 * <pre>
 * java -jar picard.jar FastqToSam \
 *      F1=input_reads.fastq \
 *      O=unaligned_reads.bam \
 *      SM=sample001 \
 *      RG=rg0013
 * </pre>
 * <h4>Example 2:</h4>
 * <p>
 *     Similar to example 1 above, but for paired-end sequencing.
 * </p>
 * <pre>
 * java -jar picard.jar FastqToSam \
 *      F1=forward_reads.fastq \
 *      F2=reverse_reads.fastq \
 *      O=unaligned_read_pairs.bam \
 *      SM=sample001 \
 *      RG=rg0013
 * </pre>
 */
@CommandLineProgramProperties(
        summary = "<p>" + FastqToSam.USAGE_SUMMARY + ".</p>" + FastqToSam.USAGE_DETAILS,
        oneLineSummary = FastqToSam.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature
public class FastqToSam extends CommandLineProgram {
    static final String USAGE_SUMMARY =
            "Converts a FASTQ file to an unaligned BAM or SAM file";
    static final String USAGE_DETAILS =
            "<p>Output read records will contain the original base calls and quality scores will be " +
                    "translated depending on the base quality score encoding: FastqSanger, FastqSolexa and FastqIllumina.</p>" +
        "<p>There are also arguments to provide values for SAM header and read attributes that are not present in FASTQ " +
        "(e.g see RG or SM below).</p>" +
        "<h3>Inputs</h3>" +
        "<p>One FASTQ file name for single-end or two for pair-end sequencing input data. " +
        "These files might be in gzip compressed format (when file name is ending with \".gz\").</p>" +
        "<p>Alternatively, for larger inputs you can provide a collection of FASTQ files indexed by their name " +
        "(see USE_SEQUENCIAL_FASTQ for details below).</p>" +
        "<p>By default, this tool will try to guess the base quality score encoding. However you can indicate it explicitly " +
        "using the QUALITY_FORMAT argument.</p>" +
        "<h3>Output</h3>" +
        "<p>A single unaligned BAM or SAM file. By default, the records are sorted by query (read) name.</p>" +
        "<h3>Usage examples</h3>" +
        "<h4>Example 1:</h4>" +
        "<p>Single-end sequencing FASTQ file conversion. All reads are annotated " +
        "as belonging to the \"rg0013\" read group that in turn is part of the sample \"sample001\".</p>" +
        "<pre>java -jar picard.jar FastqToSam \\\n" +
        "        F1=input_reads.fastq \\\n" +
        "        O=unaligned_reads.bam \\\n" +
        "        SM=sample001 \\\n" +
        "        RG=rg0013</pre>" +
        "<h4>Example 2:</h4>" +
        "<p>Similar to example 1 above, but for paired-end sequencing.</p>" +
        "<pre>java -jar picard.jar FastqToSam \\\n" +
        "       F1=forward_reads.fastq \\\n" +
        "       F2=reverse_reads.fastq \\\n" +
        "       O=unaligned_read_pairs.bam \\\n" +
        "       SM=sample001 \\\n" +
        "       RG=rg0013</pre><hr />";

    private static final Log LOG = Log.getInstance(FastqToSam.class);

    @Argument(shortName="F1", doc="Input fastq file (optionally gzipped) for single end data, or first read in paired end data.")
    public PicardHtsPath FASTQ;

    @Argument(shortName="F2", doc="Input fastq file (optionally gzipped) for the second read of paired end data.", optional=true)
    public PicardHtsPath FASTQ2;

    @Argument(doc="Use sequential fastq files with the suffix <prefix>_###.fastq or <prefix>_###.fastq.gz." +
            "The files should be named:\n" +
            "    <prefix>_001.<extension>, <prefix>_002.<extension>, ..., <prefix>_XYZ.<extension>\n" +
            " The base files should be:\n" +
            "    <prefix>_001.<extension>\n" +
            " An example would be:\n" +
            "    RUNNAME_S8_L005_R1_001.fastq\n" +
            "    RUNNAME_S8_L005_R1_002.fastq\n" +
            "    RUNNAME_S8_L005_R1_003.fastq\n" +
            "    RUNNAME_S8_L005_R1_004.fastq\n" +
            "RUNNAME_S8_L005_R1_001.fastq should be provided as FASTQ.", optional=true)
    public boolean USE_SEQUENTIAL_FASTQS = false;

    @Argument(shortName="V", doc="A value describing how the quality values are encoded in the input FASTQ file.  " +
            "Either Solexa (phred scaling + 66), Illumina (phred scaling + 64) or Standard (phred scaling + 33).  " +
            "If this value is not specified, the quality format will be detected automatically.", optional = true)
    public FastqQualityFormat QUALITY_FORMAT;

    @Argument(doc="Output BAM/SAM/CRAM file. ", shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT ;

    @Argument(shortName="RG", doc="Read group name")
    public String READ_GROUP_NAME = "A";

    @Argument(shortName="SM", doc="Sample name to insert into the read group header")
    public String SAMPLE_NAME;

    @Argument(shortName="LB", doc="The library name to place into the LB attribute in the read group header", optional=true)
    public String LIBRARY_NAME;

    @Argument(shortName="PU", doc="The platform unit (often run_barcode.lane) to insert into the read group header", optional=true)
    public String PLATFORM_UNIT;

    @Argument(shortName="PL", doc="The platform type (e.g. ILLUMINA, SOLID) to insert into the read group header", optional=true)
    public String PLATFORM;

    @Argument(shortName="CN", doc="The sequencing center from which the data originated", optional=true)
    public String SEQUENCING_CENTER;

    @Argument(shortName = "PI", doc = "Predicted median insert size, to insert into the read group header", optional = true)
    public Integer PREDICTED_INSERT_SIZE;

    @Argument(shortName = "PG", doc = "Program group to insert into the read group header.", optional=true)
    public String PROGRAM_GROUP;

    @Argument(shortName = "PM", doc = "Platform model to insert into the group header (free-form text providing further details of the platform/technology used)", optional=true)
    public String PLATFORM_MODEL;

    @Argument(doc="Comment(s) to include in the merged output file's header.", optional=true, shortName="CO")
    public List<String> COMMENT = new ArrayList<>();

    @Argument(shortName = "DS", doc = "Inserted into the read group header", optional = true)
    public String DESCRIPTION;

    @Argument(shortName = "DT", doc = "Date the run was produced, to insert into the read group header", optional = true)
    public Iso8601Date RUN_DATE;

    @Argument(shortName="SO", doc="The sort order for the output BAM/SAM/CRAM file.")
    public SortOrder SORT_ORDER = SortOrder.queryname;

    @Argument(doc="Minimum quality allowed in the input fastq.  An exception will be thrown if a quality is less than this value.")
    public int MIN_Q = 0;

    @Argument(doc="Maximum quality allowed in the input fastq.  An exception will be thrown if a quality is greater than this value.")
    public int MAX_Q = SAMUtils.MAX_PHRED_SCORE;

    @Deprecated
    @Argument(doc="Deprecated (No longer used). If true and this is an unpaired fastq any occurrence of '/1' or '/2' will be removed from the end of a read name.")
    public Boolean STRIP_UNPAIRED_MATE_NUMBER = false;

    @Argument(doc="Allow (and ignore) empty lines")
    public Boolean ALLOW_AND_IGNORE_EMPTY_LINES = false;

    private static final SolexaQualityConverter solexaQualityConverter = SolexaQualityConverter.getSingleton();

    /**
     * Looks at fastq input(s) and attempts to determine the proper quality format
     *
     * Closes the reader(s) by side effect
     *
     * @param reader1 The first fastq input
     * @param reader2 The second fastq input, if necessary. To not use this input, set it to null
     * @param expectedQuality If provided, will be used for sanity checking. If left null, autodetection will occur
     */
    public static FastqQualityFormat determineQualityFormat(final FastqReader reader1, final FastqReader reader2, final FastqQualityFormat expectedQuality) {
        final QualityEncodingDetector detector = new QualityEncodingDetector();

        if (reader2 == null) {
            detector.add(QualityEncodingDetector.DEFAULT_MAX_RECORDS_TO_ITERATE, reader1);
        } else {
            detector.add(QualityEncodingDetector.DEFAULT_MAX_RECORDS_TO_ITERATE, reader1, reader2);
            reader2.close();
        }

        reader1.close();

        final FastqQualityFormat qualityFormat =  detector.generateBestGuess(QualityEncodingDetector.FileContext.FASTQ, expectedQuality);
        if (detector.isDeterminationAmbiguous()) {
            LOG.warn("Making ambiguous determination about fastq's quality encoding; more than one format possible based on observed qualities.");
        }
        LOG.info(String.format("Auto-detected quality format as: %s.", qualityFormat));

        return qualityFormat;
    }


    /**
     * Get a list of FASTQs that are sequentially numbered based on the first (base) fastq.
     * The files should be named:
     *   <prefix>_001.<extension>, <prefix>_002.<extension>, ..., <prefix>_XYZ.<extension>
     * The base files should be:
     *   <prefix>_001.<extension>
     * An example would be:
     *   RUNNAME_S8_L005_R1_001.fastq
     *   RUNNAME_S8_L005_R1_002.fastq
     *   RUNNAME_S8_L005_R1_003.fastq
     *   RUNNAME_S8_L005_R1_004.fastq
     * where `baseFastq` is the first in that list.
     */
    protected static List<Path> getSequentialFileList(final Path baseFastq) {
        final List<Path> files = new ArrayList<>();
        files.add(baseFastq);

        // Find the correct extension used in the base FASTQ
        FastqExtensions fastqExtensions = null;
        String suffix = null; // store the suffix including the extension
        for (final FastqExtensions ext : FastqExtensions.values()) {
            suffix = "_001" + ext.getExtension();
            if (baseFastq.toString().endsWith(suffix)) {
                fastqExtensions = ext;
                break;
            }
        }
        if (null == fastqExtensions) {
            throw new PicardException(String.format("Could not parse the FASTQ extension (expected '_001' + '%s'): %s", Arrays.toString(FastqExtensions.values()), baseFastq));
        }

        // Find all the files
        for (int idx = 2; true; idx++) {
            String fastq = baseFastq.toAbsolutePath().toString();
            fastq = String.format("%s_%03d%s", fastq.substring(0, fastq.length() - suffix.length()), idx, fastqExtensions.getExtension());
            try {
                IOUtil.assertFileIsReadable(Paths.get(fastq));
            } catch (final SAMException e) { // the file is not readable, so do not continue
                break;
            }
            files.add(Paths.get(fastq));
        }

        return files;
    }

    /* Simply invokes the right method for unpaired or paired data. */
    protected int doWork() {
        IOUtil.assertFileIsReadable(FASTQ.toPath());
        if (FASTQ2 != null) {
            IOUtil.assertFileIsReadable(FASTQ2.toPath());
        }
        IOUtil.assertFileIsWritable(OUTPUT);

        final SAMFileHeader header = createSamFileHeader();
        final SAMFileWriter writer = new SAMFileWriterFactory().makeWriter(header, false, OUTPUT, REFERENCE_SEQUENCE);

        // Set the quality format
        QUALITY_FORMAT = FastqToSam.determineQualityFormat(fileToFastqReader(FASTQ.toPath()),
                (FASTQ2 == null) ? null : fileToFastqReader(FASTQ2.toPath()),
                QUALITY_FORMAT);

        // Lists for sequential files, but also used when not sequential
        final List<FastqReader> readers1 = new ArrayList<>();
        final List<FastqReader> readers2 = new ArrayList<>();

        if (USE_SEQUENTIAL_FASTQS) {
            // Get all the files
            for (final Path fastq : getSequentialFileList(FASTQ.toPath())) {
                readers1.add(fileToFastqReader(fastq));
            }
            if (null != FASTQ2) {
                for (final Path fastq : getSequentialFileList(FASTQ2.toPath())) {
                    readers2.add(fileToFastqReader(fastq));
                }
                if (readers1.size() != readers2.size()) {
                    throw new PicardException(String.format("Found %d files for FASTQ and %d files for FASTQ2.", readers1.size(), readers2.size()));
                }
            }
        }
        else {
            readers1.add(fileToFastqReader(FASTQ.toPath()));
            if (FASTQ2 != null) {
                readers2.add(fileToFastqReader(FASTQ2.toPath()));
            }
        }

        // Loop through the FASTQs
        for (int idx = 0; idx < readers1.size(); idx++) {
            makeItSo(readers1.get(idx),
                    (readers2.isEmpty()) ? null : readers2.get(idx),
                    writer);
        }

        // Close all the things
        for (final FastqReader reader : readers1) reader.close();
        for (final FastqReader reader : readers2) reader.close();
        writer.close();

        return 0;
    }

    /**
     * Handles the FastqToSam execution on the FastqReader(s).
     *
     * In some circumstances it might be useful to circumvent the command line based instantiation of this
     * class, however note that there is no handholding or guardrails to running in this manner.
     *
     * It is the caller's responsibility to close the reader(s)
     *
     * @param reader1 The FastqReader for the first fastq file
     * @param reader2 The second FastqReader if applicable. Pass in null if only using a single reader
     * @param writer The SAMFileWriter where the new SAM file is written
     *
     */
    public void makeItSo(final FastqReader reader1, final FastqReader reader2, final SAMFileWriter writer) {
        final int readCount = (reader2 == null) ?  doUnpaired(reader1, writer) : doPaired(reader1, reader2, writer);
        LOG.info("Processed " + readCount + " fastq reads");
    }

    /** Creates a simple SAM file from a single fastq file. */
    protected int doUnpaired(final FastqReader freader, final SAMFileWriter writer) {
        int readCount = 0;
        final ProgressLogger progress = new ProgressLogger(LOG);
        for ( ; freader.hasNext()  ; readCount++) {
            final FastqRecord frec = freader.next();
            final SAMRecord srec = createSamRecord(writer.getFileHeader(), SequenceUtil.getSamReadNameFromFastqHeader(frec.getReadHeader()) , frec, false) ;
            srec.setReadPairedFlag(false);
            writer.addAlignment(srec);
            progress.record(srec);
        }

        return readCount;
    }

    /** More complicated method that takes two fastq files and builds pairing information in the SAM. */
    protected int doPaired(final FastqReader freader1, final FastqReader freader2, final SAMFileWriter writer) {
        int readCount = 0;
        final ProgressLogger progress = new ProgressLogger(LOG);
        for ( ; freader1.hasNext() && freader2.hasNext() ; readCount++) {
            final FastqRecord frec1 = freader1.next();
            final FastqRecord frec2 = freader2.next();

            final String frec1Name = SequenceUtil.getSamReadNameFromFastqHeader(frec1.getReadHeader());
            final String frec2Name = SequenceUtil.getSamReadNameFromFastqHeader(frec2.getReadHeader());
            final String baseName = getBaseName(frec1Name, frec2Name, freader1, freader2);

            final SAMRecord srec1 = createSamRecord(writer.getFileHeader(), baseName, frec1, true) ;
            srec1.setFirstOfPairFlag(true);
            srec1.setSecondOfPairFlag(false);
            writer.addAlignment(srec1);
            progress.record(srec1);

            final SAMRecord srec2 = createSamRecord(writer.getFileHeader(), baseName, frec2, true) ;
            srec2.setFirstOfPairFlag(false);
            srec2.setSecondOfPairFlag(true);
            writer.addAlignment(srec2);
            progress.record(srec2);
        }

        if (freader1.hasNext() || freader2.hasNext()) {
            throw new PicardException("Input paired fastq files must be the same length");
        }

        return readCount;
    }

    private FastqReader fileToFastqReader(final Path path) throws PicardException {
        return new FastqReader(null, IOUtil.openFileForBufferedReading(path), ALLOW_AND_IGNORE_EMPTY_LINES);
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
    public SAMFileHeader createSamFileHeader() {
        final SAMReadGroupRecord rgroup = new SAMReadGroupRecord(this.READ_GROUP_NAME);
        rgroup.setSample(this.SAMPLE_NAME);
        if (this.LIBRARY_NAME != null) rgroup.setLibrary(this.LIBRARY_NAME);
        if (this.PLATFORM != null) rgroup.setPlatform(this.PLATFORM);
        if (this.PLATFORM_UNIT != null) rgroup.setPlatformUnit(this.PLATFORM_UNIT);
        if (this.SEQUENCING_CENTER != null) rgroup.setSequencingCenter(SEQUENCING_CENTER);
        if (this.PREDICTED_INSERT_SIZE != null) rgroup.setPredictedMedianInsertSize(PREDICTED_INSERT_SIZE);
        if (this.DESCRIPTION != null) rgroup.setDescription(this.DESCRIPTION);
        if (this.RUN_DATE != null) rgroup.setRunDate(this.RUN_DATE);
        if (this.PLATFORM_MODEL != null) rgroup.setPlatformModel(this.PLATFORM_MODEL);
        if (this.PROGRAM_GROUP != null) rgroup.setProgramGroup(this.PROGRAM_GROUP);

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

        final int idx = readName.lastIndexOf('/');
        final String[] result = new String[2];

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

    @Override
    protected String[] customCommandLineValidation() {
        if (MIN_Q < 0) return new String[]{"MIN_Q must be >= 0"};
        if (MAX_Q > SAMUtils.MAX_PHRED_SCORE) return new String[]{"MAX_Q must be <= " + SAMUtils.MAX_PHRED_SCORE};
        return null;
    }
}
