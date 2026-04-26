package picard.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.IntervalsManipulationProgramGroup;
import picard.util.IntervalFileReader.FormatDetectionResult;
import picard.util.IntervalFileReader.IntervalFileFormat;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

/**
 * @author nhomer
 */
@CommandLineProgramProperties(
        summary = BedToIntervalList.USAGE_SUMMARY + BedToIntervalList.USAGE_DETAILS,
        oneLineSummary = BedToIntervalList.USAGE_SUMMARY,
        programGroup = IntervalsManipulationProgramGroup.class
)
@DocumentedFeature
public class BedToIntervalList extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Converts a BED file to a Picard Interval List.  " ;
    static final String USAGE_DETAILS = "This tool provides easy conversion from BED to the Picard interval_list format which is " +
            "required by many Picard processing tools. Note that the coordinate system of BED files is such that the first base or " +
            "position in a sequence is numbered \"0\", while in interval_list files it is numbered \"1\"." +
            "<br /><br />" +
            "BED files contain sequence data displayed in a flexible format that includes nine optional fields, " +
            "in addition to three required fields within the annotation tracks. The required fields of a BED file include:" +
            "<pre>" +
            "     chrom - The name of the chromosome (e.g. chr20) or scaffold (e.g. scaffold10671) <br />"   +
            "     chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered \"0\" <br />"   +
            "     chromEnd - The ending position of the feature in the chromosome or scaffold.  The chromEnd base is not" +
            " included in the display of the feature. For example, the first 100 bases of a " +
            "chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99." +
            "</pre>" +
            "In each annotation track, the number of fields per line must be consistent throughout a data set. " +
            "For additional information regarding BED files and the annotation field options, please see:" +
            " http://genome.ucsc.edu/FAQ/FAQformat.html#format1." +
            "<br /> <br /> " +
            "Interval_list files contain sequence data distributed into intervals. The interval_list file format is relatively simple " +
            "and reflects the SAM alignment format to a degree.  A SAM style header must be present in the file that lists the sequence " +
            "records against which the intervals are described.  After the header, the file then contains records, one per line in plain " +
            "text format with the following values tab-separated::" +
            "<pre> " +
            "     -Sequence name (SN) - The name of the sequence in the file for identification purposes, can be chromosome number e.g. chr20 <br /> " +
            "     -Start position - Interval start position (starts at +1) <br /> " +
            "     -End position - Interval end position (1-based, end inclusive) <br /> " +
            "     -Strand - Indicates +/- strand for the interval (either + or -) <br /> " +
            "     -Interval name - (Each interval should have a unique name) " +
            "</pre>" +
            "<br/>" +
            "This tool requires a sequence dictionary, provided with the SEQUENCE_DICTIONARY or SD argument. " +
            "The value given to this argument can be any of the following:" +
            "<pre>" +
            "    - A file with .dict extension generated using Picard's CreateSequenceDictionaryTool</br>" +
            "    - A reference.fa or reference.fasta file with a reference.dict in the same directory</br>" +
            "    - Another IntervalList with @SQ lines in the header from which to generate a dictionary</br>" +
            "    - A VCF that contains #contig lines from which to generate a sequence dictionary</br>" +
            "    - A SAM or BAM file with @SQ lines in the header from which to generate a dictionary</br>" +
            "</pre>" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar BedToIntervalList \\<br />" +
            "      I=input.bed \\<br />" +
            "      O=list.interval_list \\<br />" +
            "      SD=reference_sequence.dict" +
            "</pre>" +
            "<br /> <br /> "+
            "<hr />"
            ;
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input BED file")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME,
            doc = "The sequence dictionary, or BAM/VCF/IntervalList from which a dictionary can be extracted.")
    public File SEQUENCE_DICTIONARY;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output Picard Interval List")
    public File OUTPUT;

    @Argument(doc="If true, sort the output interval list before writing it.")
    public boolean SORT = true;

    @Argument(doc="If true, unique the output interval list by merging overlapping regions, before writing it (implies sort=true).")
    public boolean UNIQUE = false;

    @Hidden
    @Argument(doc = "If true, entries that are on contig-names that are missing from the provided dictionary will be dropped.")
    public boolean DROP_MISSING_CONTIGS = false;

    @Argument(doc = "If true, write length zero intervals in input bed file to resulting interval list file.")
    public boolean KEEP_LENGTH_ZERO_INTERVALS = false;

    private final Log LOG = Log.getInstance(getClass());

    @Override
    protected int doWork() {
        // Only assert readability for regular files; FIFOs, named pipes, /dev/stdin,
        // and process-substitution paths (/dev/fd/N) all return false from isFile().
        if (INPUT.isFile()) {
            IOUtil.assertFileIsReadable(INPUT);
        }
        IOUtil.assertFileIsReadable(SEQUENCE_DICTIONARY);
        IOUtil.assertFileIsWritable(OUTPUT);

        try {
            final SAMFileHeader header = new SAMFileHeader();
            final SAMSequenceDictionary samSequenceDictionary = SAMSequenceDictionaryExtractor.extractDictionary(SEQUENCE_DICTIONARY.toPath());
            header.setSequenceDictionary(samSequenceDictionary);
            header.setSortOrder(SAMFileHeader.SortOrder.coordinate);

            // AbstractFeatureReader requires a real file path, so buffer /dev/stdin to a temp file.
            final File bedFile;
            File tempFile = null;
            if (INPUT.getPath().equals("/dev/stdin")) {
                tempFile = File.createTempFile("bed_to_interval_list_stdin_", ".bed");
                tempFile.deleteOnExit();
                Files.copy(System.in, tempFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                bedFile = tempFile;
            } else {
                bedFile = INPUT;
            }

            // Sniff the format before parsing; reject anything that isn't BED.
            try (final BufferedReader reader = new BufferedReader(new FileReader(bedFile))) {
                reader.mark(8 * 1024);
                final FormatDetectionResult detected = IntervalFileReader.detectIntervalFormat(reader);
                if (detected.format() != IntervalFileFormat.BED) {
                    final String hint = detected.format() == IntervalFileFormat.INTERVAL_LIST
                            ? " Input appears to be an interval_list file; supply a BED file instead."
                            : (detected.firstLine() != null ? " First data line: " + detected.firstLine() : " File appears to be empty or contain only headers.");
                    throw new PicardException("BedToIntervalList requires BED format input." + hint);
                }
            }

            IntervalList out = IntervalFileReader.fromBed(bedFile, header, DROP_MISSING_CONTIGS, KEEP_LENGTH_ZERO_INTERVALS);
            if (SORT) out = out.sorted();
            if (UNIQUE) out = out.uniqued();
            out.write(OUTPUT);
            LOG.info(String.format("Wrote %d intervals spanning a total of %d bases",
                    out.getIntervals().size(), out.getBaseCount()));
        } catch (final IOException e) {
            throw new RuntimeException(e);
        }

        return 0;
    }
}
