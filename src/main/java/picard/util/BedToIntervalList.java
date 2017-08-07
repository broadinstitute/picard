package picard.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Intervals;

import java.io.File;
import java.io.IOException;

/**
 * @author nhomer
 */
@CommandLineProgramProperties(
        summary = BedToIntervalList.USAGE_SUMMARY + BedToIntervalList.USAGE_DETAILS,
        oneLineSummary = BedToIntervalList.USAGE_SUMMARY,
        programGroup = Intervals.class
)
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

    final Log LOG = Log.getInstance(getClass());

    // Stock main method
    public static void main(final String[] args) {
        new BedToIntervalList().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(SEQUENCE_DICTIONARY);
        IOUtil.assertFileIsWritable(OUTPUT);
        try {
            // create a new header that we will assign the dictionary provided by the SAMSequenceDictionaryExtractor to.
            final SAMFileHeader header = new SAMFileHeader();
            final SAMSequenceDictionary samSequenceDictionary = SAMSequenceDictionaryExtractor.extractDictionary(SEQUENCE_DICTIONARY);
            header.setSequenceDictionary(samSequenceDictionary);
            // set the sort order to be sorted by coordinate, which is actually done below
            // by getting the .uniqued() intervals list before we write out the file
            header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
            final IntervalList intervalList = new IntervalList(header);

            /**
             * NB: BED is zero-based, but a BEDCodec by default (since it is returns tribble Features) has an offset of one,
             * so it returns 1-based starts.  Ugh.  Set to zero.
             */
            final FeatureReader<BEDFeature> bedReader = AbstractFeatureReader.getFeatureReader(INPUT.getAbsolutePath(), new BEDCodec(BEDCodec.StartOffset.ZERO), false);
            final CloseableTribbleIterator<BEDFeature> iterator = bedReader.iterator();
            final ProgressLogger progressLogger = new ProgressLogger(LOG, (int) 1e6);

            while (iterator.hasNext()) {
                final BEDFeature bedFeature = iterator.next();
                final String sequenceName = bedFeature.getContig();
                /**
                 * NB: BED is zero-based, so we need to add one here to make it one-based.  Please observe we set the start
                 * offset to zero when creating the BEDCodec.
                 */
                final int start = bedFeature.getStart() + 1;
                /**
                 * NB: BED is 0-based OPEN (which, for the end is equivalent to 1-based closed).
                 */
                final int end = bedFeature.getEnd();
                // NB: do not use an empty name within an interval
                String name = bedFeature.getName();
                if (name.isEmpty()) name = null;

                final SAMSequenceRecord sequenceRecord = header.getSequenceDictionary().getSequence(sequenceName);

                // Do some validation
                if (null == sequenceRecord) {
                    throw new PicardException(String.format("Sequence '%s' was not found in the sequence dictionary", sequenceName));
                } else if (start < 1) {
                    throw new PicardException(String.format("Start on sequence '%s' was less than one: %d", sequenceName, start));
                } else if (sequenceRecord.getSequenceLength() < start) {
                    throw new PicardException(String.format("Start on sequence '%s' was past the end: %d < %d", sequenceName, sequenceRecord.getSequenceLength(), start));
                } else if (end < 1) {
                    throw new PicardException(String.format("End on sequence '%s' was less than one: %d", sequenceName, end));
                } else if (sequenceRecord.getSequenceLength() < end) {
                    throw new PicardException(String.format("End on sequence '%s' was past the end: %d < %d", sequenceName, sequenceRecord.getSequenceLength(), end));
                } else if (end < start - 1) {
                    throw new PicardException(String.format("On sequence '%s', end < start-1: %d <= %d", sequenceName, end, start));
                }

                final boolean isNegativeStrand = bedFeature.getStrand() == Strand.NEGATIVE;
                final Interval interval = new Interval(sequenceName, start, end, isNegativeStrand, name);
                intervalList.add(interval);

                progressLogger.record(sequenceName, start);
            }
            CloserUtil.close(bedReader);

            // Sort and write the output
            IntervalList out = intervalList;
            if (SORT) out = out.sorted();
            if (UNIQUE) out = out.uniqued();
            out.write(OUTPUT);

        } catch (final IOException e) {
            throw new RuntimeException(e);
        }

        return 0;
    }
}
