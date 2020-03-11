package picard.sam;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * <p> Extracts read sequences and qualities from the input SAM/BAM file and SAM/BAM tags and writes them into
 * output files in Sanger FASTQ format.
 * See <a href="http://maq.sourceforge.net/fastq.shtml">MAQ FASTQ specification</a> for details.
 * <br />
 * <h4>Usage example:</h4>
 * <pre>
 * java -jar picard.jar SamToFastqWithTags
 *     I=input.bam
 *     FASTQ=output.fastq
 *     SEQUENCE_TAG_GROUP="CR"
 *     QUALITY_TAG_GROUP="CY"
 *     SEQUENCE_TAG_GROUP="CB,UR"
 *     QUALITY_TAG_GROUP="CY,UY"
 * </pre>
 * <hr />
 */
@CommandLineProgramProperties(
        summary = SamToFastqWithTags.USAGE_SUMMARY + SamToFastqWithTags.USAGE_DETAILS,
        oneLineSummary = SamToFastqWithTags.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
public class SamToFastqWithTags extends SamToFastq {
    static final String USAGE_SUMMARY = "Converts a SAM or BAM file to FASTQ alongside FASTQs created from tags.";
    static final String USAGE_DETAILS = " Extracts read sequences and qualities from the input SAM/BAM file and SAM/BAM tags" +
            " and writes them into the output file in Sanger FASTQ format." +
            " See <a href=\"http://maq.sourceforge.net/fastq.shtml\">MAQ FASTQ specification</a> for details.<br /> <br />" +
            "The following example will create two FASTQs from tags.  One will be converted with the base sequence coming from " +
            "the \"CR\" tag and base quality from the \"CY\" tag.  The other fastq will be converted with the base sequence coming" +
            " from the \"CB\" and \"UR\" tags concatenated together with no separator (not specified on command line) with the base" +
            " qualities coming from the \"CY\" and \"UY\" tags concatenated together.  The two files will be named CR.fastq" +
            " and CB_UR.fastq." +
            "<br />" +
            "<pre>" +
            "java -jar picard.jar SamToFastqWithTags <br />" +
            "     I=input.bam<br />" +
            "     FASTQ=output.fastq<br />" +
            "     SEQUENCE_TAG_GROUP=CR<br />" +
            "     QUALITY_TAG_GROUP=CY<br />" +
            "     SEQUENCE_TAG_GROUP=\"CB,UR\"<br />" +
            "     QUALITY_TAG_GROUP=\"CY,UY\"" +
            "</pre>" +
            "<hr />";

    @Argument(shortName = "STG", doc = "List of comma separated tag values to extract from Input SAM/BAM to be used as read sequence", minElements = 1)
    public List<String> SEQUENCE_TAG_GROUP;

    @Argument(shortName = "QTG", doc = "List of comma separated tag values to extract from Input SAM/BAM to be used as read qualities", optional = true)
    public List<String> QUALITY_TAG_GROUP;

    @Argument(shortName = "SEP", doc = "List of any sequences (e.g. 'AACCTG`) to put in between each comma separated list of sequence tags in each SEQUENCE_TAG_GROUP (STG)", optional = true)
    public List<String> TAG_GROUP_SEPERATOR;

    @Argument(shortName = "GZOPTG", doc = "Compress output FASTQ files per Tag grouping using gzip and append a .gz extension to the file names.")
    public Boolean COMPRESS_OUTPUTS_PER_TAG_GROUP = false;

    private final Log log = Log.getInstance(SamToFastqWithTags.class);

    private static final String TAG_SPLIT_DEFAULT_SEP = "";
    private static final String TAG_SPLIT_QUAL = "~";

    private ArrayList<String[]> SPLIT_SEQUENCE_TAGS;
    private ArrayList<String[]> SPLIT_QUALITY_TAGS;
    private ArrayList<String> SPLIT_SEPARATOR_TAGS;

    @Override
    protected void initializeAdditionalWriters() {
        setupTagSplitValues();
    }

    @Override
    protected void handleAdditionalRecords(SAMRecord currentRecord, Map<SAMReadGroupRecord, List<FastqWriter>> tagWriters, SAMRecord read1, SAMRecord read2) {
        final List<FastqWriter> rgTagWriters = tagWriters.get(currentRecord.getReadGroup());
        if (currentRecord.getReadPairedFlag()) {
            if (read1 != null && read2 !=null) {
                writeTagRecords(read1, 1, rgTagWriters);
                writeTagRecords(read2, 2, rgTagWriters);
            }
        } else {
            writeTagRecords(currentRecord, null, rgTagWriters);
        }
    }

    @Override
    protected Map<SAMReadGroupRecord, List<FastqWriter>> generateAdditionalWriters(List<SAMReadGroupRecord> readGroups, FastqWriterFactory factory) {
        return generateTagWriters(readGroups, factory);
    }

    // generate writers
    private Map<SAMReadGroupRecord, List<FastqWriter>> generateTagWriters(final List<SAMReadGroupRecord> samReadGroupRecords,
                                                                          final FastqWriterFactory factory) {
        final Map<SAMReadGroupRecord, List<FastqWriter>> writerMap = new HashMap<>();

        if (!OUTPUT_PER_RG) {
            /* Prepare tag writers based on sequence tag groups provided in command line */

            final List<FastqWriter> tagFastqWriters = makeTagWriters(null, factory);

            writerMap.put(null, tagFastqWriters);
            for (final SAMReadGroupRecord rg : samReadGroupRecords) {
                writerMap.put(rg, tagFastqWriters);
            }
        } else {
            /* prepare tag writers based on readgroup names */
            for (final SAMReadGroupRecord rg : samReadGroupRecords) {
                final List<FastqWriter> tagWriters = makeTagWriters(rg, factory);

                writerMap.put(rg, tagWriters);
            }
        }
        return writerMap;
    }

    /**
     *     Creates fastq writers based on readgroup passed in and sequence tag groupings from command line
     */
    private List<FastqWriter> makeTagWriters(final SAMReadGroupRecord readGroup, final FastqWriterFactory factory) {
        String baseFilename = null;
        if (readGroup != null) {
            if (RG_TAG.equalsIgnoreCase("PU")) {
                baseFilename = readGroup.getPlatformUnit() + "_";
            } else if (RG_TAG.equalsIgnoreCase("ID")) {
                baseFilename = readGroup.getReadGroupId() + "_";
            }
            if (baseFilename == null) {
                throw new PicardException("The selected RG_TAG: " + RG_TAG + " is not present in the bam header.");
            }
        } else {
            baseFilename = "";
        }

        List<File> tagFiles = new ArrayList<>();
        for (String tagSplit : SEQUENCE_TAG_GROUP) {
            String fileName = baseFilename + tagSplit.replace(",", "_");
            fileName = IOUtil.makeFileNameSafe(fileName);
            fileName += COMPRESS_OUTPUTS_PER_TAG_GROUP ? ".fastq.gz" : ".fastq";

            final File result = (OUTPUT_DIR != null)
                    ? new File(OUTPUT_DIR, fileName)
                    : new File(FASTQ.getParent(), fileName);
            IOUtil.assertFileIsWritable(result);
            tagFiles.add(result);
        }
        return tagFiles.stream().map(factory::newWriter).collect(Collectors.toList());
    }

    // Sets up the Groupings of Sequence Tags, Quality Tags, and Separator Strings so we dont have to calculate them for every loop
    private void setupTagSplitValues() {
        SPLIT_SEQUENCE_TAGS = new ArrayList<>();
        SPLIT_QUALITY_TAGS = new ArrayList<>();
        SPLIT_SEPARATOR_TAGS = new ArrayList<>();

        for (int i = 0; i < SEQUENCE_TAG_GROUP.size(); i++) {
            SPLIT_SEQUENCE_TAGS.add(SEQUENCE_TAG_GROUP.get(i).trim().split(","));
            SPLIT_QUALITY_TAGS.add(QUALITY_TAG_GROUP.isEmpty() ? null : QUALITY_TAG_GROUP.get(i).trim().split(","));
            SPLIT_SEPARATOR_TAGS.add(TAG_GROUP_SEPERATOR.isEmpty() ? TAG_SPLIT_DEFAULT_SEP : TAG_GROUP_SEPERATOR.get(i));
        }
    }

    private void writeTagRecords(final SAMRecord read, final Integer mateNumber, final List<FastqWriter> tagWriters) {
        if (SEQUENCE_TAG_GROUP.isEmpty()) {
            return;
        }

        final String seqHeader = mateNumber == null ? read.getReadName() : read.getReadName() + "/" + mateNumber;

        for (int i = 0; i < SEQUENCE_TAG_GROUP.size(); i++) {
            final String tmpTagSep = SPLIT_SEPARATOR_TAGS.get(i);
            final String[] sequenceTagsToWrite = SPLIT_SEQUENCE_TAGS.get(i);
            final String newSequence = String.join(tmpTagSep, Arrays.stream(sequenceTagsToWrite)
                    .map(tag -> assertTagExists(read, tag))
                    .collect(Collectors.toList()));

            final String tmpQualSep = StringUtils.repeat(TAG_SPLIT_QUAL, tmpTagSep.length());
            final String[] qualityTagsToWrite = SPLIT_QUALITY_TAGS.get(i);
            final String newQual = QUALITY_TAG_GROUP.isEmpty() ? StringUtils.repeat(TAG_SPLIT_QUAL, newSequence.length()) :
                    String.join(tmpQualSep, Arrays.stream(qualityTagsToWrite)
                            .map(tag -> assertTagExists(read, tag))
                            .collect(Collectors.toList()));
            FastqWriter writer = tagWriters.get(i);
            writer.write(new FastqRecord(seqHeader, newSequence, "", newQual));
        }
    }

    private String assertTagExists(final SAMRecord record, final String tag) {
        String value = record.getStringAttribute(tag);
        if (value == null) {
            throw new PicardException("Record: " + record.getReadName() + " does have a value for tag: " + tag);
        }
        return value;
    }

    @Override
    protected String[] customCommandLineValidation() {
        List<String> errors = new ArrayList<>();

        if (!QUALITY_TAG_GROUP.isEmpty() && SEQUENCE_TAG_GROUP.size() != QUALITY_TAG_GROUP.size()) {
            errors.add("QUALITY_TAG_GROUP size must be equal to SEQUENCE_TAG_GROUP or not be specified at all.");
        }

        if (!TAG_GROUP_SEPERATOR.isEmpty() && SEQUENCE_TAG_GROUP.size() != TAG_GROUP_SEPERATOR.size()) {
            errors.add("TAG_GROUP_SEPERATOR size must be equal to SEQUENCE_TAG_GROUP or not be specified at all.");
        }

        return errors.isEmpty() ? super.customCommandLineValidation() : errors.toArray(new String[errors.size()]);
    }
}
