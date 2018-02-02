package picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = SamToFastqWithTags.USAGE_SUMMARY + SamToFastqWithTags.USAGE_DETAILS,
        oneLineSummary = SamToFastqWithTags.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
public class SamToFastqWithTags extends SamToFastq {
    static final String USAGE_SUMMARY = "Converts a SAM or BAM file to FASTQ.";
    static final String USAGE_DETAILS =
            "This will create the fastq from read sequences alongside two fastq files.  One will be converted from the " +
            " \"CR\" tag with quality from the \"CY\" tag.  The other fastq will be converted from the \"CB\" and \"UR\" " +
            "tags concatenated together with no separator (not specified on command line) with the qualities coming from " +
            "the \"CY\" and \"UY\" tags concatenated together." +
            "<br />" +
            "<pre>" +
            "java -jar picard.jar SamToFastq <br />" +
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

    @Argument(shortName = "SEP", doc = "List of sequences to put in between each comma separated list of sequence tags in each SEQUENCE_TAG_GROUP (STG)", optional = true)
    public List<String> TAG_GROUP_SEPERATOR;

    @Argument(shortName = "GZOPTG", doc = "Compress output FASTQ files per Tag grouping using gzip and append a .gz extension to the file names.")
    public Boolean COMPRESS_OUTPUTS_PER_TAG_GROUP = false;

    private final Log log = Log.getInstance(SamToFastqWithTags.class);

    private final static String TAG_SPLIT_DEFAULT_SEP = "";
    private final static String TAG_SPLIT_DEFAULT_QUAL = "~";

    private ArrayList<String[]> SPLIT_SEQUENCE_TAGS;
    private ArrayList<String[]> SPLIT_QUALITY_TAGS;
    private ArrayList<String> SPLIT_SEPARATOR_TAGS;

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        setupTagSplitValues();
        final FastqWriterFactory factory = new FastqWriterFactory();
        factory.setCreateMd5(CREATE_MD5_FILE);
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        final Map<String, SAMRecord> firstSeenMates = new HashMap<>();

        final Map<SAMReadGroupRecord, FastqWriters> writers = super.generateWriters(reader.getFileHeader().getReadGroups(), factory);
        final Map<SAMReadGroupRecord, List<FastqWriter>> tagWriters = generateTagWriters(reader.getFileHeader().getReadGroups(), factory);

        if (writers.isEmpty()) {
            final String msgBase = INPUT + " does not contain Read Groups";
            final String msg = OUTPUT_PER_RG ? msgBase + ", consider not using the OUTPUT_PER_RG option" : msgBase;
            throw new PicardException(msg);
        }

        final ProgressLogger progress = new ProgressLogger(log);

        for (final SAMRecord currentRecord : reader) {
            super.handleRecord(currentRecord, writers, firstSeenMates);
            handleTagRecord(currentRecord, tagWriters, firstSeenMates);
            progress.record(currentRecord);
        }


        CloserUtil.close(reader);
        for (final FastqWriters writerMapping : new HashSet<>(writers.values())) {
            writerMapping.closeAll();

        }

        for (final List<FastqWriter> writersToClose : new HashSet<>(tagWriters.values())){
            for (FastqWriter writer : writersToClose){
                writer.close();
            }
        }

        if (!firstSeenMates.isEmpty()) {
            SAMUtils.processValidationError(new SAMValidationError(SAMValidationError.Type.MATE_NOT_FOUND,
                    "Found " + firstSeenMates.size() + " unpaired mates", null), VALIDATION_STRINGENCY);
        }

        return 0;
    }


    private Map<SAMReadGroupRecord, List<FastqWriter>> generateTagWriters(final List<SAMReadGroupRecord> samReadGroupRecords,
                                                                  final FastqWriterFactory factory) {

        final Map<SAMReadGroupRecord, List<FastqWriter>> writerMap = new HashMap<>();

        final FastqWriters fastqWriters;
        if (!OUTPUT_PER_RG) {
            /* Prepare tag writers if tag groupings are provided to the tool */

            final List<FastqWriter> tagFastqWriters = makeTagWriters(null).stream().map(factory::newWriter).collect(Collectors.toList());

            writerMap.put(null, tagFastqWriters);
            for (final SAMReadGroupRecord rg : samReadGroupRecords) {
                writerMap.put(rg, tagFastqWriters);
            }
        } else {
            // When we're creating a fastq-group per readgroup, by convention we do not emit a special fastq for unpaired reads.
            for (final SAMReadGroupRecord rg : samReadGroupRecords) {
                List<FastqWriter> tagWriters = null;

                /* Prepare tag writers if tag groupings are provided to the tool */
                tagWriters = makeTagWriters(rg).stream().map(factory::newWriter).collect(Collectors.toList());

                writerMap.put(rg, tagWriters);
            }
        }
        return writerMap;
    }

    private void handleTagRecord(final SAMRecord currentRecord,final Map<SAMReadGroupRecord, List<FastqWriter>> tagWriters, final Map<String, SAMRecord> firstSeenMates ){
        if (currentRecord.isSecondaryOrSupplementary() && !INCLUDE_NON_PRIMARY_ALIGNMENTS)
            return;

        // Skip non-PF reads as necessary
        if (currentRecord.getReadFailsVendorQualityCheckFlag() && !INCLUDE_NON_PF_READS)
            return;

        final List<FastqWriter> rgTagWriters = tagWriters.get(currentRecord.getReadGroup());
        if (currentRecord.getReadPairedFlag()) {
            final String currentReadName = currentRecord.getReadName();
            final SAMRecord firstRecord = firstSeenMates.remove(currentReadName);
            if (firstRecord == null) {
                firstSeenMates.put(currentReadName, currentRecord);
            } else {
                super.assertPairedMates(firstRecord, currentRecord);

                final SAMRecord read1 =
                        currentRecord.getFirstOfPairFlag() ? currentRecord : firstRecord;
                final SAMRecord read2 =
                        currentRecord.getFirstOfPairFlag() ? firstRecord : currentRecord;
                writeTagRecords(read1, 1, rgTagWriters);
                writeTagRecords(read2, 2, rgTagWriters);
            }
        } else {
            writeTagRecords(currentRecord, null, rgTagWriters);
        }
    }

    private List<File> makeTagWriters(final SAMReadGroupRecord readGroup) {
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
            String fileName = baseFilename;

            fileName += tagSplit.replace(",", "_");
            fileName = IOUtil.makeFileNameSafe(fileName);

            fileName += COMPRESS_OUTPUTS_PER_TAG_GROUP ? ".fastq.gz" : ".fastq";

            final File result = (OUTPUT_DIR != null)
                    ? new File(OUTPUT_DIR, fileName)
                    : new File(fileName);
            IOUtil.assertFileIsWritable(result);
            tagFiles.add(result);
        }
        return tagFiles;
    }

    // Setting up the Groupings of Sequence Tags, Quality Tags, and Separator Strings so we dont have to calculate them for every loop
    private void setupTagSplitValues() {
        if (SEQUENCE_TAG_GROUP.isEmpty()) return;

        SPLIT_SEQUENCE_TAGS = new ArrayList<>();
        SPLIT_QUALITY_TAGS = new ArrayList<>();
        SPLIT_SEPARATOR_TAGS = new ArrayList<>();

        for (int i = 0; i < SEQUENCE_TAG_GROUP.size(); i ++){
            SPLIT_SEQUENCE_TAGS.add(SEQUENCE_TAG_GROUP.get(i).trim().split(","));
            SPLIT_QUALITY_TAGS.add(QUALITY_TAG_GROUP.isEmpty() ? null : QUALITY_TAG_GROUP.get(i).trim().split(","));
            SPLIT_SEPARATOR_TAGS.add(TAG_GROUP_SEPERATOR.isEmpty() ? TAG_SPLIT_DEFAULT_SEP : TAG_GROUP_SEPERATOR.get(i));
        }
    }

    private void writeTagRecords (final SAMRecord read, final Integer mateNumber, final List<FastqWriter> tagWriters){
        if (SEQUENCE_TAG_GROUP.isEmpty()) return;

        final String seqHeader = mateNumber == null ? read.getReadName() : read.getReadName() + "/" + mateNumber;

        for (int i = 0; i < SEQUENCE_TAG_GROUP.size(); i ++){
            final String tmpTagSep = SPLIT_SEPARATOR_TAGS.get(i);
            final String[] sequenceTagsToWrite = SPLIT_SEQUENCE_TAGS.get(i);
            final String newSequence = String.join(tmpTagSep, Arrays.stream(sequenceTagsToWrite)
                    .map(tag -> assertTagExists(read, tag))
                    .collect(Collectors.toList()));

            final String tmpQualSep = StringUtils.repeat(TAG_SPLIT_DEFAULT_QUAL, tmpTagSep.length());
            final String[] qualityTagsToWrite = SPLIT_QUALITY_TAGS.get(i);
            final String newQual = QUALITY_TAG_GROUP.isEmpty() ? StringUtils.repeat(TAG_SPLIT_DEFAULT_QUAL, newSequence.length()):
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
            throw new PicardException("Record: " + record.getReadName() + " does have a value for tag: " + tag );
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


        if (!errors.isEmpty()) return errors.toArray(new String[errors.size()]);

        return super.customCommandLineValidation();
    }
}
