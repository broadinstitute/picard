package picard.arrays.illumina;

import htsjdk.samtools.Defaults;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.argumentcollections.ReferenceArgumentCollection;
import picard.vcf.ByIntervalListVariantContextIterator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Create an Extended Illumina Manifest by performing a liftover to Build 37.
 */
@CommandLineProgramProperties(
        summary = CreateExtendedIlluminaManifest.USAGE_DETAILS,
        oneLineSummary = "Create an Extended Illumina Manifest for usage by the Picard tool GtcToVcf",
        programGroup = picard.cmdline.programgroups.GenotypingArraysProgramGroup.class
)
public class CreateExtendedIlluminaManifest extends CommandLineProgram {

    static final String USAGE_DETAILS =
            "CreateExtendedIlluminaManifest takes an Illumina manifest file (this is the text version of an Illumina '.bpm' file) " +
                    "And creates an 'extended' version of this text file by adding fields that facilitate VCF generation by downstream tools. " +
                    "As part of generating this extended version of the manifest, the tool may mark loci as 'FAIL' if they do not pass validation. " +
                    "<h4>Usage example:</h4>" +
                    "<pre>" +
                    "java -jar picard.jar CreateExtendedIlluminaManifest \\<br />" +
                    "      --INPUT illumina_chip_manifest.csv \\<br />" +
                    "      --OUTPUT illumina_chip_manifest.extended.csv \\<br />" +
                    "      --REPORT_FILE illumina_chip_manifest.report.txt \\<br />" +
                    "      --CLUSTER_FILE illumina_chip_manifest.egt \\<br />" +
                    "      --REFERENCE_SEQUENCE reference.fasta \\<br />" +
                    "      --TB 37 \\<br />" +
                    "</pre>" +
                    "Some Illumina manifest files have records that are not consistently on the the build that this tool supports " +
                    "(currently Build 37).  To assist with migrating these records to Build 37, you can provide a liftover chain file " +
                    "and CreateExtendedIlluminaManifest will attempt to lift these records from the indicated build to Build 37. " +
                    "If you do not provide a liftover file, or there are records on builds other than the liftover file that you have " +
                    "provided, then those records will be marked as 'FAIL' in the extended manifest. " +
                    "<h4>Usage example with liftover:<h4>" +
                    "<pre>" +
                    "java -jar picard.jar CreateExtendedIlluminaManifest \\<br />" +
                    "      --INPUT illumina_chip_manifest.csv \\<br />" +
                    "      --OUTPUT illumina_chip_manifest.extended.csv \\<br />" +
                    "      --REPORT_FILE illumina_chip_manifest.report.txt \\<br />" +
                    "      --CLUSTER_FILE illumina_chip_manifest.egt \\<br />" +
                    "      --REFERENCE_SEQUENCE reference.fasta \\<br />" +
                    "      --TB 37 \\<br />" +
                    "      --SB 36 \\<br />" +
                    "      --SR build36_reference.fasta \\<br />" +
                    "      --SC build36ToBuild37_liftover.chain \\<br />" +
                    "</pre>" +
                    "  (that will lifover any records found on build 36 to build 37 using the build36ToBuild37 liftover file ";

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "This is the text version of the Illumina .bpm file")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The name of the extended manifest to be written.")
    public File OUTPUT;

    @Argument(shortName = "BAF", doc = "The name of the the 'bad assays file'. This is a subset version of the extended manifest, " +
            "containing only unmappable assays", optional = true)
    public File BAD_ASSAYS_FILE;

    @Argument(shortName = "RF", doc = "The name of the the report file")
    public File REPORT_FILE;

    @Argument(shortName = "FD", doc = "Flag duplicates in the extended manifest.  " +
            "If this is set and there are multiple passing assays at the same site (same locus and alleles) " +
            "then all but one will be marked with the 'DUPE' flag in the extended manifest. " +
            "The one that is not marked as 'DUPE' will be the one with the highest Gentrain score as read from the cluster file.", optional = true)
    public Boolean FLAG_DUPLICATES = true;

    @Argument(shortName = "CF", doc = "The Standard (Hapmap-trained) cluster file (.egt) from Illumina. " +
            "If there are duplicate assays at a site, this is used to decide which is the 'best' (non-filtered in generated VCFs) " +
            "by choosing the assay with the best GenTrain scores)", optional = true)
    public File CLUSTER_FILE;

    @Argument(shortName = "DBSNP", doc = "Reference dbSNP file in VCF format.", optional = true)
    public File DBSNP_FILE;

    @Argument(shortName = "TB", doc = "The target build.  This specifies the reference for which the extended manifest will be generated. " +
            "Currently this tool only supports Build 37 (Genome Reference Consortium Human Build 37 (GRCh37)). " +
            "If entries are found in the Illumina manifest that are on this build they will be used with the coordinate specified in the manifest, " +
            "If there are entries found on other builds, they will be marked as failed in the extended manifest UNLESS the " +
            "build and liftover information (SUPPORTED_BUILD, SUPPORTED_REFERENCE_FILE, and SUPPORTED_CHAIN_FILE) is supplied.")
    public String TARGET_BUILD = "37";

    @Argument(shortName = "SB", doc = "A supported build. The order of the input must match the order for SUPPORTED_REFERENCE_FILE and SUPPORTED_CHAIN_FILE. " +
            "This is the name of the build as specified in the 'GenomeBuild' column of the Illumina manifest file.",
            optional = true)
    public List<String> SUPPORTED_BUILD;

    @Argument(shortName = "SR", doc = "A reference file for the provided SUPPORTED_BUILD. " +
            "This is the reference file that corresponds to the 'SUPPORTED_BUILD' as specified above.",
            optional = true)
    public List<File> SUPPORTED_REFERENCE_FILE;

    @Argument(shortName = "SC", doc = "A chain file that maps from SUPPORTED_BUILD -> TARGET_BUILD. Must provide a corresponding supported reference file.",
            optional = true)
    public List<File> SUPPORTED_CHAIN_FILE;

    private static final Log log = Log.getInstance(CreateExtendedIlluminaManifest.class);

    private static final String VERSION = "2.0";

    @Override
    protected ReferenceArgumentCollection makeReferenceArgumentCollection() {
        return new ReferenceArgumentCollection() {
            @Argument(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, common=false,
                    doc = "The reference sequence (fasta) for the TARGET genome build.")
            public final File REFERENCE_SEQUENCE = Defaults.REFERENCE_FASTA;

            @Override
            public File getReferenceFile() {
                return REFERENCE_SEQUENCE;
            }
        };
    }

    @Override
    protected int doWork() {
        try {
            // Load the sequence dictionary from the Target Reference file
            final SAMSequenceDictionary sequenceDictionary = SAMSequenceDictionaryExtractor.extractDictionary(REFERENCE_SEQUENCE.toPath());

            ProgressLogger logger = new ProgressLogger(log, 10000);
            final Map<String, ReferenceSequenceFile> referenceSequenceMap = new HashMap<>();
            final Map<String, File> chainFilesMap = new HashMap<>();

            referenceSequenceMap.put(TARGET_BUILD, ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE));

            for (int i = 0; i < SUPPORTED_BUILD.size(); i++) {
                referenceSequenceMap.put(SUPPORTED_BUILD.get(i), ReferenceSequenceFileFactory.getReferenceSequenceFile(SUPPORTED_REFERENCE_FILE.get(i)));
                chainFilesMap.put(SUPPORTED_BUILD.get(i), SUPPORTED_CHAIN_FILE.get(i));
            }

            // Open the Original Illumina Manifest
            final IlluminaManifest manifestFile = new IlluminaManifest(INPUT);

            IOUtil.assertFileIsWritable(OUTPUT);
            if (BAD_ASSAYS_FILE != null) {
                IOUtil.assertFileIsWritable(BAD_ASSAYS_FILE);
            }
            IOUtil.assertFileIsWritable(REPORT_FILE);

            IntervalList manifestSnpIntervals = new IntervalList(sequenceDictionary);
            IntervalList manifestIndelIntervals = new IntervalList(sequenceDictionary);

            Build37ExtendedIlluminaManifestRecordCreator creator = new Build37ExtendedIlluminaManifestRecordCreator(TARGET_BUILD, referenceSequenceMap, chainFilesMap);

            // first iteration through the manifest to find all dupes
            log.info("Phase 1.  First Pass through the manifest.  Build coordinate map for dupe flagging and make SNP and indel-specific interval lists for parsing dbSnp");
            final Iterator<IlluminaManifestRecord> firstPassIterator = manifestFile.iterator();

            List<Build37ExtendedIlluminaManifestRecord> records = new ArrayList<>();

            while (firstPassIterator.hasNext()) {
                logger.record("0", 0);
                final IlluminaManifestRecord record = firstPassIterator.next();

                // Create an ExtendedIlluminaManifestRecord here so that we can get the (potentially lifted over) coordinates
                final Build37ExtendedIlluminaManifestRecord rec = creator.createRecord(record);
                records.add(rec);

                if (!rec.isFail()) {
                    final int length = Integer.max(rec.getAlleleA().length(), rec.getAlleleB().length());
                    Interval interval = new Interval(rec.getB37Chr(), rec.getB37Pos(), rec.getB37Pos() + length);
                    if (rec.isSnp()) {
                        manifestSnpIntervals.add(interval);
                    } else {
                        manifestIndelIntervals.add(interval);
                    }
                }
            }

            // Generate a sorted set of the variants in the Illumina Manifest so that we can check them
            // Against the (sorted) dbSnp Vcf.
            manifestSnpIntervals = manifestSnpIntervals.sorted();
            manifestIndelIntervals = manifestIndelIntervals.sorted();

            log.info("Phase 2.  Parse dbSnpVCF and build SNP and indel-specific locus to rsId maps");
            Map<String, String> snpLocusToRsId = new HashMap<>();
            Map<String, String> indelLocusToRsId = new HashMap<>();
            if (DBSNP_FILE != null) {
                // Because dbSnp can contain both SNPs and indels which may be at the same locus,
                // We do two passes through dbSnpVcf to build separate maps.
                log.info("SNP-specific");
                snpLocusToRsId = generateLocusToRsidMap(DBSNP_FILE, manifestSnpIntervals);
                log.info("indel-specific");
                indelLocusToRsId = generateLocusToRsidMap(DBSNP_FILE, manifestIndelIntervals);
            }

            Set<Integer> dupeIndices = null;
            if (FLAG_DUPLICATES) {
                dupeIndices = flagDuplicates(records);
            }

            final BufferedWriter out = new BufferedWriter(new FileWriter(OUTPUT, false));
            writeExtendedIlluminaManifestHeaders(manifestFile, out);

            // second iteration to write all records after dupe evaluation
            log.info("Phase 3.  Generate the Extended Illumina Manifest");
            logger = new ProgressLogger(log, 10000);
            ManifestStatistics manifestStatistics = new ManifestStatistics(TARGET_BUILD);

            List<Build37ExtendedIlluminaManifestRecord> badRecords = new ArrayList<>();
            for (Build37ExtendedIlluminaManifestRecord record: records) {
                logger.record("0", 0);
                final String locus = record.getChr() + "." + record.getPosition();
                String rsId;
                if (record.isSnp()) {
                    rsId = snpLocusToRsId.get(locus);
                } else {
                    rsId = indelLocusToRsId.get(locus);
                }
                record.setRsId(rsId);
                if (record.isFail()) {
                    badRecords.add(record);
                } else {
                    if (dupeIndices != null) {
                        record.setDupe(dupeIndices.contains(record.getIndex()));
                    }
                }
                manifestStatistics.updateStatistics(record);
                out.write(record.getLine());
                out.newLine();
            }

            out.flush();
            out.close();

            if (BAD_ASSAYS_FILE != null) {
                writeBadAssaysFile(BAD_ASSAYS_FILE, badRecords);
            }
            StringBuilder sb = new StringBuilder();
            sb.append("CreateExtendedIlluminaManifest (version: ").append(VERSION).append(") Report For: ").append(OUTPUT.getName()).append("\n");
            sb.append("Generated on: ").append(new Date()).append("\n");
            sb.append("Using Illumina Manifest: ").append(INPUT.getAbsolutePath()).append("\n");
            if (FLAG_DUPLICATES) {
                sb.append("Duplicates were flagged\n");
            }
            if (CLUSTER_FILE != null) {
                sb.append("Using Illumina EGT: ").append(CLUSTER_FILE.getAbsolutePath()).append("\n");
            }
            sb.append("\n");
            if (!creator.isRefStrandDefinedInManifest() || (!creator.getUnsupportedBuilds().isEmpty())) {
                sb.append("NOTES / Warnings:\n");
                if (!creator.isRefStrandDefinedInManifest()) {
                    sb.append("REF_STRAND was NOT defined in the manifest.  We have inferred it from sequence / strand information in the manifest.\n");
                }
                if (!creator.getUnsupportedBuilds().isEmpty()) {
                    sb.append("Records were found within the manifest on Genome Builds for which you have not provided liftover information.\n");
                    sb.append(" They have been failed with the flag: UNSUPPORTED_GENOME_BUILD\n");
                    sb.append(" The following unexpected Genome Builds were found: ").append(StringUtils.join(creator.getUnsupportedBuilds(), ", ")).append("\n");
                }
                sb.append("\n");
            }
            manifestStatistics.logStatistics(REPORT_FILE, sb.toString());
        } catch (IOException e) {
            throw new PicardException(e.getMessage(), e);
        }

        return 0;
    }

    /**
     * This method goes through a list of records and if it finds duplicates, it flags the non-best
     * with the duplicate flag.
     *
     * It iterates over all of the records, looking for records that share the same locus (chrom/position) and alleles.
     * If it finds more than one record that share these attributes, it will flag all but one as 'duplicate'
     *  this flag will later be used to set the 'DUPE' filter in the GtcToVcf tool.
     * The decision as to which record will be flagged as duplicate is based upon the GenTrain score
     *  (cluster quality from the GenTrain clustering algorithm) which is pulled from the Illumina cluster file.
     * The record with the highest GenTrain score will NOT be flagged as a duplicate, all others will.
     *
     * @param records A list of records to search through and flag for duplicates.
     * @return a set of integers indicating the indices of records that are to be flagged as duplicates.
     */
    private Set<Integer> flagDuplicates(final List<Build37ExtendedIlluminaManifestRecord> records) {
        // Load the cluster file to get the GenTrain scores
        // load the egt first, and create a map of ilmnid to gentrain score.  Save that and use it for deduplicating.
        log.info("Loading the egt file for duplicate resolution");
        final InfiniumEGTFile infiniumEGTFile;
        try {
            infiniumEGTFile = new InfiniumEGTFile(CLUSTER_FILE);
        } catch (IOException e) {
            throw new PicardException("Error reading cluster file '" + CLUSTER_FILE.getAbsolutePath() + "'", e);
        }
        final Map<String, Float> nameToGenTrainScore = new HashMap<>();
        for (String rsName : infiniumEGTFile.rsNameToIndex.keySet()) {
            nameToGenTrainScore.put(rsName, infiniumEGTFile.totalScore[infiniumEGTFile.rsNameToIndex.get(rsName)]);
        }

        return flagDuplicates(records, nameToGenTrainScore);
    }

    /**
     * This method goes through a list of records and if it finds duplicates, it flags the non-best
     * with the duplicate flag.
     *
     * It iterates over all of the records, looking for records that share the same locus (chrom/position) and alleles.
     * If it finds more than one record that share these attributes, it will flag all but one as 'duplicate'
     *  this flag will later be used to set the 'DUPE' filter in the GtcToVcf tool.
     * The decision as to which record will be flagged as duplicate is based upon the GenTrain score
     *  (cluster quality from the GenTrain clustering algorithm) which is pulled from the Illumina cluster file.
     * The record with the highest GenTrain score will NOT be flagged as a duplicate, all others will.
     *
     * This version was separated out to facilitate testing (without having to load an entire cluster file)
     *
     * @param records A list of records to search through and flag for duplicates.
     * @param nameToGenTrainScore A map of name to gentrain score.
     * @return a set of integers indicating the indices of records that are to be flagged as duplicates.
     */
    Set<Integer> flagDuplicates(final List<Build37ExtendedIlluminaManifestRecord> records, final Map<String, Float> nameToGenTrainScore) {
        Map<String, List<Build37ExtendedIlluminaManifestRecord>> coordinateMap = new HashMap<>();
        for (Build37ExtendedIlluminaManifestRecord record : records) {

            String key = record.getB37Chr() + ":" + record.getB37Pos() + "." + record.getSnpRefAllele();
            if (!record.getSnpAlleleA().equals(record.getSnpRefAllele())) {
                key += "." + record.getSnpAlleleA();
            }
            if (!record.getSnpAlleleB().equals(record.getSnpAlleleA()) && !(record.getSnpAlleleB().equals(record.getSnpRefAllele()))) {
                key += "." + record.getSnpAlleleB();
            }

            if (!record.isFail()) {
                if (coordinateMap.containsKey(key)) {
                    coordinateMap.get(key).add(record);
                } else {
                    List<Build37ExtendedIlluminaManifestRecord> newList = new ArrayList<>();
                    newList.add(record);
                    coordinateMap.put(key, newList);
                }
            }
        }

        // filter out all unique coordinates
        Map<String, List<Build37ExtendedIlluminaManifestRecord>> dupeMap = coordinateMap.entrySet().stream()
                .filter(map -> map.getValue().size() > 1)
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
        coordinateMap.clear();

        // evaluate each coordinate assay and remove the assay with the best GenTrain score (all remaining are dupes)
        dupeMap.entrySet().forEach(entry ->
                entry.getValue().remove(entry.getValue().stream().max(Comparator.comparingDouble(assay ->
                        nameToGenTrainScore.get(assay.getName()))).get()));

        // we really only need the list of indices for the dupes
        return dupeMap.entrySet().stream()
                .flatMapToInt(entry ->
                        entry.getValue().stream()
                                .mapToInt(IlluminaManifestRecord::getIndex))
                .boxed().collect(Collectors.toSet());
    }

    private void writeBadAssaysFile(File badAssaysFile, List<Build37ExtendedIlluminaManifestRecord> badRecords) throws IOException {
        BufferedWriter badAssaysFileWriter = new BufferedWriter(new FileWriter(badAssaysFile, false));
        badAssaysFileWriter.write("## The following assays were marked by CreateExtendedIlluminaManifest as Unparseable (input file: " + INPUT.getAbsolutePath() + ")");
        badAssaysFileWriter.newLine();
        badAssaysFileWriter.write("#IlmnId,Name,GenomeBuild,Chr,MapInfo,FailureFlag");
        badAssaysFileWriter.newLine();

        for (Build37ExtendedIlluminaManifestRecord record : badRecords) {
            final List<String> badRecord = java.util.Arrays.asList(record.getIlmnId(), record.getName(), record.getGenomeBuild(), record.getChr(), "" + record.getPosition(), record.getFlag().toString());
            badAssaysFileWriter.write(StringUtils.join(badRecord, ","));
            badAssaysFileWriter.newLine();
        }
        badAssaysFileWriter.flush();
        badAssaysFileWriter.close();
    }

    private static class ManifestStatistics {
        private final String targetBuild;

        int numAssays;
        int numAssaysFlagged;
        int numAssaysDuplicated;        // The number of passing assays which are flagged as duplicates

        int numSnps;
        int numSnpsDuplicated;          // The number of passing SNP assays which are flagged as duplicates
        int numSnpsFlagged;
        int numSnpsIlluminaFlagged;
        int numSnpProbeSequenceMismatch;
        int numSnpMissingAlleleBProbeSequence;

        int numAmbiguousSnpsOnPosStrand;
        int numAmbiguousSnpsOnNegStrand;

        int numIndels;
        int numIndelsDuplicated;        // The number of passing SNP assays which are flagged as duplicates
        int numIndelsFlagged;
        int numIndelsIlluminaFlagged;
        int numIndelProbeSequenceMismatch;
        int numIndelSourceSequenceInvalid;
        int numIndelsNotFound;
        int numIndelConfict;

        int numOnTargetBuild;
        Map<String, Integer> numOnOtherBuild;
        int numOnUnsupportedGenomeBuild;
        int numLiftoverFailed;
        int numRefStrandMismatch;

        public ManifestStatistics(final String targetBuild) {
            this.targetBuild = targetBuild;
            this.numOnOtherBuild = new TreeMap<>();
        }

        void updateStatistics(Build37ExtendedIlluminaManifestRecord rec) {
            numAssays++;
            if (rec.isSnp()) {
                numSnps++;
            } else {
                numIndels++;
            }
            if (rec.getMajorGenomeBuild().equals(targetBuild)) {
                numOnTargetBuild++;
            } else {
                Integer num = numOnOtherBuild.get(rec.getMajorGenomeBuild());
                if (num == null) {
                    num = 0;
                }
                numOnOtherBuild.put(rec.getMajorGenomeBuild(), ++num);
            }
            if (rec.getFlag().equals(Build37ExtendedIlluminaManifestRecord.Flag.UNSUPPORTED_GENOME_BUILD)) {
                numOnUnsupportedGenomeBuild++;
            }
            if (rec.getFlag().equals(Build37ExtendedIlluminaManifestRecord.Flag.LIFTOVER_FAILED)) {
                numLiftoverFailed++;
            }
            if (!rec.isFail()) {
                if (rec.isDupe()) {
                    numAssaysDuplicated++;
                    if (rec.isSnp()) {
                        numSnpsDuplicated++;
                    } else {
                        numIndelsDuplicated++;
                    }
                }
                if (rec.isAmbiguous()) {
                    if (rec.getRefStrand() == Strand.NEGATIVE) {
                        numAmbiguousSnpsOnNegStrand++;
                    }
                    if (rec.getRefStrand() == Strand.POSITIVE) {
                        numAmbiguousSnpsOnPosStrand++;
                    }
                }
            } else {
                numAssaysFlagged++;
                if (rec.isSnp()) {
                    numSnpsFlagged++;
                    switch (rec.getFlag()) {
                        case ILLUMINA_FLAGGED:
                            numSnpsIlluminaFlagged++;
                            break;
                        case PROBE_SEQUENCE_MISMATCH:
                            numSnpProbeSequenceMismatch++;
                            break;
                        case MISSING_ALLELE_B_PROBESEQ:
                            numSnpMissingAlleleBProbeSequence++;
                            break;
                        case UNSUPPORTED_GENOME_BUILD:
                        case LIFTOVER_FAILED:
                            break;          // These are covered above
                        default:
                            throw new PicardException("Unhandled Flag: " + rec.getFlag());
                    }
                } else {
                    numIndelsFlagged++;
                    switch (rec.getFlag()) {
                        case ILLUMINA_FLAGGED:
                            numIndelsIlluminaFlagged++;
                            break;
                        case PROBE_SEQUENCE_MISMATCH:
                            numIndelProbeSequenceMismatch++;
                            break;
                        case SOURCE_SEQUENCE_INVALID:
                            numIndelSourceSequenceInvalid++;
                            break;
                        case INDEL_NOT_FOUND:
                            numIndelsNotFound++;
                            break;
                        case INDEL_CONFLICT:
                            numIndelConfict++;
                            break;
                        case UNSUPPORTED_GENOME_BUILD:
                        case LIFTOVER_FAILED:
                            break;          // These are covered above
                        default:
                            throw new PicardException("Unhandled Flag: " + rec.getFlag());
                    }
                }
            }
        }

        void logStatistics(File output, final String header) throws IOException {
            try (BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(
                    new FileOutputStream(output), StandardCharsets.UTF_8))) {
                writer.write(header);

                writer.write("Total Number of Assays: " + numAssays);
                writer.newLine();

                writer.write("Number of Assays on Build " + targetBuild + ": " + numOnTargetBuild);
                writer.newLine();
                int numOnOtherBuilds = 0;
                for (final String build : numOnOtherBuild.keySet()) {
                    writer.write("Number of Assays on Build " + build + ": " + numOnOtherBuild.get(build));
                    writer.newLine();
                    numOnOtherBuilds += numOnOtherBuild.get(build);
                }
                writer.write("Number of Assays on unsupported genome build: " + numOnUnsupportedGenomeBuild);
                writer.newLine();
                writer.write("Number of Assays failing liftover: " + numLiftoverFailed);
                writer.newLine();
                writer.newLine();

                writer.write("Number of Assays on Build " + targetBuild + " or successfully lifted over: " + (numOnTargetBuild + (numOnOtherBuilds - numOnUnsupportedGenomeBuild - numLiftoverFailed)));
                writer.newLine();
                writer.write("Number of Passing Assays: " + (numAssays - numAssaysFlagged));
                writer.newLine();
                writer.write("Number of Duplicated Assays: " + numAssaysDuplicated);
                writer.newLine();
                writer.write("Number of Failing Assays: " + numAssaysFlagged);
                writer.newLine();
                writer.newLine();
                writer.write("Number of SNPs: " + numSnps);
                writer.newLine();
                writer.write("Number of Passing SNPs: " + (numSnps - numSnpsFlagged));
                writer.newLine();
                writer.write("Number of Duplicated SNPs: " + numSnpsDuplicated);
                writer.newLine();
                writer.newLine();
                writer.write("Number of Failing SNPs: " + numSnpsFlagged);
                writer.newLine();
                writer.write("Number of SNPs failed by Illumina: " + numSnpsIlluminaFlagged);
                writer.newLine();
                // Note - currently this is only calculated for SNPs - that's why it's not in the indel section too.
                writer.write("Number of SNPs failed for refStrand mismatch: "  + numRefStrandMismatch);
                writer.newLine();
                writer.write("Number of SNPs failed for missing AlleleB ProbeSeq: "  + numSnpMissingAlleleBProbeSequence);
                writer.newLine();
                writer.write("Number of SNPs failed for alleleA probe sequence mismatch: " + numSnpProbeSequenceMismatch);
                writer.newLine();
                writer.write("Number of ambiguous SNPs on Positive Strand: " + numAmbiguousSnpsOnPosStrand);
                writer.newLine();
                writer.write("Number of ambiguous SNPs on Negative Strand: " + numAmbiguousSnpsOnNegStrand);
                writer.newLine();
                writer.newLine();

                writer.write("Number of Indels: " + numIndels);
                writer.newLine();
                writer.write("Number of Passing Indels: " + (numIndels - numIndelsFlagged));
                writer.newLine();
                writer.write("Number of Duplicated Indels: " + numIndelsDuplicated);
                writer.newLine();
                writer.newLine();
                writer.write("Number of Failing Indels: " + numIndelsFlagged);
                writer.newLine();
                writer.write("Number of Indels failed by Illumina: " + numIndelsIlluminaFlagged);
                writer.newLine();
                writer.write("Number of Indels failed for probe sequence mismatch: " + numIndelProbeSequenceMismatch);
                writer.newLine();
                writer.write("Number of Indels failed for source sequence invalid: " + numIndelSourceSequenceInvalid);
                writer.newLine();
                writer.write("Number of Indels not found: " + numIndelsNotFound);
                writer.newLine();
                writer.write("Number of Indels failed for conflict: " + numIndelConfict);
                writer.newLine();
            }
        }
    }

    @Override
    protected String[] customCommandLineValidation() {
        IOUtil.assertFileIsReadable(INPUT);
        if (CLUSTER_FILE != null) {
            IOUtil.assertFileIsReadable(CLUSTER_FILE);
        }
        if (DBSNP_FILE != null) {
            IOUtil.assertFileIsReadable(DBSNP_FILE);
        }

        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
        for (File f : SUPPORTED_REFERENCE_FILE) IOUtil.assertFileIsReadable(f);
        for (File f : SUPPORTED_CHAIN_FILE) IOUtil.assertFileIsReadable(f);

        final List<String> errors = new ArrayList<>();

        if (!TARGET_BUILD.equals(Build37ExtendedIlluminaManifestRecordCreator.BUILD_37)) {
            errors.add("Currently this tool only supports Build 37");
        }

        if (FLAG_DUPLICATES && CLUSTER_FILE == null) {
            errors.add("In order to flag duplicates, a CLUSTER_FILE must be supplied");
        }

        if ((!SUPPORTED_BUILD.isEmpty()) || (!SUPPORTED_REFERENCE_FILE.isEmpty()) || (!SUPPORTED_CHAIN_FILE.isEmpty())) {
            if ((SUPPORTED_BUILD.isEmpty()) || (SUPPORTED_REFERENCE_FILE.isEmpty()) || (SUPPORTED_CHAIN_FILE.isEmpty())) {
                errors.add("Parameters for 'SUPPORTED_BUILD', 'SUPPORTED_REFERENCE_FILE', and 'SUPPORTED_CHAIN_FILE' must ALL be specified or not at all.");
            } else {
                if (SUPPORTED_BUILD.size() != SUPPORTED_REFERENCE_FILE.size()) {
                    errors.add("The number of inputs for 'SUPPORTED_BUILD' does not match the number of inputs for 'SUPPORTED_REFERENCE_FILE'");
                }

                if (SUPPORTED_BUILD.size() != SUPPORTED_CHAIN_FILE.size()) {
                    errors.add("The number of inputs for 'SUPPORTED_BUILD' does not match the number of inputs for 'SUPPORTED_CHAIN_FILE'");
                }
            }
        }

        return (errors.size() > 0)
                ? errors.toArray(new String[errors.size()])
                : null;
    }

    /**
     * Generates a mapping of locus (contig.posn) in the manifest file to rsId.
     * Uses the passed interval list to selectively parse dbSnpVcf.
     * Returns a map of locus to rsId
     * @param dbSnpFile the dbSnp file to parse.
     * @param intervals interval list for which intervals to parse out of the dbSnp file
     * @return mapping of locus in the manifest file (contig.posn) to rsId
     */
    Map<String, String> generateLocusToRsidMap(File dbSnpFile, IntervalList intervals) {
        ProgressLogger logger = new ProgressLogger(log, 10000);

        Map<String, String> manifestLocusToRsId = new HashMap<>();

        final VCFFileReader dbSnpReader = new VCFFileReader(dbSnpFile, true);
        final Iterator<VariantContext> dbSnpIterator = new ByIntervalListVariantContextIterator(dbSnpReader, intervals);
        while (dbSnpIterator.hasNext()) {
            VariantContext variantContext = dbSnpIterator.next();
            logger.record(variantContext.getContig(), variantContext.getStart());

            for (int posn = variantContext.getStart(); posn <= variantContext.getEnd(); posn++) {
                final String locus = variantContext.getContig() + "." + posn;
                manifestLocusToRsId.put(locus, variantContext.getID());
            }
        }

        return manifestLocusToRsId;
    }

    void writeExtendedIlluminaManifestHeaders(final IlluminaManifest manifest, final BufferedWriter output) throws IOException {
        int numColumns = -1;
        List<String[]> currentHeader = manifest.getHeaderContents();
        String[] lastRowInHeader = currentHeader.get(currentHeader.size() - 1); // "Loci Count" which needs to be last to terminate the header...
        for (int i = 0; i < currentHeader.size() - 1; i++) {
            String[] rowValues = currentHeader.get(i);
            if (numColumns == -1) {
                numColumns = rowValues.length;
            }
            addHeaderLine(output, numColumns, rowValues);
        }
        addHeaderLine(output, numColumns, Build37ExtendedIlluminaManifest.EXTENDED_MANIFEST_VERSION_HEADER_NAME, VERSION);
        addHeaderLine(output, numColumns, Build37ExtendedIlluminaManifest.EXTENDED_MANIFEST_TARGET_BUILD_HEADER_NAME, TARGET_BUILD);
        addHeaderLine(output, numColumns, Build37ExtendedIlluminaManifest.EXTENDED_MANIFEST_TARGET_REFERENCE_HEADER_NAME, REFERENCE_SEQUENCE.getAbsolutePath());
        if (CLUSTER_FILE != null) {
            addHeaderLine(output, numColumns, Build37ExtendedIlluminaManifest.EXTENDED_MANIFEST_CLUSTER_FILE_HEADER_NAME, CLUSTER_FILE.getAbsolutePath());
        }
        if (DBSNP_FILE != null) {
            addHeaderLine(output, numColumns, Build37ExtendedIlluminaManifest.EXTENDED_MANIFEST_DBSNP_FILE_HEADER_NAME, DBSNP_FILE.getAbsolutePath());
        }

        if (!SUPPORTED_BUILD.isEmpty()) {
            final String[] supportedBuildsFields = new String[SUPPORTED_BUILD.size() + 1];
            final String[] supportedReferenceFileFields = new String[SUPPORTED_BUILD.size() + 1];
            final String[] supportedChainFileFields = new String[SUPPORTED_BUILD.size() + 1];
            supportedBuildsFields[0] = Build37ExtendedIlluminaManifest.EXTENDED_MANIFEST_SUPPORTED_BUILD_HEADER_NAME;
            supportedReferenceFileFields[0] = Build37ExtendedIlluminaManifest.EXTENDED_MANIFEST_SUPPORTED_REFERENCE_HEADER_NAME;
            supportedChainFileFields[0] = Build37ExtendedIlluminaManifest.EXTENDED_MANIFEST_SUPPORTED_CHAIN_FILE_HEADER_NAME;
            for (int i = 0; i < SUPPORTED_BUILD.size(); i++) {
                supportedBuildsFields[i + 1] = SUPPORTED_BUILD.get(i);
                supportedReferenceFileFields[i + 1] = SUPPORTED_REFERENCE_FILE.get(i).getAbsolutePath();
                supportedChainFileFields[i + 1] = SUPPORTED_CHAIN_FILE.get(i).getAbsolutePath();
            }
            addHeaderLine(output, numColumns, supportedBuildsFields);
            addHeaderLine(output, numColumns, supportedReferenceFileFields);
            addHeaderLine(output, numColumns, supportedChainFileFields);
        }

        addHeaderLine(output, numColumns, lastRowInHeader);

        addHeaderLine(output, numColumns, "[Assay]");

        // write the extended headers
        final String[] extendedHeader = ArrayUtils.addAll(manifest.getManifestFileHeaderNames(), Build37ExtendedIlluminaManifest.EXTENDED_MANIFEST_HEADERS);
        output.write(StringUtils.join(extendedHeader, ","));
        output.newLine();
    }

    private void addHeaderLine(final BufferedWriter out, final int numColumns, final String... fields) throws IOException {
        String[] rowValues = new String[numColumns];
        System.arraycopy(fields, 0, rowValues, 0, fields.length);
        out.write(StringUtils.join(rowValues, ","));
        out.newLine();
    }
}


