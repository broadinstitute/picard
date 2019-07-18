/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
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

package picard.arrays;

import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.arrays.illumina.ArraysControlInfo;
import picard.arrays.illumina.Build37ExtendedIlluminaManifest;
import picard.arrays.illumina.Build37ExtendedIlluminaManifestRecord;
import picard.arrays.illumina.IlluminaManifestRecord;
import picard.arrays.illumina.InfiniumEGTFile;
import picard.arrays.illumina.InfiniumGTCFile;
import picard.arrays.illumina.InfiniumGTCRecord;
import picard.arrays.illumina.InfiniumNormalizationManifest;
import picard.arrays.illumina.InfiniumVcfFields;
import picard.pedigree.Sex;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFRecordCodec;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Class to convert a GTC file and a BPM file to a VCF file.
 */
@CommandLineProgramProperties(
        summary = GtcToVcf.USAGE_DETAILS,
        oneLineSummary = "Program to convert a GTC file to a VCF",
        programGroup = picard.cmdline.programgroups.GenotypingArraysProgramGroup.class
)
@DocumentedFeature
public class GtcToVcf extends CommandLineProgram {


    static final String USAGE_DETAILS =
            "GtcToVcf takes an Illumina GTC file and converts it to a VCF file using several supporting files. " +
                    "A GTC file is an Illumina-specific file containing called genotypes in AA/AB/BB format. " +
                    "<a href='https://github.com/Illumina/BeadArrayFiles/blob/develop/docs/GTC_File_Format_v5.pdf'></a> " +
                    "A VCF, aka Variant Calling Format, is a text file for storing how a sequenced sample differs from the reference genome. " +
                    "<a href='http://software.broadinstitute.org/software/igv/book/export/html/184'></a>" +
                    "<h4>Usage example:</h4>" +
                    "<pre>" +
                    "java -jar picard.jar GtcToVcf \\<br />" +
                    "      INPUT=input.gtc \\<br />" +
                    "      REFERENCE_SEQUENCE=reference.fasta \\<br />" +
                    "      OUTPUT=output.vcf \\<br />" +
                    "      EXTENDED_ILLUMINA_MANIFEST=chip_name.extended.csv \\<br />" +
                    "      CLUSTER_FILE=chip_name.egt \\<br />" +
                    "      ILLUMINA_NORMALIZATION_MANIFEST=chip_name.bpm.csv \\<br />" +
                    "      SAMPLE_ALIAS=my_sample_alias \\<br />" +
                    "</pre>";

    private final static Log log = Log.getInstance(GtcToVcf.class);

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "GTC file to be converted")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output VCF file to write.")
    public File OUTPUT;

    @Argument(shortName = "MANIFEST", doc = "An Extended Illumina Manifest file (csv).  This is an extended version of the Illumina manifest" +
            " it contains additional reference-specific fields")
    public File EXTENDED_ILLUMINA_MANIFEST;

    @Argument(shortName = "CF", doc = "An Illumina cluster file (egt)")
    public File CLUSTER_FILE;

    @Argument(shortName = "NORM_MANIFEST", doc = "An Illumina bead pool manifest (a manifest containing the Illumina normalization ids) (bpm.csv)")
    public File ILLUMINA_NORMALIZATION_MANIFEST;

    @Argument(shortName = "E_GENDER", doc = "The expected gender for this sample.", optional = true)
    public String EXPECTED_GENDER;

    @Argument(doc = "The sample alias")
    public String SAMPLE_ALIAS;

    @Argument(doc = "The analysis version of the data used to generate this VCF", optional = true)
    public Integer ANALYSIS_VERSION_NUMBER;

    @Argument(shortName = "G_GTC", doc = "An optional GTC file that was generated by calling the chip using a cluster file designed to optimize gender calling.", optional = true)
    public File GENDER_GTC;

    @Argument(shortName = "FP_VCF", doc = "The fingerprint VCF for this sample", optional = true)
    public File FINGERPRINT_GENOTYPES_VCF_FILE;

    @Argument(doc = "Causes the program to fail if it finds a case where there is a call on an assay that is flagged as 'zeroed-out' in the Illumina cluster file.", optional = true)
    public boolean DO_NOT_ALLOW_CALLS_ON_ZEROED_OUT_ASSAYS = false;

    static final List<Allele> NO_CALL_ALLELES = Collections.unmodifiableList(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL));

    // This file gets initialized during customCommandLineValidation.
    // It is a static member so we don't have to parse the file twice.
    private static InfiniumGTCFile infiniumGTCFile;

    private static InfiniumEGTFile infiniumEGTFile;

    private static String gtcGender = null;

    private static Sex fingerprintGender;

    private static ReferenceSequenceFile refSeq;

    private static final DecimalFormat df = new DecimalFormat();

    private static final String DOT = ".";

    static {
        df.setMaximumFractionDigits(3);
    }

    @Override
    protected boolean requiresReference() {
        return true;
    }

    @Override
    protected int doWork() {
        final Build37ExtendedIlluminaManifest manifest = setupAndGetManifest();

        final VCFHeader vcfHeader = createVCFHeader(manifest, infiniumGTCFile, gtcGender, CLUSTER_FILE,
                REFERENCE_SEQUENCE, refSeq.getSequenceDictionary());

        // Setup a collection that will sort contexts properly
        // Necessary because input GTC file is not sorted
        final SortingCollection<VariantContext> contexts =
                SortingCollection.newInstance(
                        VariantContext.class,
                        new VCFRecordCodec(vcfHeader),
                        new VariantContextComparator(refSeq.getSequenceDictionary()),
                        MAX_RECORDS_IN_RAM,
                        TMP_DIR.stream().map(File::toPath).toArray(Path[]::new));

        // fill the sorting collection
        fillContexts(contexts, infiniumGTCFile, manifest, infiniumEGTFile);

        writeVcf(contexts, OUTPUT, refSeq.getSequenceDictionary(), vcfHeader);

        return 0;
    }

    @Override
    protected String[] customCommandLineValidation() {

        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(EXTENDED_ILLUMINA_MANIFEST);
        IOUtil.assertFileIsWritable(OUTPUT);
        refSeq = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);
        final SAMSequenceDictionary sequenceDictionary = refSeq.getSequenceDictionary();
        final String assembly = sequenceDictionary.getSequence(0).getAssembly();
        if (!assembly.equals("GRCh37")) {
            return new String[]{"The selected reference sequence ('" + assembly + "') is not supported.  This tool is currently only implemented to support NCBI Build 37 / HG19 Reference Sequence."};
        }

        if (FINGERPRINT_GENOTYPES_VCF_FILE != null) {
            IOUtil.assertFileIsReadable(FINGERPRINT_GENOTYPES_VCF_FILE);
        }
        if (GENDER_GTC != null) {
            IOUtil.assertFileIsReadable(GENDER_GTC);
        }

        return super.customCommandLineValidation();
    }

    private Build37ExtendedIlluminaManifest setupAndGetManifest() {
        fingerprintGender = getFingerprintSex(FINGERPRINT_GENOTYPES_VCF_FILE);
        final InfiniumNormalizationManifest infiniumNormalizationManifest = new InfiniumNormalizationManifest(ILLUMINA_NORMALIZATION_MANIFEST);
        try (final DataInputStream gtcInputStream = new DataInputStream(new FileInputStream(INPUT))) {

            infiniumEGTFile = new InfiniumEGTFile(CLUSTER_FILE);
            infiniumGTCFile = new InfiniumGTCFile(gtcInputStream, infiniumNormalizationManifest);
            final Build37ExtendedIlluminaManifest manifest = new Build37ExtendedIlluminaManifest(EXTENDED_ILLUMINA_MANIFEST);

            if (GENDER_GTC != null) {
                try (DataInputStream genderGtcStream = new DataInputStream(new FileInputStream(GENDER_GTC))) {
                    gtcGender = new InfiniumGTCFile(genderGtcStream, infiniumNormalizationManifest).getGender();
                }
            }

            final String gtcManifestName = FilenameUtils.removeExtension(infiniumGTCFile.getSnpManifest());
            final String illuminaManifestName = FilenameUtils.removeExtension(manifest.getDescriptorFileName());

            if (!gtcManifestName.equalsIgnoreCase(illuminaManifestName)) {
                throw new PicardException("The GTC's manifest name " + gtcManifestName +
                        " does not match the Illumina manifest name " + illuminaManifestName);
            }

            if (infiniumGTCFile.getNumberOfSnps() != manifest.getNumAssays()) {
                log.warn("The number of Assays in the GTC file: " + infiniumGTCFile.getNumberOfSnps() +
                        " does not equal the number of Assays in the Illumina manifest file: " + manifest.getNumAssays());
            }
            return manifest;
        } catch (IOException e) {
            throw new PicardException("Error during setup", e);
        }
    }

    static Sex getFingerprintSex(final File file) {
        if (file != null) {
            try (VCFFileReader reader = new VCFFileReader(file, false)) {
                final VCFHeader header = reader.getFileHeader();
                final VCFHeaderLine gender = header.getMetaDataLine("gender");
                if (gender != null) {
                    return Sex.valueOf(gender.getValue());
                }
            }
        }
        return Sex.Unknown;
    }

    private void fillContexts(final SortingCollection<VariantContext> contexts, final InfiniumGTCFile gtcFile,
                              final Build37ExtendedIlluminaManifest manifest, final InfiniumEGTFile egtFile) {
        final ProgressLogger progressLogger = new ProgressLogger(log, 100000, "sorted");

        final Iterator<Build37ExtendedIlluminaManifestRecord> iterator = manifest.extendedIterator();
        int gtcIndex = 0;

        int numVariantsWritten = 0;

        while (iterator.hasNext()) {
            final Build37ExtendedIlluminaManifestRecord record = iterator.next();

            if (!record.isBad()) {
                InfiniumGTCRecord gtcRecord = gtcFile.getRecord(gtcIndex);
                VariantContext context = makeVariantContext(record, gtcRecord, egtFile, progressLogger);
                numVariantsWritten++;
                contexts.add(context);
            }
            gtcIndex++;
        }

        log.info(numVariantsWritten + " Variants were written to file");
        log.info(gtcFile.getNumberOfSnps() + " SNPs in the GTC file");
        log.info(manifest.getNumAssays() + " Variants on the " + manifest.getDescriptorFileName() + " genotyping array manifest file");
    }

    private VariantContext makeVariantContext(Build37ExtendedIlluminaManifestRecord record, final InfiniumGTCRecord gtcRecord,
                                              final InfiniumEGTFile egtFile, final ProgressLogger progressLogger) {
        // If the record is not flagged as errant in the manifest we include it in the VCF
        Allele A = record.getAlleleA();
        Allele B = record.getAlleleB();
        Allele ref = record.getRefAllele();

        final String chr = record.getB37Chr();
        final Integer position = record.getB37Pos();
        final Integer endPosition = position + ref.length() - 1;
        Integer egtIndex = egtFile.rsNameToIndex.get(record.getName());
        if (egtIndex == null) {
            throw new PicardException("Found no record in cluster file for manifest entry '" + record.getName() + "'");
        }

        progressLogger.record(chr, position);

        // Create list of unique alleles
        final List<Allele> assayAlleles = new ArrayList<>();
        assayAlleles.add(ref);

        if (A.equals(B)) {
            throw new PicardException("Found same allele (" + A.getDisplayString() + ") for A and B ");
        }

        if (!ref.equals(A, true)) {
            assayAlleles.add(A);
        }

        if (!ref.equals(B, true)) {
            assayAlleles.add(B);
        }

        final String sampleName = FilenameUtils.removeExtension(INPUT.getName());
        final Genotype genotype = getGenotype(sampleName, gtcRecord, record, A, B);

        final VariantContextBuilder builder = new VariantContextBuilder();

        builder.source(record.getName());
        builder.chr(chr);
        builder.start(position);
        builder.stop(endPosition);
        builder.alleles(assayAlleles);
        builder.log10PError(VariantContext.NO_LOG10_PERROR);
        builder.id(record.getName());
        builder.genotypes(genotype);

        VariantContextUtils.calculateChromosomeCounts(builder, false);

        //custom info fields
        builder.attribute(InfiniumVcfFields.ALLELE_A, record.getAlleleA());
        builder.attribute(InfiniumVcfFields.ALLELE_B, record.getAlleleB());
        builder.attribute(InfiniumVcfFields.ILLUMINA_STRAND, record.getIlmnStrand());
        builder.attribute(InfiniumVcfFields.PROBE_A, record.getAlleleAProbeSeq());
        builder.attribute(InfiniumVcfFields.PROBE_B, record.getAlleleBProbeSeq());
        builder.attribute(InfiniumVcfFields.BEADSET_ID, record.getBeadSetId());
        builder.attribute(InfiniumVcfFields.ILLUMINA_CHR, record.getChr());
        builder.attribute(InfiniumVcfFields.ILLUMINA_POS, record.getPosition());
        builder.attribute(InfiniumVcfFields.ILLUMINA_BUILD, record.getGenomeBuild());
        builder.attribute(InfiniumVcfFields.SOURCE, record.getSource().replace(' ', '_'));
        builder.attribute(InfiniumVcfFields.GC_SCORE, formatFloatForVcf(egtFile.totalScore[egtIndex]));

        for (InfiniumVcfFields.GENOTYPE_VALUES gtValue : InfiniumVcfFields.GENOTYPE_VALUES.values()) {
            final int ordinalValue = gtValue.ordinal();
            builder.attribute(InfiniumVcfFields.N[ordinalValue], egtFile.n[egtIndex][ordinalValue]);
            builder.attribute(InfiniumVcfFields.DEV_R[ordinalValue], formatFloatForVcf(egtFile.devR[egtIndex][ordinalValue]));
            builder.attribute(InfiniumVcfFields.MEAN_R[ordinalValue], formatFloatForVcf(egtFile.meanR[egtIndex][ordinalValue]));
            builder.attribute(InfiniumVcfFields.DEV_THETA[ordinalValue], formatFloatForVcf(egtFile.devTheta[egtIndex][ordinalValue]));
            builder.attribute(InfiniumVcfFields.MEAN_THETA[ordinalValue], formatFloatForVcf(egtFile.meanTheta[egtIndex][ordinalValue]));

            final EuclideanValues genotypeEuclideanValues = polarToEuclidean(egtFile.meanR[egtIndex][ordinalValue], egtFile.devR[egtIndex][ordinalValue],
                    egtFile.meanTheta[egtIndex][ordinalValue], egtFile.devTheta[egtIndex][ordinalValue]);
            builder.attribute(InfiniumVcfFields.DEV_X[ordinalValue], formatFloatForVcf(genotypeEuclideanValues.devX));
            builder.attribute(InfiniumVcfFields.MEAN_X[ordinalValue], formatFloatForVcf(genotypeEuclideanValues.meanX));
            builder.attribute(InfiniumVcfFields.DEV_Y[ordinalValue], formatFloatForVcf(genotypeEuclideanValues.devY));
            builder.attribute(InfiniumVcfFields.MEAN_Y[ordinalValue], formatFloatForVcf(genotypeEuclideanValues.meanY));
        }


        final String rsid = record.getRsId();
        if (StringUtils.isNotEmpty(rsid)) {
            builder.attribute(InfiniumVcfFields.RS_ID, rsid);
        }
        if (egtFile.totalScore[egtIndex] == 0.0) {
            builder.filter(InfiniumVcfFields.ZEROED_OUT_ASSAY);
            if (genotype.isCalled()) {
                if (DO_NOT_ALLOW_CALLS_ON_ZEROED_OUT_ASSAYS) {
                    throw new PicardException("Found a call (genotype: " + genotype + ") on a zeroed out Assay!!");
                } else {
                    log.warn("Found a call (genotype: " + genotype + ") on a zeroed out Assay. " +
                            "This could occur if you called genotypes on a different cluster file than used here.");
                }
            }
        }
        if (record.isDupe()) {
            builder.filter(InfiniumVcfFields.DUPE);
        }
        return builder.make();
    }

    //Uses Manhattan distance conversion
    EuclideanValues polarToEuclidean(float r, float rDeviation, float theta, float thetaDeviation) {
        //calculate variance (deviation^2)
        final double thetaVariance = Math.pow(thetaDeviation, 2.0);
        final double rVariance = Math.pow(rDeviation, 2.0);

        final double halfPi = Math.PI / 2;
        // Note that normalizedTheta= is a normal angle measured in radians,
        // while theta has been divided by pi/2 so that it goes from 0 to 1 as normalizedTheta goes from 0 to pi/2
        final double normalizedTheta = halfPi * theta;
        final double rOverX = (1 + Math.tan(normalizedTheta));

        //calculate X and Y variances from R and Theta variances
        final double thetaVarianceFactorX = -1 * (halfPi * r) * Math.pow(rOverX  * Math.cos(normalizedTheta), -2);
        final double rVarianceFactorX = 1 / rOverX;
        final double varianceX = (Math.pow(thetaVarianceFactorX, 2) * thetaVariance) + (Math.pow(rVarianceFactorX, 2) * rVariance);
        final double thetaVarianceFactorY = -1 * thetaVarianceFactorX;
        final double rVarianceFactorY = 1 - rVarianceFactorX;
        final double varianceY = (Math.pow(thetaVarianceFactorY, 2) * thetaVariance) + (Math.pow(rVarianceFactorY, 2) * rVariance);

        /*
            Theta quantifies the relative amount of signal measured by the A and B intensities, defined by the equation:
            (pi/2) * arctan(Y/X). R is a measurement of the total intensity observed from the A and B signals, defined as: R = A+B
            Illumina uses Manhattan distance https://en.wikipedia.org/wiki/Taxicab_geometry which is why R is A+B and not sqrt(A^2 + B^2)
            So Theta = (2/pi) * arctan(Y/X) and R = X + Y
         */

        final double meanX = r / rOverX;
        final double meanY = r - meanX;
        final double devX = Math.pow(varianceX, 0.5);
        final double devY = Math.pow(varianceY, 0.5);

        return new EuclideanValues((float) meanX, (float) meanY, (float) devX, (float) devY);
    }

    class EuclideanValues {
        final float meanX, meanY, devX, devY;

        EuclideanValues(float meanX, float meanY, float devX, float devY) {
            this.meanX = meanX;
            this.meanY = meanY;
            this.devX = devX;
            this.devY = devY;
        }
    }

    public Genotype getGenotype(final String sampleName,
                                final InfiniumGTCRecord infiniumGtcRecord,
                                final IlluminaManifestRecord record,
                                final Allele A,
                                final Allele B) {

        // The Sample Alleles
        final List<Allele> alleles;

        if (infiniumGtcRecord.genotype == InfiniumGTCFile.NO_CALL) alleles = NO_CALL_ALLELES;
        else if (infiniumGtcRecord.genotype == InfiniumGTCFile.AA_CALL) alleles = Arrays.asList(A, A);
        else if (infiniumGtcRecord.genotype == InfiniumGTCFile.AB_CALL) alleles = Arrays.asList(A, B);
        else if (infiniumGtcRecord.genotype == InfiniumGTCFile.BB_CALL) alleles = Arrays.asList(B, B);
        else {
            throw new PicardException("Unexpected genotype call [" + infiniumGtcRecord.genotype + "]" + " for SNP: " + record.getName());
        }

        final Map<String, Object> attributes = new HashMap<>();
        attributes.put(InfiniumVcfFields.IGC, formatFloatForVcf(infiniumGtcRecord.genotypeScore));
        attributes.put(InfiniumVcfFields.X, infiniumGtcRecord.rawXIntensity);
        attributes.put(InfiniumVcfFields.Y, infiniumGtcRecord.rawYIntensity);
        attributes.put(InfiniumVcfFields.NORMX, formatFloatForVcf(infiniumGtcRecord.normalizedXIntensity));
        attributes.put(InfiniumVcfFields.NORMY, formatFloatForVcf(infiniumGtcRecord.normalizedYIntensity));
        attributes.put(InfiniumVcfFields.R, formatFloatForVcf(infiniumGtcRecord.RIlmn));
        attributes.put(InfiniumVcfFields.THETA, formatFloatForVcf(infiniumGtcRecord.thetaIlmn));
        attributes.put(InfiniumVcfFields.BAF, formatFloatForVcf(infiniumGtcRecord.bAlleleFreq));
        attributes.put(InfiniumVcfFields.LRR, formatFloatForVcf(infiniumGtcRecord.logRRatio));

        return GenotypeBuilder.create(sampleName, alleles, attributes);
    }

    public static String formatFloatForVcf(final float value) {
        if (Float.isNaN(value)) {
            return DOT;
        }
        return df.format(value);
    }

    /**
     * Writes out the VariantContext objects in the order presented to the supplied output file
     * in VCF format.
     */
    private void writeVcf(final SortingCollection<VariantContext> variants,
                          final File output,
                          final SAMSequenceDictionary dict,
                          final VCFHeader vcfHeader) {

        try (final VariantContextWriter writer = new VariantContextWriterBuilder()
                .setOutputFile(output)
                .setReferenceDictionary(dict)
                .setOptions(VariantContextWriterBuilder.DEFAULT_OPTIONS)
                .build()) {

            writer.writeHeader(vcfHeader);

            for (final VariantContext variant : variants) {
                if (variant.getAlternateAlleles().size() > 1) {
                    variant.getCommonInfo().addFilter(InfiniumVcfFields.TRIALLELIC);
                }

                writer.add(variant);
            }
        }
    }

    private VCFHeader createVCFHeader(final Build37ExtendedIlluminaManifest manifest,
                                      final InfiniumGTCFile gtcFile,
                                      final String gtcGender,
                                      final File clusterFile,
                                      final File reference,
                                      final SAMSequenceDictionary dict) {
        final String inputName = INPUT.getName();
        final String chipWellBarcode = inputName.substring(0, inputName.lastIndexOf('.'));

        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFHeaderLine("fileDate", new Date().toString()));
        lines.add(new VCFHeaderLine("source", "GtcToVcf"));
        final String descriptorFileName = manifest.getDescriptorFileName();
        lines.add(new VCFHeaderLine(InfiniumVcfFields.ARRAY_TYPE, descriptorFileName.substring(0, descriptorFileName.lastIndexOf(DOT))));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.EXTENDED_ILLUMINA_MANIFEST_FILE, EXTENDED_ILLUMINA_MANIFEST.getName()));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.EXTENDED_ILLUMINA_MANIFEST_VERSION, manifest.getExtendedManifestVersion()));

        lines.add(new VCFHeaderLine(InfiniumVcfFields.CHIP_WELL_BARCODE, chipWellBarcode));
        if (ANALYSIS_VERSION_NUMBER != null) {
            lines.add(new VCFHeaderLine(InfiniumVcfFields.ANALYSIS_VERSION_NUMBER, ANALYSIS_VERSION_NUMBER.toString()));
        }
        lines.add(new VCFHeaderLine(InfiniumVcfFields.SAMPLE_ALIAS, SAMPLE_ALIAS));
        if (EXPECTED_GENDER != null) {
            lines.add(new VCFHeaderLine(InfiniumVcfFields.EXPECTED_GENDER, EXPECTED_GENDER));
        }
        //add control codes
        final int measurementCount = gtcFile.getRawControlXIntensities().length / ArraysControlInfo.CONTROL_INFO.length;
        for (int i = 0; i < ArraysControlInfo.CONTROL_INFO.length; i++) {
            final int offset = i * measurementCount;
            ArraysControlInfo controlInfo = ArraysControlInfo.CONTROL_INFO[i];
            final int redIntensity = gtcFile.getRawControlXIntensity(offset);
            final int greenIntensity = gtcFile.getRawControlYIntensity(offset);
            lines.add(new VCFHeaderLine(controlInfo.getControl(), controlInfo.toString() + "|" + redIntensity + "|" + greenIntensity));
        }
        lines.add(new VCFHeaderLine(InfiniumVcfFields.FINGERPRINT_GENDER, fingerprintGender.name()));
        if (gtcGender != null) {
            lines.add(new VCFHeaderLine(InfiniumVcfFields.AUTOCALL_GENDER, gtcGender));
        } else {
            lines.add(new VCFHeaderLine(InfiniumVcfFields.AUTOCALL_GENDER, gtcFile.getGender()));
        }
        lines.add(new VCFHeaderLine(InfiniumVcfFields.AUTOCALL_DATE, gtcFile.getAutoCallDate()));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.IMAGING_DATE, gtcFile.getImagingDate()));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.CLUSTER_FILE, clusterFile.getName()));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.MANIFEST_FILE, descriptorFileName));
        lines.add(new VCFHeaderLine("content", manifest.getManifestFile().getName()));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.AUTOCALL_VERSION, gtcFile.getAutoCallVersion()));
        lines.add(new VCFHeaderLine("reference", reference.getAbsolutePath()));
        lines.add(new VCFHeaderLine("picardVersion", this.getVersion()));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.P_95_RED, String.valueOf(gtcFile.getP95Red())));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.P_95_GREEN, String.valueOf(gtcFile.getP95Green())));
        lines.add(new VCFHeaderLine(InfiniumVcfFields.SCANNER_NAME, gtcFile.getScannerName()));

        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));
        lines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
        lines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));
        lines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));

        lines.add(new VCFFormatHeaderLine(InfiniumVcfFields.IGC, 1, VCFHeaderLineType.Float, "Illumina GenCall Confidence Score"));
        lines.add(new VCFFormatHeaderLine(InfiniumVcfFields.X, 1, VCFHeaderLineType.Integer, "Raw X intensity"));
        lines.add(new VCFFormatHeaderLine(InfiniumVcfFields.Y, 1, VCFHeaderLineType.Integer, "Raw Y intensity"));
        lines.add(new VCFFormatHeaderLine(InfiniumVcfFields.NORMX, 1, VCFHeaderLineType.Float, "Normalized X intensity"));
        lines.add(new VCFFormatHeaderLine(InfiniumVcfFields.NORMY, 1, VCFHeaderLineType.Float, "Normalized Y intensity"));
        lines.add(new VCFFormatHeaderLine(InfiniumVcfFields.R, 1, VCFHeaderLineType.Float, "Normalized R value"));
        lines.add(new VCFFormatHeaderLine(InfiniumVcfFields.THETA, 1, VCFHeaderLineType.Float, "Normalized Theta value"));
        lines.add(new VCFFormatHeaderLine(InfiniumVcfFields.BAF, 1, VCFHeaderLineType.Float, "B Allele Frequency"));
        lines.add(new VCFFormatHeaderLine(InfiniumVcfFields.LRR, 1, VCFHeaderLineType.Float, "Log R Ratio"));

        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.ALLELE_A, 1, VCFHeaderLineType.String, "A allele"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.ALLELE_B, 1, VCFHeaderLineType.String, "B allele"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.ILLUMINA_STRAND, 1, VCFHeaderLineType.String, "Probe strand"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.PROBE_A, 1, VCFHeaderLineType.String, "Probe base pair sequence"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.PROBE_B, 1, VCFHeaderLineType.String, "Probe base pair sequence; not missing for strand-ambiguous SNPs"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.BEADSET_ID, 1, VCFHeaderLineType.Integer, "Bead set ID for normalization"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.ILLUMINA_CHR, 1, VCFHeaderLineType.String, "Chromosome in Illumina manifest"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.ILLUMINA_POS, 1, VCFHeaderLineType.Integer, "Position in Illumina manifest"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.ILLUMINA_BUILD, 1, VCFHeaderLineType.String, "Genome Build in Illumina manifest"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.SOURCE, 1, VCFHeaderLineType.String, "Probe source"));
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.GC_SCORE, 1, VCFHeaderLineType.Float, "Gentrain Score"));
        for (InfiniumVcfFields.GENOTYPE_VALUES gtValue : InfiniumVcfFields.GENOTYPE_VALUES.values()) {
            final int ordinalValue = gtValue.ordinal();
            lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.N[ordinalValue], 1, VCFHeaderLineType.Integer, "Number of " + gtValue.name() +" calls in training set"));
            lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.DEV_R[ordinalValue], 1, VCFHeaderLineType.Float, "Standard deviation of normalized R for " + gtValue.name() + " cluster"));
            lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.DEV_THETA[ordinalValue], 1, VCFHeaderLineType.Float, "Standard deviation of normalized THETA for " + gtValue.name() + " cluster"));
            lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.DEV_X[ordinalValue], 1, VCFHeaderLineType.Float, "Standard deviation of normalized X for " + gtValue.name() +" cluster"));
            lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.DEV_Y[ordinalValue], 1, VCFHeaderLineType.Float, "Standard deviation of normalized Y for " + gtValue.name() +" cluster"));
            lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.MEAN_R[ordinalValue], 1, VCFHeaderLineType.Float, "Mean of normalized R for " + gtValue.name() + " cluster"));
            lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.MEAN_THETA[ordinalValue], 1, VCFHeaderLineType.Float, "Mean of normalized THETA for " + gtValue.name() + " cluster"));
            lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.MEAN_X[ordinalValue], 1, VCFHeaderLineType.Float, "Mean of normalized X for " + gtValue.name() +" cluster"));
            lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.MEAN_Y[ordinalValue], 1, VCFHeaderLineType.Float, "Mean of normalized Y for " + gtValue.name() +" cluster"));
        }
        lines.add(new VCFInfoHeaderLine(InfiniumVcfFields.RS_ID, 1, VCFHeaderLineType.String, "dbSNP rsID"));

        lines.add(new VCFFilterHeaderLine(InfiniumVcfFields.DUPE, "Duplicate assays position."));
        lines.add(new VCFFilterHeaderLine(InfiniumVcfFields.TRIALLELIC, "Tri-allelic assay."));
        lines.add(new VCFFilterHeaderLine(InfiniumVcfFields.FAIL_REF, "Assay failed to map to reference."));
        lines.add(new VCFFilterHeaderLine(InfiniumVcfFields.ZEROED_OUT_ASSAY, "Assay Zeroed out (marked as uncallable) in the Illumina Cluster File"));

        final VCFHeader header = new VCFHeader(lines, Collections.singletonList(chipWellBarcode));
        header.setSequenceDictionary(dict);
        return header;
    }
}