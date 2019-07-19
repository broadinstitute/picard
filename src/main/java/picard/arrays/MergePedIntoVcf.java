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

import picard.arrays.illumina.InfiniumVcfFields;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import static htsjdk.variant.variantcontext.Allele.SPAN_DEL_STRING;
import static htsjdk.variant.vcf.VCFConstants.MISSING_VALUE_v4;
import static htsjdk.variant.vcf.VCFConstants.UNPHASED;
import static picard.arrays.MergePedIntoVcf.USAGE_DETAILS;
import static picard.arrays.illumina.InfiniumVcfFields.ALLELE_A;
import static picard.arrays.illumina.InfiniumVcfFields.ALLELE_B;

/**
 * Class to take genotype calls from a ped file output from zCall and merge them into a vcf from autocall.
 */
@CommandLineProgramProperties(
        summary = USAGE_DETAILS,
        oneLineSummary = "Program to merge a single-sample ped file from zCall into a single-sample VCF.",
        programGroup = picard.cmdline.programgroups.GenotypingArraysProgramGroup.class
)
public class MergePedIntoVcf extends CommandLineProgram {

    static final String USAGE_DETAILS = "MergePedIntoVcf takes a single-sample ped file output from zCall and merges " +
            "into a single-sample vcf file using several supporting files." +
            "A VCF, aka Variant Calling Format, is a text file for storing how a sequenced sample differs from the reference genome. " +
            "<a href='https://samtools.github.io/hts-specs/VCFv4.2.pdf'></a>" +
            "A PED file is a whitespace-separated text file for storing genotype information. " +
            "<a href='http://zzz.bwh.harvard.edu/plink/data.shtml#ped'></a>" +
            "A MAP file is a whitespace-separated text file for storing information about genetic distance. " +
            "<a href='http://zzz.bwh.harvard.edu/plink/data.shtml#map'></a>" +
            "A zCall thresholds file is a whitespace-separated text file for storing the thresholds for " +
            "genotype clustering for a SNP as determined by zCall." +
            "<a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3463112/#SEC2title'></a>" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar MergePedIntoVcf \\<br />" +
            "      VCF=input.vcf \\<br />" +
            "      PED=zcall.output.ped \\<br />" +
            "      MAP=zcall.output.map \\<br />" +
            "      ZCALL_T_FILE=zcall.thresholds.7.txt \\<br />" +
            "      OUTPUT=output.vcf <br />" +
            "</pre>";

    private final Log log = Log.getInstance(MergePedIntoVcf.class);

    private final ProgressLogger logger = new ProgressLogger(log, 10000);

    @Argument(shortName = "VCF", doc = "The vcf containing the original autocall genotypes.")
    public File ORIGINAL_VCF;

    @Argument(shortName = "PED", doc = "PED file to be merged into VCF.")
    public File PED_FILE;

    @Argument(shortName = "MAP", doc = "MAP file for the PED file.")
    public File MAP_FILE;

    @Argument(shortName = "ZCALL_T_FILE", doc = "The zcall thresholds file.")
    public File ZCALL_THRESHOLDS_FILE = null;

    @Argument(shortName = "ZCALL_VERSION", doc = "The version of zcall used")
    public String ZCALL_VERSION;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output VCF file to write with merged genotype calls.")
    public File OUTPUT;

    private static HashMap<String, String[]> zCallThresholds = new HashMap<>();

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(PED_FILE);
        IOUtil.assertFileIsReadable(MAP_FILE);
        IOUtil.assertFileIsReadable(ZCALL_THRESHOLDS_FILE);
        IOUtil.assertFileIsWritable(OUTPUT);
        try {
            parseZCallThresholds();
            final ZCallPedFile zCallPedFile = ZCallPedFile.fromFile(PED_FILE, MAP_FILE);
            final VCFFileReader vcfFileReader = new VCFFileReader(ORIGINAL_VCF, false);
            final VCFHeader header = vcfFileReader.getFileHeader();
            if (header.getGenotypeSamples().size() > 1) {
                throw new PicardException("MergePedIntoVCF only works with single-sample VCFs.");
            }
            addAdditionalHeaderFields(header);
            writeVcf(vcfFileReader.iterator(), OUTPUT, vcfFileReader.getFileHeader().getSequenceDictionary(), header, zCallPedFile);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return 0;
    }

    private void parseZCallThresholds() {
        final String notApplicable = "NA";
        try (Stream<String> stream = Files.lines(ZCALL_THRESHOLDS_FILE.toPath())) {

            stream.forEach(line -> {
                final String[] tokens = line.split("\t");
                if (tokens[1].equals(notApplicable) || tokens[2].equals(notApplicable)) {
                    if (!tokens[1].equals(notApplicable) || !tokens[2].equals(notApplicable)) {
                        throw new PicardException("Thresholds should either both exist or both not exist.");
                    }
                    zCallThresholds.put(tokens[0], new String[]{MISSING_VALUE_v4, MISSING_VALUE_v4});
                } else {
                    zCallThresholds.put(tokens[0], new String[]{tokens[1], tokens[2]});
                }
            });

        } catch (IOException e) {
            throw new PicardException("Error parsing ZCall Thresholds File", e);
        }
    }

    private void addAdditionalHeaderFields(VCFHeader header) {
        header.addMetaDataLine(new VCFHeaderLine(InfiniumVcfFields.ZCALL_VERSION, ZCALL_VERSION));
        header.addMetaDataLine(new VCFHeaderLine(InfiniumVcfFields.ZCALL_THRESHOLDS, ZCALL_THRESHOLDS_FILE.getName()));
        header.addMetaDataLine(new VCFInfoHeaderLine(InfiniumVcfFields.ZTHRESH_X, 1, VCFHeaderLineType.Float, "zCall X threshold"));
        header.addMetaDataLine(new VCFInfoHeaderLine(InfiniumVcfFields.ZTHRESH_Y, 1, VCFHeaderLineType.Float, "zCall Y threshold"));
        header.addMetaDataLine(new VCFFormatHeaderLine(InfiniumVcfFields.GTA, 1, VCFHeaderLineType.String, "Illumina Autocall Genotype"));
        header.addMetaDataLine(new VCFFormatHeaderLine(InfiniumVcfFields.GTZ, 1, VCFHeaderLineType.String, "zCall Genotype"));
    }

    /**
     * Writes out the VariantContext objects in the order presented to the supplied output file
     * in VCF format.
     */
    private void writeVcf(final CloseableIterator<VariantContext> variants,
                          final File output,
                          final SAMSequenceDictionary dict,
                          final VCFHeader vcfHeader, ZCallPedFile zCallPedFile) {

        try(VariantContextWriter writer = new VariantContextWriterBuilder()
                .setOutputFile(output)
                .setReferenceDictionary(dict)
                .setOptions(VariantContextWriterBuilder.DEFAULT_OPTIONS)
                .build()) {
            writer.writeHeader(vcfHeader);
            while (variants.hasNext()) {
                final VariantContext context = variants.next();

                final VariantContextBuilder builder = new VariantContextBuilder(context);
                if (zCallThresholds.containsKey(context.getID())) {
                    final String[] zThresh = zCallThresholds.get(context.getID());
                    builder.attribute(InfiniumVcfFields.ZTHRESH_X, zThresh[0]);
                    builder.attribute(InfiniumVcfFields.ZTHRESH_Y, zThresh[1]);
                }
                final Genotype originalGenotype = context.getGenotype(0);
                final Map<String, Object> newAttributes = originalGenotype.getExtendedAttributes();
                final VCFEncoder vcfEncoder = new VCFEncoder(vcfHeader, false, false);
                final Map<Allele, String> alleleMap = vcfEncoder.buildAlleleStrings(context);

                final String zCallAlleles = zCallPedFile.getAlleles(context.getID());
                if (zCallAlleles == null) {
                    throw new PicardException("No zCall alleles found for snp " + context.getID());
                }
                final List<Allele> zCallPedFileAlleles = buildNewAllelesFromZCall(zCallAlleles, context.getAttributes());
                newAttributes.put(InfiniumVcfFields.GTA, alleleMap.get(originalGenotype.getAllele(0)) + UNPHASED + alleleMap.get(originalGenotype.getAllele(1)));
                newAttributes.put(InfiniumVcfFields.GTZ, alleleMap.get(zCallPedFileAlleles.get(0)) + UNPHASED + alleleMap.get(zCallPedFileAlleles.get(1)));

                final Genotype newGenotype = GenotypeBuilder.create(originalGenotype.getSampleName(), zCallPedFileAlleles,
                        newAttributes);
                builder.genotypes(newGenotype);
                logger.record("0", 0);
                // AC, AF, and AN are recalculated here
                VariantContextUtils.calculateChromosomeCounts(builder, false);
                final VariantContext newContext = builder.make();
                writer.add(newContext);
            }
        }
    }

    private List<Allele> buildNewAllelesFromZCall(final String zCallPedFileAlleles,
                                                  final Map<String, Object> newAttributes) {
        final char allele1 = zCallPedFileAlleles.charAt(0);
        final char allele2 = zCallPedFileAlleles.charAt(1);
        final List<Allele> newAlleles = new ArrayList<>();
        final String alleleA = String.valueOf(newAttributes.get(ALLELE_A));
        final String alleleB = String.valueOf(newAttributes.get(ALLELE_B));
        newAlleles.add(translateAllele(alleleA, alleleB, allele1));
        newAlleles.add(translateAllele(alleleA, alleleB, allele2));
        return newAlleles;
    }

    private Allele translateAllele(final String alleleA,
                                   final String alleleB,
                                   final char allele) {
        if (allele == 'A') {
            return formatAllele(alleleA);
        } else if (allele == 'B') {
            return formatAllele(alleleB);
        } else if (allele == '0') {
            return Allele.NO_CALL;
        } else {
            throw new PicardException("Illegal allele: " + allele);
        }
    }

    private Allele formatAllele(final String alleleA) {
        if (alleleA.equals(SPAN_DEL_STRING)) {
            return Allele.SPAN_DEL;
        } else if (alleleA.contains(SPAN_DEL_STRING)) {
            return Allele.create(alleleA.replace(SPAN_DEL_STRING.charAt(0), ' ').trim(), true);
        } else {
            return Allele.create(alleleA, false);
        }
    }
}
