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
import htsjdk.variant.vcf.VCFFilterHeaderLine;
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

/**
 * Class to take genotype calls from a ped file output from zCall and merge them into a vcf from autocall.
 */
@CommandLineProgramProperties(
        summary = "Program to merge a ped file from zCall into a VCF from autocall.",
        oneLineSummary = "Program to merge a ped file from zCall into a VCF.",
        programGroup = picard.cmdline.programgroups.GenotypingArraysProgramGroup.class
)
public class MergePedIntoVcf extends CommandLineProgram {

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

    @Argument(shortName = "ZCALL_VERSION", doc = "The version of zcall used", optional = true)
    public String ZCALL_VERSION = null;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output VCF file to write with merged genotype calls.")
    public File OUTPUT;

    private static final String DOT = ".";

    private static HashMap<String, String[]> zCallThresholds = new HashMap<>();

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(PED_FILE);
        IOUtil.assertFileIsReadable(MAP_FILE);
        IOUtil.assertFileIsReadable(ZCALL_THRESHOLDS_FILE);
        IOUtil.assertFileIsWritable(OUTPUT);
        try {
            parseZCallThresholds();
            ZCallPedFile zCallPedFile = ZCallPedFile.fromFile(PED_FILE, MAP_FILE);
            VCFFileReader vcfFileReader = new VCFFileReader(ORIGINAL_VCF, false);
            VCFHeader header = vcfFileReader.getFileHeader();
            addAdditionalHeaderFields(header);
            writeVcf(vcfFileReader.iterator(), OUTPUT, vcfFileReader.getFileHeader().getSequenceDictionary(), header, zCallPedFile);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        return 0;
    }

    private void parseZCallThresholds() {
        final String notApplicable = "NA";
        try (Stream<String> stream = Files.lines(ZCALL_THRESHOLDS_FILE.toPath())) {

            stream.forEach(line -> {
                String[] tokens = line.split("\t");
                if ((!tokens[1].equals(notApplicable)) && (!tokens[2].equals(notApplicable))) {
                    zCallThresholds.put(tokens[0], new String[]{tokens[1], tokens[2]});
                } else {
                    zCallThresholds.put(tokens[0], new String[]{DOT, DOT});
                }
            });

        } catch (IOException e) {
            throw new PicardException("Error parsing ZCall Thresholds File", e);
        }
    }

    private void addAdditionalHeaderFields(VCFHeader header) {
        header.addMetaDataLine(new VCFHeaderLine(InfiniumVcfFields.ZCALL_VERSION, ZCALL_VERSION));
        header.addMetaDataLine(new VCFHeaderLine(InfiniumVcfFields.ZCALL_THRESHOLDS, ZCALL_THRESHOLDS_FILE.getName()));
        header.addMetaDataLine(new VCFFilterHeaderLine(InfiniumVcfFields.ZCALL_DIFF));
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

        final VariantContextWriter writer = new VariantContextWriterBuilder()
                .setOutputFile(output)
                .setReferenceDictionary(dict)
                .setOptions(VariantContextWriterBuilder.DEFAULT_OPTIONS)
                .build();

        writer.writeHeader(vcfHeader);
        while (variants.hasNext()) {
            VariantContext context = variants.next();

            final VariantContextBuilder builder = new VariantContextBuilder(context);
            if (zCallThresholds.containsKey(context.getID())) {
                String[] zThresh = zCallThresholds.get(context.getID());
                builder.attribute(InfiniumVcfFields.ZTHRESH_X, zThresh[0]);
                builder.attribute(InfiniumVcfFields.ZTHRESH_Y, zThresh[1]);
            }
            final Genotype originalGenotype = context.getGenotype(0);
            final Map<String, Object> newAttributes = originalGenotype.getExtendedAttributes();
            VCFEncoder vcfEncoder = new VCFEncoder(vcfHeader, false, false);
            Map<Allele, String> alleleMap = vcfEncoder.buildAlleleStrings(context);

            String zCallAlleles = zCallPedFile.getAlleles(context.getID());
            if (zCallAlleles == null) {
                throw new PicardException("No zCall alleles found for snp " + context.getID());
            }
            List<Allele> zCallPedFileAlleles = buildNewAllelesFromZCall(zCallAlleles, context.getAttributes());
            newAttributes.put(InfiniumVcfFields.GTA, alleleMap.get(originalGenotype.getAllele(0)) + "/" + alleleMap.get(originalGenotype.getAllele(1)));
            newAttributes.put(InfiniumVcfFields.GTZ, alleleMap.get(zCallPedFileAlleles.get(0)) + "/" + alleleMap.get(zCallPedFileAlleles.get(1)));

            final Genotype newGenotype = GenotypeBuilder.create(originalGenotype.getSampleName(), zCallPedFileAlleles,
                    newAttributes);
            builder.genotypes(newGenotype);
            logger.record("0", 0);
            // AC, AF, and AN are recalculated here
            VariantContextUtils.calculateChromosomeCounts(builder, false);
            VariantContext newContext = builder.make();
            if (!zCallPedFileAlleles.equals(originalGenotype.getAlleles())) {
                newContext.getCommonInfo().addFilter(InfiniumVcfFields.ZCALL_DIFF);
            }
            writer.add(newContext);
        }

        writer.close();
    }

    private List<Allele> buildNewAllelesFromZCall(String zCallPedFileAlleles, Map<String, Object> newAttributes) {
        char allele1 = zCallPedFileAlleles.charAt(0);
        char allele2 = zCallPedFileAlleles.charAt(1);
        List<Allele> newAlleles = new ArrayList<>();
        String alleleA = String.valueOf(newAttributes.get("ALLELE_A"));
        String alleleB = String.valueOf(newAttributes.get("ALLELE_B"));
        newAlleles.add(translateAllele(alleleA, alleleB, allele1));
        newAlleles.add(translateAllele(alleleA, alleleB, allele2));
        return newAlleles;
    }

    private Allele translateAllele(String alleleA, String alleleB, char allele) {
        if (allele == 'A') {
            return convertIndels(alleleA);
        } else if (allele == 'B') {
            return convertIndels(alleleB);
        } else if (allele == '0') {
            return Allele.NO_CALL;
        } else {
            throw new PicardException("Illegal allele: " + allele);
        }
    }

    private Allele convertIndels(String alleleA) {
        if (alleleA.equals("*")) {
            return Allele.SPAN_DEL;
        } else if (alleleA.contains("*")) {
            return Allele.create(alleleA.replace('*', ' ').trim(), true);
        } else {
            return Allele.create(alleleA, false);
        }
    }
}
