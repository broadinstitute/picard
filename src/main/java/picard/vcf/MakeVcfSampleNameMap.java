package picard.vcf;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VariantManipulationProgramGroup;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

/**
 * Creates a TSV from sample name to VCF/GVCF path, with one line per input.
 *
 * <h3>Summary</h3>
 * Creates a TSV from sample name to VCF/GVCF path, with one line per input.
 * Input VCF/GVCFs must contain a header describing exactly one sample.
 * <br/>
 * <h3>Usage example:</h3>
 * <pre>
 *     java -jar picard.jar MakeVcfSampleNameMap \
 *      INPUT=sample1.vcf.gz \
 *      INPUT=sample2.vcf.gz \
 *      OUTPUT=cohort.sample_map
 * </pre>
 *
 * @author Daniel Moran
 */
@CommandLineProgramProperties(
        summary = MakeVcfSampleNameMap.SUMMARY,
        oneLineSummary = MakeVcfSampleNameMap.SHORT_SUMMARY,
        programGroup = VariantManipulationProgramGroup.class
)
public class MakeVcfSampleNameMap extends CommandLineProgram {

    static final String SHORT_SUMMARY = "Creates a TSV from sample name to VCF/GVCF path, with one line per input.";
    static final String SUMMARY = SHORT_SUMMARY + "\n" +
            "Input VCF/GVCFs must contain a header describing exactly one sample.\n" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "    java -jar picard.jar MakeVcfSampleNameMap \\\n" +
            "        INPUT=sample1.vcf.gz \\\n" +
            "        INPUT=sample2.vcf.gz \\\n" +
            "        OUTPUT=cohort.sample_map" +
            "</pre>";

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "One or more input VCFs to extract sample names from.", minElements = 1)
    public List<String> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output file to write the sample-name map to.")
    public File OUTPUT;

    private static final Log log = Log.getInstance(MakeVcfSampleNameMap.class);

    @Override
    protected int doWork() {
        IOUtil.assertFileIsWritable(OUTPUT);

        // For building the output-map without sample-name collisions.
        final Map<String, String> pathToNameMap = new HashMap<>(INPUT.size());
        // For detecting sample-name collisions, so we can warn about them.
        final Map<String, String> nameToPathMap = new HashMap<>(INPUT.size());

        for (final String variantPathString : INPUT) {
            final Path variantPath = getVariantPath(variantPathString);
            final VCFHeader header = getHeaderFromPath(variantPath);

            final int numberOfSamples = header.getNGenotypeSamples();
            if (numberOfSamples != 1) {
                throw new PicardException("Input: " + variantPathString + " was expected to contain a single sample" +
                        " but actually contained " + numberOfSamples + " samples.");
            }

            final String sampleName = header.getGenotypeSamples().get(0);
            pathToNameMap.put(variantPathString, sampleName);
            final String previousPath = nameToPathMap.put(sampleName, variantPathString);
            if (previousPath != null) {
                log.warn("Duplicate sample: " + sampleName + ". Sample was found in both " +
                        previousPath + " and in " + variantPathString);
            }
        }

        final List<String> mapLines = new ArrayList<>(INPUT.size());
        pathToNameMap.forEach((path, name) -> mapLines.add(path + "\t" + name));

        try {
            Files.write(OUTPUT.toPath(), mapLines);
        } catch (IOException e) {
            throw new PicardException("Error while writing output", e);
        }

        return 0;
    }

    private static Path getVariantPath(final String pathString) {
        try {
            return IOUtil.getPath(pathString);
        } catch (IOException e) {
            throw new PicardException("Error while converting input " + pathString + " to Path", e);
        }
    }

    private static VCFHeader getHeaderFromPath(final Path variantPath) {
        try(final AbstractFeatureReader<VariantContext, LineIterator> reader = getReaderFromPath(variantPath)) {
            final VCFHeader header = (VCFHeader) reader.getHeader();
            if (header == null) {
                throw new PicardException("Null header found in " + variantPath.toUri() + ".");
            }
            return header;
        } catch (final IOException e) {
            throw new PicardException("Error while reading VCF header from " + variantPath.toUri(), e);
        }
    }

    private static AbstractFeatureReader<VariantContext, LineIterator> getReaderFromPath(final Path variantPath) {
        final String variantURI = variantPath.toAbsolutePath().toUri().toString();
        try {
            return AbstractFeatureReader.getFeatureReader(variantURI, null, new VCFCodec(),
                    false, Function.identity(), Function.identity());
        } catch (final TribbleException e) {
            throw new PicardException("Failed to create reader from " + variantURI, e);
        }
    }
}
