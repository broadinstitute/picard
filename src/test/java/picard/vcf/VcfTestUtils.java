package picard.vcf;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import java.io.File;
import java.io.IOException;
import java.util.EnumSet;

public class VcfTestUtils {
    /**
     * Useful test method.  Creates a (temporary) indexed VCF so that we don't have to store the index file in the testdata set.
     * @param vcfFile the vcf file to index
     * @return File a vcf file (index file is created in same path).
     */
    public static File createIndexedVcf(final File vcfFile, final String tempFilePrefix) throws IOException {
        final File output = File.createTempFile(tempFilePrefix, ".vcf");
        output.deleteOnExit();
        final File indexFile = new File(output.getAbsolutePath() + ".idx");
        indexFile.deleteOnExit();
        final VCFFileReader in = new VCFFileReader(vcfFile, false);
        final VCFHeader header = in.getFileHeader();

        final VariantContextWriter out = new VariantContextWriterBuilder().
                setReferenceDictionary(header.getSequenceDictionary()).
                setOptions(EnumSet.of(Options.INDEX_ON_THE_FLY)).
                setOutputFile(output).build();
        out.writeHeader(header);
        for (final VariantContext ctx : in) {
            out.add(ctx);
        }
        out.close();
        in.close();
        return output;
    }
}
