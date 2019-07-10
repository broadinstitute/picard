package picard.util;

import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.testng.Assert;

import java.io.File;
import java.io.IOException;
import java.util.EnumSet;
import java.util.Iterator;

public class VcfTestUtils {
    /**
     * This method makes a copy of the input VCF and creates an index file for it in the same location.
     * This is done so that we don't need to store the index file in the same repo
     * The copy of the input is done so that it and its index are in the same directory which is typically required.
     * @param vcfFile the vcf file to index
     * @return File a vcf file (index file is created in same path).
     */
    public static File createTemporaryIndexedVcfFromInput(final File vcfFile, final String tempFilePrefix) throws IOException {
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

    // These four utility functions have been copied from Picard. It would be best if they were moved from the testing
    // tree to the main tree so that we can use them here...
    public static void assertVcfFilesAreEqual(final File actual, final File expected) throws IOException {

        final File indexedActual = createTemporaryIndexedVcfFromInput(actual, "assert");
        final File indexedExpected = createTemporaryIndexedVcfFromInput(expected, "assert");

        try (final VCFFileReader vcfReaderActual   = new VCFFileReader(indexedActual);
             final VCFFileReader vcfReaderExpected = new VCFFileReader(indexedExpected)) {
            assertEquals(vcfReaderActual, vcfReaderExpected);
            assertEquals(vcfReaderActual.getFileHeader(), vcfReaderExpected.getFileHeader());
        }
    }

    private static void assertEquals(final VariantContext actual, final VariantContext expected) {

        if (expected == null) {
            org.testng.Assert.assertNull(actual);
            return;
        }

        org.testng.Assert.assertNotNull(actual, "null status");
        org.testng.Assert.assertEquals(actual.getContig(), expected.getContig(), "Different contigs: ");
        org.testng.Assert.assertEquals(actual.getStart(), expected.getStart(), "Different starts: ");
        org.testng.Assert.assertEquals(actual.getEnd(), expected.getEnd(), "Different ends: ");

        org.testng.Assert.assertTrue(actual.hasSameAllelesAs(expected), "Alleles differ between " + actual + " and " + expected + ": ");
        assertEquals(actual.getGenotypes(), expected.getGenotypes());

        org.testng.Assert.assertEquals(actual.getID(), expected.getID());
        org.testng.Assert.assertEquals(actual.getFilters(), expected.getFilters());
        org.testng.Assert.assertEquals(actual.getAttributes(), expected.getAttributes(), "");
    }

    private static void assertEquals(final GenotypesContext actual, final GenotypesContext expected) {
        if (expected == null) {
            org.testng.Assert.assertNull(actual);
            return;
        }
        org.testng.Assert.assertEquals(actual.getSampleNamesOrderedByName(), expected.getSampleNamesOrderedByName(), "Sample names differ");

        for (final String name : expected.getSampleNamesOrderedByName()) {
            org.testng.Assert.assertEquals(actual.get(name).getAlleles(), expected.get(name).getAlleles(), "Alleles differ for sample " + name);
            org.testng.Assert.assertEquals(actual.get(name).getAD(), expected.get(name).getAD());
            org.testng.Assert.assertEquals(actual.get(name).getPL(), expected.get(name).getPL());
        }
    }

    private static void assertEquals(final Iterable<VariantContext> actual, final Iterable<VariantContext> expected) {

        final Iterator<VariantContext> actualI = actual.iterator();
        final Iterator<VariantContext> expectedI = expected.iterator();

        org.testng.Assert.assertNotNull(actualI);
        org.testng.Assert.assertNotNull(expectedI);

        while (actualI.hasNext() && expectedI.hasNext()) {
            assertEquals(actualI.next(), expectedI.next());
        }

        org.testng.Assert.assertFalse(actualI.hasNext(), "Found an unexpected variant: " + actualI.next());
        org.testng.Assert.assertFalse(expectedI.hasNext(), "Didn't find an expected variant: " + expectedI.next());
    }

    private static void assertEquals(final VCFHeader actual, final VCFHeader expected) {
        for (VCFHeader.HEADER_FIELDS field: VCFHeader.HEADER_FIELDS.values()) {
            org.testng.Assert.assertEquals(actual.getMetaDataLine(field.name()), expected.getMetaDataLine(field.name()));
        }
        Assert.assertEquals(actual.getMetaDataLine("gender"), expected.getMetaDataLine("gender"));
    }
}
