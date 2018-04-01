package picard.vcf;

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
     * This method creates a temporary VCF or Bam file and its appropriately named index file, and will delete them on exit.
     * @param prefix - The prefix string to be used in generating the file's name; must be at least three characters long
     * @param suffix - The suffix string to be used in generating the file's name; may be null, in which case the suffix ".tmp" will be used
     * @return A File object referencing the newly created temporary file
     * @throws IOException - if a file could not be created.
     */
    public static File createTemporaryIndexedFile(final String prefix, final String suffix) throws IOException {
        final File out = File.createTempFile(prefix, suffix);
        out.deleteOnExit();
        String indexFileExtension = null;
        if (suffix.endsWith("vcf.gz")) {
            indexFileExtension = ".tbi";
        } else if (suffix.endsWith("vcf")) {
            indexFileExtension = ".idx";
        } else if (suffix.endsWith(".bam")) {
            indexFileExtension = ".bai";
        }

        if (indexFileExtension != null) {
            final File indexOut = new File(out.getAbsolutePath() + indexFileExtension);
            indexOut.deleteOnExit();
        }
        return out;
    }

    /**
     *
     * @param vcfFile
     * @param tempFilePrefix
     * @return
     * @throws IOException
     *
     * @deprecated use createTemporaryIndexedFile instead
     */
    @Deprecated
    public static File createTemporaryIndexedVcfFromInput(final File vcfFile, final String tempFilePrefix) throws IOException {
        return createTemporaryIndexedVcfFromInput(vcfFile, tempFilePrefix, null);
    }

        /**
         * This method makes a copy of the input VCF and creates an index file for it in the same location.
         * This is done so that we don't need to store the index file in the same repo
         * The copy of the input is done so that it and its index are in the same directory which is typically required.
         * @param vcfFile the vcf file to index
         * @return File a vcf file (index file is created in same path).
         */
    public static File createTemporaryIndexedVcfFromInput(final File vcfFile, final String tempFilePrefix, final String suffix) throws IOException {
        final String extension;

        if (suffix != null) {
            extension = suffix;
        } else if (vcfFile.getAbsolutePath().endsWith(".vcf") ) {
            extension = ".vcf";
        } else if (vcfFile.getAbsolutePath().endsWith(".vcf.gz") ) {
            extension = ".vcf.gz";
        } else {
            extension = "nope!";
        }

        if (!extension.equals(".vcf") && !extension.equals(".vcf.gz"))
            throw new IllegalArgumentException("couldn't find a .vcf or .vcf.gz ending for input file " + vcfFile.getAbsolutePath());

        File output = createTemporaryIndexedFile(tempFilePrefix, extension);

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

    public static void assertEquals(final Iterable<VariantContext> actual, final Iterable<VariantContext> expected) {

        final Iterator<VariantContext> actualI = actual.iterator();
        final Iterator<VariantContext> expectedI = expected.iterator();

        Assert.assertNotNull(actualI);
        Assert.assertNotNull(expectedI);

        while (actualI.hasNext() && expectedI.hasNext()) {
            assertEquals(actualI.next(), expectedI.next());
        }

        Assert.assertFalse(actualI.hasNext(), "Found an unexpected variant: " + actualI.next());
        Assert.assertFalse(expectedI.hasNext(), "Didn't find an expected variant: " + expectedI.next());
    }

    public static void assertEquals(final VariantContext actual, final VariantContext expected) {

        if (expected == null) {
            Assert.assertNull(actual);
            return;
        }

        Assert.assertNotNull(actual, "null status");
        Assert.assertEquals(actual.getContig(), expected.getContig(), "Different contigs: ");
        Assert.assertEquals(actual.getStart(), expected.getStart(), "Different starts: ");
        Assert.assertEquals(actual.getEnd(), expected.getEnd(), "Different ends: ");

        Assert.assertTrue(actual.hasSameAllelesAs(expected), "Alleles differ between " + actual + " and " + expected + ": ");
        assertEquals(actual.getGenotypes(), expected.getGenotypes());

        Assert.assertEquals(actual.getID(), expected.getID());
        Assert.assertEquals(actual.getFilters(), expected.getFilters());
        Assert.assertEquals(actual.getAttributes(), expected.getAttributes(), "");
    }

    public static void assertEquals(final GenotypesContext actual, final GenotypesContext expected) {
        if (expected == null) {
            Assert.assertNull(actual);
            return;
        }
        Assert.assertEquals(actual.getSampleNamesOrderedByName(), expected.getSampleNamesOrderedByName(), "Sample names differ");

        for (final String name : expected.getSampleNamesOrderedByName()) {
            Assert.assertEquals(actual.get(name).getAlleles(), expected.get(name).getAlleles(), "Alleles differ for sample " + name);
            Assert.assertEquals(actual.get(name).getAD(), expected.get(name).getAD());
            Assert.assertEquals(actual.get(name).getPL(), expected.get(name).getPL());
        }
    }

    public static void assertVcfFilesAreEqual(final File actual, final File expected) throws IOException {

        final File indexedActual = createTemporaryIndexedVcfFromInput(actual, "assert");
        final File indexedExpected = createTemporaryIndexedVcfFromInput(expected, "assert");

        try (
                final VCFFileReader vcfReaderActual = new VCFFileReader(indexedActual);
                final VCFFileReader vcfReaderExpected = new VCFFileReader(indexedExpected)) {
            VcfTestUtils.assertEquals(vcfReaderActual, vcfReaderExpected);
        }
    }
}
