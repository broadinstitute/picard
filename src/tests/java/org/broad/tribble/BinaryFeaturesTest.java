package org.broad.tribble;

import org.broad.tribble.bed.BEDCodec;
import org.broad.tribble.example.ExampleBinaryCodec;
import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.readers.LineReader;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;


public class BinaryFeaturesTest {
    @DataProvider(name = "BinaryFeatureSources")
    public Object[][] createData1() {
        return new Object[][] {
                { new File(TestUtils.DATA_DIR + "test.bed"),  new BEDCodec() },
                { new File(TestUtils.DATA_DIR + "bed/Unigene.sample.bed"),  new BEDCodec() },
                { new File(TestUtils.DATA_DIR + "bed/NA12878.deletions.10kbp.het.gq99.hand_curated.hg19_fixed.bed"),  new BEDCodec() },
        };
    }

    @Test(enabled = true, dataProvider = "BinaryFeatureSources")
    public void testBinaryCodec(final File source, final FeatureCodec<Feature, LineIterator> codec) throws IOException {
        final File tmpFile = File.createTempFile("testBinaryCodec", ".binary.bed");
        ExampleBinaryCodec.convertToBinaryTest(source, tmpFile, codec);
        tmpFile.deleteOnExit();

        final FeatureReader<Feature> originalReader = AbstractFeatureReader.getFeatureReader(source.getAbsolutePath(), codec, false);
        final FeatureReader<Feature> binaryReader = AbstractFeatureReader.getFeatureReader(tmpFile.getAbsolutePath(), new ExampleBinaryCodec(), false);

        // make sure the header is what we expect
        final List<String> header = (List<String>) binaryReader.getHeader();
        Assert.assertEquals(header.size(), 1, "We expect exactly one header line");
        Assert.assertEquals(header.get(0), ExampleBinaryCodec.HEADER_LINE, "Failed to read binary header line");

        final Iterator<Feature> oit = originalReader.iterator();
        final Iterator<Feature> bit = binaryReader.iterator();
        while ( oit.hasNext() ) {
            final Feature of = oit.next();

            Assert.assertTrue(bit.hasNext(), "Original iterator has items, but there's no items left in binary iterator");
            final Feature bf = bit.next();

            Assert.assertEquals(bf.getChr(), of.getChr(), "Chr not equal between original and binary encoding");
            Assert.assertEquals(bf.getStart(), of.getStart(), "Start not equal between original and binary encoding");
            Assert.assertEquals(bf.getEnd(), of.getEnd(), "End not equal between original and binary encoding");
        }
        Assert.assertTrue(! bit.hasNext(), "Original iterator is done, but there's still some data in binary iterator");

        originalReader.close();
        binaryReader.close();
    }
}
