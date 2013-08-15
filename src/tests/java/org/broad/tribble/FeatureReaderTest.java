package org.broad.tribble;

import net.sf.samtools.seekablestream.SeekableFileStream;
import net.sf.samtools.util.CloserUtil;
import org.broad.tribble.bed.BEDCodec;
import org.broad.tribble.example.ExampleBinaryCodec;
import org.broad.tribble.index.Block;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.readers.LocationAware;
import org.broad.tribble.util.ParsingUtils;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;


public class FeatureReaderTest {
    private final static File asciiBedFile = new File(TestUtils.DATA_DIR + "test.bed");
    private final static File binaryBedFile = new File(TestUtils.OUTPUT_DIR + "test.binary.bed");
    private final static File tabixBedFile = new File(TestUtils.DATA_DIR + "test.tabix.bed.gz");

    @BeforeClass
    public void setup() throws IOException {
        ExampleBinaryCodec.convertToBinaryTest(asciiBedFile, binaryBedFile, new BEDCodec());
        binaryBedFile.deleteOnExit();
    }

    @AfterClass
    public void tearDown() throws Exception {
        binaryBedFile.delete();
    }

    @DataProvider(name = "indexProvider")
    public Object[][] createData1() {
        return new Object[][]{
                {asciiBedFile, IndexFactory.IndexType.LINEAR, new BEDCodec()},
                {asciiBedFile, IndexFactory.IndexType.INTERVAL_TREE, new BEDCodec()},
                {tabixBedFile, IndexFactory.IndexType.TABIX, new BEDCodec()},
                {binaryBedFile, IndexFactory.IndexType.LINEAR, new ExampleBinaryCodec()},
                {binaryBedFile, IndexFactory.IndexType.INTERVAL_TREE, new ExampleBinaryCodec()},
        };
    }

    @Test(dataProvider = "indexProvider")
    public void testBedQuery(final File featureFile, final IndexFactory.IndexType indexType, final FeatureCodec<Feature, LocationAware> codec) throws IOException {
        final AbstractFeatureReader<Feature, ?> reader = getReader(featureFile, indexType, codec);

        // Query
        testQuery(reader, "chr1", 1, 500, 3);
        testQuery(reader, "chr1", 1, 200, 1);
        testQuery(reader, "chr1", 1, 201, 2);
        testQuery(reader, "chr1", 500, 600, 0);
        testQuery(reader, "chr1", 100000, 100010, 1);
        testQuery(reader, "chr1", 100000, 100000, 0);
        testQuery(reader, "chr1", 100001, 100001, 1);
        testQuery(reader, "chr1", 100005, 100006, 1);
        testQuery(reader, "chr1", 100009, 100011, 1);
        testQuery(reader, "chr1", 100010, 100010, 1);
        testQuery(reader, "chr1", 100011, 100011, 0);
        testQuery(reader, "chr2", 1, 100, 2);
        testQuery(reader, "chr2", 1, 10, 1);
        testQuery(reader, "chr2", 15, 16, 0);
        testQuery(reader, "chr3", 1, 201, 0);

        // Close reader
        reader.close();
    }

    @Test(dataProvider = "indexProvider")
    public void testLargeNumberOfQueries(final File featureFile, final IndexFactory.IndexType indexType, final FeatureCodec<Feature, LocationAware> codec) throws IOException {
        final AbstractFeatureReader<Feature, LocationAware> reader = getReader(featureFile, indexType, codec);
        for (int i = 0; i < 2000; i++) {
            for (final int start : Arrays.asList(500, 200, 201, 600, 100000)) {
                for (final String chr : Arrays.asList("chr1", "chr2", "chr3")) {
                    CloseableTribbleIterator<Feature> iter = null;
                    try {
                        iter = reader.query(chr, start, start + 1);
                        Assert.assertNotNull(iter, "Failed to create non-null iterator");
                    } finally {
                        CloserUtil.close(iter);
                    }
                }
            }
        }

        // Close reader
        reader.close();
    }

    private void testQuery(final AbstractFeatureReader<Feature, ?> reader, final String chr, final int start, final int stop, final int expectedNumRecords) throws IOException {
        final Iterator<Feature> iter = reader.query(chr, start, stop);
        int count = 0;
        while (iter.hasNext()) {
            final Feature f = iter.next();
            Assert.assertTrue(f.getEnd() >= start && f.getStart() <= stop);
            count++;
        }
        Assert.assertEquals(count, expectedNumRecords);
    }

    @Test(dataProvider = "indexProvider")
    public void testBedNames(final File featureFile, final IndexFactory.IndexType indexType, final FeatureCodec<Feature, LocationAware> codec) throws IOException {
        final AbstractFeatureReader<Feature, ?> reader = getReader(featureFile, indexType, codec);
        final String[] expectedSequences = {"chr1", "chr2"};

        final List<String> seqNames = reader.getSequenceNames();
        Assert.assertEquals(seqNames.size(), expectedSequences.length,
                "Expected sequences " + ParsingUtils.join(",", expectedSequences) + " but saw " + ParsingUtils.join(",", seqNames));

        for (final String s : expectedSequences) {
            Assert.assertTrue(seqNames.contains(s));
        }
    }

    private static <FEATURE extends Feature, SOURCE extends LocationAware> AbstractFeatureReader<FEATURE, SOURCE> getReader(final File featureFile,
                                                                                                                            final IndexFactory.IndexType indexType,
                                                                                                                            final FeatureCodec<FEATURE, SOURCE> codec)
            throws IOException {
        if (indexType.canCreate()) {
            // for types we can create make a new index each time
            final File idxFile = Tribble.indexFile(featureFile);

            // delete an already existing index
            if (idxFile.exists()) {
                idxFile.delete();
            }
            final Index idx = IndexFactory.createIndex(featureFile, codec, indexType);
            IndexFactory.writeIndex(idx, idxFile);

            idxFile.deleteOnExit();
        } // else  let's just hope the index exists, and if so use it

        return AbstractFeatureReader.getFeatureReader(featureFile.getAbsolutePath(), codec);
    }

    @Test
    public void testReadingBeyondIntSizedBlock() throws IOException {
        final Block block = new Block(0, ((long) Integer.MAX_VALUE) * 2);
        final SeekableFileStream stream = new SeekableFileStream(new File("/dev/zero"));
        final TribbleIndexedFeatureReader.BlockStreamWrapper blockStreamWrapper = new TribbleIndexedFeatureReader.BlockStreamWrapper(stream, block);
        final int chunkSize = 100000; // 10 Mb
        final int chunksToRead = (int) Math.ceil(block.getSize() / (chunkSize * 1.0));

        final byte[] bytes = new byte[chunkSize];
        long totalRead = 0;
        for (int chunk = 0; chunk < chunksToRead; chunk++) {
            //System.out.println("Reading chunk " + chunk + " of " + chunkSize + " total read " + totalRead);
            final int nRead = blockStreamWrapper.read(bytes);
            Assert.assertTrue(nRead != -1, "Prematurely got EOF after " + totalRead + " bytes");
            totalRead += nRead;
        }

        // /dev/zero doesn't advance file stream on Linux, so reading never terminates
        // Therefore, we only require a minimum number of bytes
        Assert.assertTrue(totalRead >= block.getSize(), "Failed to read all bytes from a block with size > 2B = " + block.getSize());

    }
}

