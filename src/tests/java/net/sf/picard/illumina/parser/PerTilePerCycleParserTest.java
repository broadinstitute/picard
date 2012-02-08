package net.sf.picard.illumina.parser;

import static net.sf.picard.util.CollectionUtil.makeList;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public class PerTilePerCycleParserTest {
    public static final List<Integer> DEFAULT_TILES = makeList(1, 2, 3, 4);
    public static final int [] DEFAULT_OUTPUT_LENGTHS = new int[]{10,5,5};
    public static final int MAX_CYCLE = 20;
    public static final int NUM_CLUSTERS = 20;

    private class MockCycledIlluminaData implements IlluminaData {
        private List<String> values;
        public MockCycledIlluminaData() {
            this.values = new ArrayList<String>();
        }

        public void addValue(String value) {
            values.add(value);
        }

        public List<String> getValues() {
            return values;
        }
    }

    class MockPerTilePerCycleParser extends PerTilePerCycleParser<MockCycledIlluminaData> {
        private int [] expectedOutputLengths;

        public MockPerTilePerCycleParser(final File directory, final int lane, final CycleIlluminaFileMap tilesToCycleFiles, final int[] outputLengths) {
            super(directory, lane, tilesToCycleFiles, outputLengths);
            expectedOutputLengths = outputLengths;
        }

        @Override
        protected MockCycledIlluminaData makeData(final int[] outputLengths) {
            Assert.assertEquals(outputLengths, expectedOutputLengths);
            return new MockCycledIlluminaData();

        }

        @Override
        protected CycleFileParser<MockCycledIlluminaData> makeCycleFileParser(final File file, final int cycle) {
            return new CycleFileParser<MockCycledIlluminaData>() {
                int currentCluster = 0;
                final String filePrefix = str_del(file.getName(), cycle);

                @Override
                public void close() {
                }

                @Override
                public void next(final MockCycledIlluminaData ild) {
                    if(!hasNext()) {
                        throw new NoSuchElementException();
                    }

                    ild.addValue(str_del(filePrefix, currentCluster++));
                }

                @Override
                public boolean hasNext() {
                    return currentCluster < NUM_CLUSTERS;
                }
            };
        }

        public Set<IlluminaDataType> supportedTypes() {
            return null;
        }
    }

    //Provides an iterator over files for one cycle for one multiple tiles
    class MontonicCycleFilesIterator extends CycleFilesIterator {
        private Iterator<String> iterator;
        private final List<String> fileNames;

        public MontonicCycleFilesIterator(final List<String> fileNames) {
            super(null, 0, 0, null);
            this.iterator = fileNames.iterator();
            this.fileNames = fileNames;
        }

        @Override
        public void reset() {
            this.iterator = fileNames.iterator();
        }

        @Override
        public File next() {
            return new File(iterator.next());
        }

        @Override
        public boolean hasNext() {
            return iterator.hasNext();
        }
    }

    public List<String> getFileNames(final int tile, final int cluster, final int maxCycle) {
        final List<String> fileNames = new ArrayList<String>();
        for(int i = 1; i <= maxCycle; i++) {
            fileNames.add(str_del(tile, i ,i, cluster));
        }
        return fileNames;
    }

    public List<String> getFileNames(final List<Integer> tiles, final int numClusters) {
        final List<String> fileNames = new ArrayList<String>();
        for(final Integer tile : tiles) {
            for(int i = 0; i < numClusters; i++) {
                fileNames.addAll(getFileNames(tile, i, MAX_CYCLE));
            }
        }

        return fileNames;
    }

    public List<CycleFilesIterator> getCycleFileIterators(final List<Integer> tiles) {
        final List<CycleFilesIterator> iterators = new ArrayList<CycleFilesIterator>();
        for(final Integer tile : tiles) {
            final List<String> fileNames = new ArrayList<String>();
            for(int i = 1; i <= MAX_CYCLE; i++) {
                fileNames.add(str_del(tile,i));
            }

            iterators.add(new MontonicCycleFilesIterator(fileNames));
        }

        return iterators;
    }

    public static String str_del(Object ... objs) {
        String out = objs[0].toString();
        for(int i = 1; i < objs.length; i++) {
            out += "_" + objs[i];
        }
        return out;
    }

    @Test
    public void basicIterationTest() {
        final List<String> expectedValues = getFileNames(DEFAULT_TILES, NUM_CLUSTERS);
        final PerTilePerCycleParser<MockCycledIlluminaData> parser = makeParser();

        int index = 0;
        while(parser.hasNext()) {
            index = compareValues(parser.next().values, expectedValues, index);
        }

        Assert.assertEquals(index, expectedValues.size());
    }


    private int compareValues(final List<String> parserValues, final List<String> expectedValues, int index) {
        for(final String parserValue : parserValues) {
            Assert.assertTrue(index < expectedValues.size());
            Assert.assertEquals(parserValue, expectedValues.get(index), "With index " + index);
            ++index;
        }

        return index;
    }

    public PerTilePerCycleParser<MockCycledIlluminaData> makeParser() {
        CycleIlluminaFileMap fm = new CycleIlluminaFileMap();
        final List<CycleFilesIterator> iterators = getCycleFileIterators(DEFAULT_TILES);
        for(int i = 0; i < iterators.size(); i++) {
            fm.put(DEFAULT_TILES.get(i), iterators.get(i));
        }
        final PerTilePerCycleParser<MockCycledIlluminaData> parser = new MockPerTilePerCycleParser(new File("FakeFile"), 1, fm, DEFAULT_OUTPUT_LENGTHS);
        return parser;
    }

    @DataProvider(name="seekingTests")
    public Object[][] seekingTests() {
        return new Object[][] {
            {1,  3, null, null},
            {22, 1, null, null},
            {38, 2, null, null},
            {75, 4, null, null},
            {1,  3,   70, 1},
            {1,  3,   45, 2},
            {12, 2,   59, 4},
            {45, 3,   70, 3},
            {14, 1,   5,  2}
        };
    }


    @Test(dataProvider = "seekingTests")
    public void seekingIterationTest(Integer seekPos1, Integer newTile1, Integer seekPos2, Integer newTile2) {
        final List<String> expectedValues = getFileNames(DEFAULT_TILES, NUM_CLUSTERS);
        final PerTilePerCycleParser<MockCycledIlluminaData> parser = makeParser();

        int index = 0;
        for(int i = 0; i <= seekPos1; i++) {
            Assert.assertTrue(parser.hasNext());
            index = compareValues(parser.next().values, expectedValues, index);
        }

        parser.seekToTile(newTile1);

        int startCluster = (newTile1-1) * NUM_CLUSTERS;
        index = startCluster * MAX_CYCLE;
        if(seekPos2 != null) {
            for(int i = startCluster; i <= seekPos2; i++) {
                Assert.assertTrue(parser.hasNext());
                index = compareValues(parser.next().values, expectedValues, index);
            }

            parser.seekToTile(newTile2);
            startCluster = (newTile2-1) * NUM_CLUSTERS;
            index = startCluster * MAX_CYCLE;
        }

        for(int i = startCluster; i < NUM_CLUSTERS * DEFAULT_TILES.size(); i++) {
            Assert.assertTrue(parser.hasNext());
            index = compareValues(parser.next().values, expectedValues, index);
        }

        Assert.assertFalse(parser.hasNext());

    }
}
