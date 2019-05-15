package picard.illumina;

import htsjdk.samtools.util.CollectionUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.function.BiFunction;

/**
 * Created by farjoun on 8/22/18.
 */
public class SelectBarcodesTest {

    private Object[][] genericTestMaker(BiFunction<Boolean, Boolean, Boolean> operator) {
        List<Object[]> tests = new ArrayList<>();
        Random rng = new Random(42);

        for (int i = 0; i < 10; i++) {
            BitSet lhs = new BitSet();
            BitSet rhs = new BitSet();
            BitSet expected = new BitSet();

            for (int j = 0; j < 30; j++) {
                final boolean leftBit = rng.nextBoolean();
                final boolean rightBit = rng.nextBoolean();
                final boolean expectedBit = operator.apply(leftBit, rightBit);
                lhs.set(j, leftBit);
                rhs.set(j, rightBit);
                expected.set(j, expectedBit);

            }
            tests.add(new Object[]{lhs, rhs, expected});
        }
        return tests.toArray(new Object[0][]);
    }

    public void testGeneric(final BitSet lhs, final BitSet rhs, final BitSet expected, BiFunction<BitSet, BitSet, BitSet> function) {
        final BitSet origLhs = BitSet.valueOf(lhs.toLongArray());
        final BitSet origRhs = BitSet.valueOf(rhs.toLongArray());

        // check that function get correct value
        Assert.assertEquals(function.apply(origLhs, origRhs), expected);

        // check that function doesn't affect inputs.
        Assert.assertEquals(origLhs, lhs);
        Assert.assertEquals(origRhs, rhs);
    }

    @DataProvider
    public Object[][] testIntersectionData() {
        return genericTestMaker((a, b) -> a & b);
    }

    @Test(dataProvider = "testIntersectionData")
    public void testIntersection(final BitSet lhs, final BitSet rhs, final BitSet expected) {
        testGeneric(lhs, rhs, expected, SelectBarcodes::intersection);
    }

    @DataProvider
    public Object[][] testUnionData() {
        return genericTestMaker((a, b) -> a | b);
    }

    @Test(dataProvider = "testUnionData")
    public void testUnion(final BitSet lhs, final BitSet rhs, final BitSet expected) {
        testGeneric(lhs, rhs, expected, SelectBarcodes::union);
    }

    @DataProvider
    public Object[][] testDifferenceData() {
        return genericTestMaker((a, b) -> a & !b);
    }

    @Test(dataProvider = "testDifferenceData")
    public void testDifference(final BitSet lhs, final BitSet rhs, final BitSet expected) {
        testGeneric(lhs, rhs, expected, SelectBarcodes::difference);
    }

    @Test
    public void test() throws IOException {
        final List<BitSet> ajacencyMatrix = new ArrayList<>();

        BitSet zero = new BitSet();
        ajacencyMatrix.add(zero);

        BitSet one = new BitSet();
        one.set(2);
        one.set(5);
        ajacencyMatrix.add(one);

        BitSet two = new BitSet();
        two.set(1);
        two.set(5);
        two.set(3);
        ajacencyMatrix.add(two);

        BitSet three = new BitSet();
        three.set(2);
        three.set(4);
        ajacencyMatrix.add(three);

        BitSet four = new BitSet();
        four.set(3);
        four.set(5);
        four.set(6);
        ajacencyMatrix.add(four);

        BitSet five = new BitSet();
        five.set(1);
        five.set(2);
        five.set(4);
        ajacencyMatrix.add(five);

        BitSet six = new BitSet();
        six.set(4);
        ajacencyMatrix.add(six);

        SelectBarcodes.barcodes.clear();
        SelectBarcodes.barcodes.addAll(CollectionUtil.makeList("zero", "one", "two", "three", "four", "five", "six"));

        SelectBarcodes.output = File.createTempFile("testing", "txt");
        SelectBarcodes.output.deleteOnExit();
        SelectBarcodes.find_cliques(ajacencyMatrix);

        List<String> result = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(SelectBarcodes.output))) {
            String line;
            while ((line = br.readLine()) != null) {
                result.add(line);
            }
        } catch (
                IOException e) {
            e.printStackTrace();
        }
        Assert.assertEquals(result, CollectionUtil.makeList("1\tone", "2\ttwo", "5\tfive"));
    }

    @Test
    public void testDistance() {
        Assert.assertEquals(SelectBarcodes.levenshtein("AATACCAT", "ATGAATTA", true, 6), 3);
    }

    @Test
    public void testDistance2() {
        Assert.assertEquals(SelectBarcodes.levenshtein("CACTTCAT", "ATGAATTA", true,6), 2);
    }
}