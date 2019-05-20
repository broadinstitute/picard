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
    public void testGraph() throws IOException {
        //graph is based on example in https://en.wikipedia.org/wiki/Clique_problem
        final BitSet[] ajacencyMatrix = new BitSet[7];

        BitSet zero = new BitSet();
        ajacencyMatrix[0]=zero;

        BitSet one = new BitSet();
        one.set(2);
        one.set(5);
        ajacencyMatrix[1]=one;

        BitSet two = new BitSet();
        two.set(1);
        two.set(5);
        two.set(3);
        ajacencyMatrix[2]=two;

        BitSet three = new BitSet();
        three.set(2);
        three.set(4);
        ajacencyMatrix[3]=three;

        BitSet four = new BitSet();
        four.set(3);
        four.set(5);
        four.set(6);
        ajacencyMatrix[4]=four;

        BitSet five = new BitSet();
        five.set(1);
        five.set(2);
        five.set(4);
        ajacencyMatrix[5]=five;

        BitSet six = new BitSet();
        six.set(4);
        ajacencyMatrix[6]=six;

        SelectBarcodes.barcodes.clear();
        SelectBarcodes.barcodes.addAll(CollectionUtil.makeList("zero", "one", "two", "three", "four", "five", "six"));

        SelectBarcodes.output = File.createTempFile("testing", "txt");
        SelectBarcodes.output.deleteOnExit();
        SelectBarcodes.find_cliques(ajacencyMatrix, new BitSet());

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
    public void testGraph2() throws IOException {
        //test based on the C(7,4) graph in https://en.wikipedia.org/wiki/Clique_problem
        final BitSet[] ajacencyMatrix = new BitSet[7];

        BitSet zero = new BitSet();
        zero.set(1);
        zero.set(2);
        zero.set(4);
        ajacencyMatrix[0]=zero;

        BitSet one = new BitSet();
        one.set(0);
        one.set(4);
        one.set(3);
        ajacencyMatrix[1]=one;

        BitSet two = new BitSet();
        two.set(0);
        two.set(3);
        two.set(5);
        two.set(6);
        ajacencyMatrix[2]=two;

        BitSet three = new BitSet();
        three.set(2);
        three.set(1);
        three.set(5);
        three.set(6);
        ajacencyMatrix[3]=three;

        BitSet four = new BitSet();
        four.set(0);
        four.set(1);
        four.set(5);
        four.set(6);
        ajacencyMatrix[4]=four;

        BitSet five = new BitSet();
        five.set(4);
        five.set(6);
        five.set(2);
        five.set(3);
        ajacencyMatrix[5]=five;

        BitSet six = new BitSet();
        six.set(5);
        six.set(3);
        six.set(2);
        six.set(4);
        ajacencyMatrix[6]=six;

        SelectBarcodes.barcodes.clear();
        SelectBarcodes.barcodes.addAll(CollectionUtil.makeList("zero", "one", "two", "three", "four", "five", "six"));

        SelectBarcodes.output = File.createTempFile("testing", "txt");
        SelectBarcodes.output.deleteOnExit();
        SelectBarcodes.find_cliques(ajacencyMatrix, new BitSet());

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
        Assert.assertEquals(result, CollectionUtil.makeList("2\ttwo","3\tthree", "5\tfive","6\tsix"));
    }

    @Test
    public void testDistance() {
        // ATGGTATT
        // |||.|||.
        // ATGAATTA
        Assert.assertEquals(SelectBarcodes.levenshtein("ATGGTATT", "ATGAATTA",  6), 2);
    }

    @Test
    public void testDistance2() {

        // ATGAAGTG
        // |||||.|.
        // ATGAATTA
        Assert.assertEquals(SelectBarcodes.levenshtein("ATGAAGTG", "ATGAATTA",  6), 2);
    }
}