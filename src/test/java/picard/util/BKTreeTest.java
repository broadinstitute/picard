package picard.util;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

import static org.testng.Assert.*;

public class BKTreeTest {

    @DataProvider(name = "testQueryDataProvider")
    public Object[][] testQueryDataProvider() {
        final Set<String> dictionary = new HashSet<>(Arrays.asList("ACGCA","GCTCG","GGCTA","CGTCG","AGTCG","GGTCG","AACCT"));
        final BKTree<String> tree = new BKTree<>((s1, s2) -> StringDistanceUtils.countMismatches(new byte[][]{s1.getBytes()}, new byte[][]{s2.getBytes()}, null, 0));
        for (final String word : dictionary) {
            tree.insert(word);
        }
        final List<Object[]> tests = new ArrayList<>();
        for (final String query : dictionary) {
            final Map<Integer, Set<String>> expectedResult= new HashMap<>();
            expectedResult.put(0, new HashSet<>(Collections.singletonList(query)));
            tests.add(new Object[]{query, 0, tree, new HashMap<>(expectedResult)});
            if (query.equals("GCTCG")) {
                expectedResult.put(1, new HashSet<>(Collections.singletonList("GGTCG")));
            } else if (query.equals("CGTCG")) {
                expectedResult.put(1, new HashSet<>(Arrays.asList("AGTCG", "GGTCG")));
            } else if (query.equals("AGTCG")) {
                expectedResult.put(1, new HashSet<>(Arrays.asList("CGTCG", "GGTCG")));
            } else if (query.equals("GGTCG")) {
                expectedResult.put(1, new HashSet<>(Arrays.asList("CGTCG", "AGTCG", "GCTCG")));
            }
            tests.add(new Object[]{query, 1, tree, new HashMap<>(expectedResult)});
        }
        return tests.toArray(new Object[][]{});
    }
    @Test(dataProvider = "testQueryDataProvider")
    public void testQuery(final String query, final int maxDist, final BKTree<String> tree, final Map<Integer, Set<String>> expectedResult) {
        final Map<Integer, LinkedList<String>> results = tree.query(query, maxDist);
        Assert.assertEquals(results.size(), expectedResult.size());
        for (final Map.Entry<Integer, LinkedList<String>> matchSet : results.entrySet()) {
            Assert.assertEquals(matchSet.getValue().size(), expectedResult.get(matchSet.getKey()).size());
            Assert.assertTrue(expectedResult.get(matchSet.getKey()).containsAll(matchSet.getValue()));
        }
    }
}