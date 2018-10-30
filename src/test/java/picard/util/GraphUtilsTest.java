package picard.util;

;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.Map;
import java.util.stream.IntStream;

import static org.testng.Assert.*;

/**
 * Created by farjoun on 11/2/16.
 */
public class GraphUtilsTest {

    @Test
    public void simpleTest() {
        GraphUtils.Graph<Integer> graph = new GraphUtils.Graph<>();

        graph.addEdge(5, 6);
        graph.addEdge(5, 8);
        graph.addEdge(7, 9);
        graph.addEdge(7, 8);

        graph.addEdge(4, 2);
        graph.addEdge(1, 2);
        graph.addEdge(3, 1);

        Map<Integer, Integer> clusters = graph.cluster();

        // 9 nodes
        assertEquals(clusters.size(), 9);

        Integer fivesCluster = clusters.get(5);
        Arrays.asList(5, 6, 7, 8, 9).stream().forEach(i -> assertEquals(clusters.get(i), fivesCluster));

        Integer foursCluster = clusters.get(4);
        Arrays.asList(1, 2, 3, 4).stream().forEach(i -> assertEquals(clusters.get(i), foursCluster));

        assertNotEquals(fivesCluster, foursCluster);
    }

    @Test
    public void secondTest() {
        GraphUtils.Graph<Integer> graph = new GraphUtils.Graph<>();

        graph.addEdge(1, 3);
        graph.addEdge(2, 3);
        graph.addEdge(0, 2);
        graph.addEdge(4, 3);

        graph.addEdge(5, 6);
        graph.addEdge(7, 6);
        graph.addNode(8);

        Map<Integer, Integer> clusters = graph.cluster();

        // 9 nodes
        assertEquals(clusters.size(), 9);

        Integer threesCluster = clusters.get(3);
        IntStream.range(0,5).forEach(i -> assertEquals(clusters.get(i), threesCluster));

        Integer fivesCluster = clusters.get(5);
        IntStream.range(5,8).forEach(i -> assertEquals(clusters.get(i), fivesCluster));

        Integer eightsCluster = clusters.get(8);
        Collections.singletonList(8).stream().forEach(i -> assertEquals(clusters.get(i), eightsCluster));

        assertNotEquals(fivesCluster, threesCluster);
        assertNotEquals(fivesCluster, eightsCluster);
        assertNotEquals(eightsCluster, threesCluster);
    }


    @Test
    public void pathologicalInputOrderingTest() {
        final GraphUtils.Graph<Integer> graph = new GraphUtils.Graph<>();

        final int a = 0, b = 1, c = 2, d = 3, e = 4, f = 5, h = 6, g = 7;
        graph.addNode(a);
        graph.addNode(b);
        graph.addNode(c);
        graph.addNode(d);
        graph.addNode(f);
        graph.addNode(h);
        graph.addNode(e);
        graph.addNode(g);

        // This graph is a hand constructed case where the algorithm produces in its final form an orphaned node that only
        // transitively points to the root for the rest of the graph (node d as it turns out)
        graph.addEdge(a, h);
        graph.addEdge(b, f);
        graph.addEdge(c, h);
        graph.addEdge(d, f);
        graph.addEdge(e, h);
        graph.addEdge(e, f);
        graph.addEdge(g, a);

        final Map<Integer, Integer> clusters = graph.cluster();

        // 9 nodes
        assertEquals(clusters.size(), 8);

        // Asserting that this connected graph has every node pointing to the same place
        final Integer repCluster = clusters.get(a);
        IntStream.range(0,7).forEach(i -> assertEquals(clusters.get(i), repCluster));
    }
}
