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
}
