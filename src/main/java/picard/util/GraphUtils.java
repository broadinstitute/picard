/*
* The MIT License
*
* Copyright (c) 2016 The Broad Institute
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
* THE SOFTWARE.
*/

package picard.util;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by farjoun on 11/2/16.
 */
public class GraphUtils {

    static public class Graph<Node extends Comparable<Node>> {
        final private List<Node> nodes;
        final private List<List<Integer>> neighbors;

        /**
         * returns the cluster map of connected components
         *
         * @return Nodes that point to the same integer are in the same cluster.
         */
        public Map<Node, Integer> cluster() {

            final int[] cluster = IntStream.range(0, nodes.size()).toArray();

            IntStream.range(0, neighbors.size()).forEach(i ->
                    neighbors.get(i).stream().forEach(j -> joinNodes(cluster, j, i)));

            return nodes.stream().collect(Collectors.toMap(n -> n, n -> cluster[nodes.indexOf(n)]));

        }

        public Graph() {
            nodes = new ArrayList<>();
            neighbors = new ArrayList<>();
        }

        /* directed and private */
        private void addNeighbor(final Integer fromNode, final Integer toNode) {
            final List<Integer> fromNodesNeighbors = neighbors.get(fromNode);

            if (!fromNodesNeighbors.contains(toNode)) {
                fromNodesNeighbors.add(toNode);
            }
        }

        public Integer addNode(final Node singleton) {

            if (!nodes.contains(singleton)) {
                nodes.add(singleton);
                neighbors.add(new ArrayList<>());
            }
            return nodes.indexOf(singleton);
        }

        /* bidirectional and public */
        public void addEdge(Node left, Node right) {

            final int leftIndex = addNode(left);
            if (left == right) return;

            final int rightIndex = addNode(right);

            addNeighbor(leftIndex, rightIndex);
            addNeighbor(rightIndex, leftIndex);
        }

        // Part of Union-Find with Path Compression that joins to nodes to be part of the same cluster.
        private static void joinNodes(int[] grouping, final int nodeId1, final int nodeId2) {
            final int repNode1 = findRepNode(grouping, nodeId1);
            final int repNode2 = findRepNode(grouping, nodeId2);
            if (repNode1 == repNode2) return;
            grouping[repNode1] = repNode2;
        }

        // Part of Union-Find with Path Compression to determine the duplicate set a particular UMI belongs to.
        private static int findRepNode(final int[] grouping, int nodeId) {
            int representativeUmi = nodeId; // All UMIs of a duplicate set will have the same reprsentativeUmi.
            while (representativeUmi != grouping[representativeUmi]) {
                representativeUmi = grouping[representativeUmi];
            }
            while (nodeId != representativeUmi) {
                int newUmiID = grouping[nodeId];
                grouping[nodeId] = representativeUmi;
                nodeId = newUmiID;
            }
            return representativeUmi;
        }
    }
}
