package picard.illumina;

import htsjdk.samtools.cram.common.IntHashMap;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.OtherProgramGroup;
import picard.sam.util.Pair;
import picard.util.StringDistanceUtils;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Program to choose barcodes from a list of candidates and a distance requirement
 */
@CommandLineProgramProperties(
        summary = "...",
        oneLineSummary = "Program to choose barcodes from a list of candidates and a distance requirement",
        programGroup = OtherProgramGroup.class
)
@DocumentedFeature
public class SelectBarcodes extends CommandLineProgram {

    @Argument
    public List<File> BARCODES_MUST_HAVE;

    @Argument
    public List<File> BARCODES_CHOOSE_FROM;

    @Argument(optional = true)
    public File SEED_BARCODES;

    @Argument()
    public Integer DISTANCE_TO_SEEDS = 2;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Argument
    public int FAR_ENOUGH = 3;

    @Argument
    public boolean COMPARE_REVCOMP = false;

    @Argument
    public boolean ALLOW_REVCOMP = false;

    @Argument
    public boolean ALLOW_REV = false;

    @Argument
    public boolean ALLOW_COMP = false;

    @Argument(optional = true)
    public File DISTANCES = null;

    private static final List<String> mustHaveBarcodes = new ArrayList<>();
    static final List<String> barcodes = new ArrayList<>();
    private final List<String> seedBarcodes = new ArrayList<>();

    private final Map<Integer, BitSet> adjacencyMatrix = new HashMap<>();
    private final Map<Integer, BitSet> seedAdjacencyMatrix = new HashMap<>();

    static File output = null;
    private static Log LOG = Log.getInstance(SelectBarcodes.class);

    // The P verticies are the "Potential" verticies, verticies that may be chosen
    private static final List<BitSet> Ps = new DefaultingArrayList<>(i->new BitSet());

    // The X verticies are the "eXcluded" verticies, those that should not be examined
    private static final List<BitSet> Xs = new DefaultingArrayList<>(i->new BitSet());

    // the R verticies are the "Required" verticies, those that must be part of the clique
    private static final List<BitSet> Rs = new DefaultingArrayList<>(i->new BitSet());

    // this is a temporary variable holding the difference between P and the neighbors of the "pivot" vertex
    private static final Map<Integer, BitSet> Diffs = new HashMap<>();
    private static int recursionLevel;


    @Override
    public int doWork() {
        output = OUTPUT;

        //open files make lists of barcodes
        openBarcodes();

        LOG.info("Opened Barcode files. ");
        // calculate distance matrix and adjacency matrix

        calculateAdjacencyMatrix();
        LOG.info("Calculated distances");

        LOG.info("there are " + mustHaveBarcodes.size() + " MUST_HAVE barcodes.");
        LOG.info("there are " + barcodes.size() + " other barcodes to choose from (after possibly rejecting some).");
        LOG.info("there are " + seedBarcodes.size() + " seed barcodes");

        try (final PrintWriter writer = new PrintWriter("all.barcodes.txt")) {
            mustHaveBarcodes.forEach(writer::println);
            barcodes.forEach(writer::println);

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        final BitSet R = new BitSet();
        final BitSet X = new BitSet();

        final BitSet nodeSubset = new BitSet();

        // add each group nodes from the seeds and
        // call BronKerbosch2 saving best selections to disk.
        for (int i = 0; i < seedBarcodes.size(); i++) {
            LOG.info("Adding " + seedAdjacencyMatrix.get(i).cardinality() + " nodes from seed " + i);
            nodeSubset.or(seedAdjacencyMatrix.get(i));

            final BitSet seedSolution = find_cliques(subsetGraph(adjacencyMatrix, nodeSubset), R);
            // add new solution to the required nodes
            R.or(seedSolution);

            // add other nodes to excluded nodes
            X.or(difference(seedSolution, nodeSubset));
        }
        // finally, add all remaining nodes

        LOG.info("Adding " + (adjacencyMatrix.size() - nodeSubset.cardinality()) + " remaining nodes.");

        final BitSet solution = find_cliques(adjacencyMatrix, R);

        LOG.info("final solution has cardinality " + solution.cardinality() + mustHaveBarcodes.size());

        return 0;
    }

    private void calculateAdjacencyMatrix() {

        final List<String> filteredBarcodes = barcodes.stream().filter(b -> {
                    final Optional<String> firstClose = mustHaveBarcodes
                            .stream()
                            .filter(m -> !areFarEnough(b, m))
                            .findAny();
                    if (firstClose.isPresent()) {

                        LOG.info(String.format("rejecting barcode: %s, it's too close to a MUST_HAVE barcode: %s.",
                                b, firstClose.get()));
                        return false;
                    }
                    return true;
                }
        ).collect(Collectors.toList());

        barcodes.clear();
        barcodes.addAll(filteredBarcodes);

        for (int ii = 0; ii < barcodes.size(); ii++) {
            final BitSet adjacency = new BitSet(barcodes.size());

            for (int jj = 0; jj < barcodes.size(); jj++) {
                adjacency.set(jj, areFarEnough(barcodes.get(ii), barcodes.get(jj)));
            }

            adjacencyMatrix.put(ii, adjacency);
        }

        if (DISTANCES != null) {
            try (final PrintWriter writer = new PrintWriter(DISTANCES)) {
                writer.append("BARCODE\t");
                writer.println(String.join("\t", barcodes));
                for (int ii = 0; ii < barcodes.size(); ii++) {
                    final BitSet adjacency = adjacencyMatrix.get(ii);

                    writer.append(barcodes.get(ii)).append('\t');
                    for (int jj = 0; jj < barcodes.size(); jj++) {
                        writer.append(adjacency.get(jj) ? "1" : "0").append('\t');
                    }
                    writer.append('\n');
                }
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }

        for (int ii = 0; ii < seedBarcodes.size(); ii++) {
            final BitSet adjacency = new BitSet(barcodes.size());

            for (int jj = 0; jj < barcodes.size(); jj++) {
                adjacency.set(jj, levenshtein(seedBarcodes.get(ii), barcodes.get(jj), DISTANCE_TO_SEEDS + 1) <= DISTANCE_TO_SEEDS);
            }

            seedAdjacencyMatrix.put(ii, adjacency);
            LOG.info("Seed " + ii + " has " + adjacency.cardinality() + " barcodes near it.");
        }
    }

    static int levenshtein(final String lhs, final String rhs, final int farEnough) {
        return StringDistanceUtils.levenshteinDistance(lhs.getBytes(), rhs.getBytes(), farEnough);
    }

    private boolean areFarEnough(final String lhs, final String rhs) {
        if (!ALLOW_REV) {
            final byte[] rev = rhs.getBytes();
            SequenceUtil.reverse(rev, 0, rev.length);
            if (lhs.equals(Arrays.toString(rev))) {
                return false;
            }
        }

        if (!ALLOW_COMP) {
            final byte[] comp = SequenceUtil.reverseComplement(rhs).getBytes();
            SequenceUtil.reverse(comp, 0, comp.length);
            if (lhs.equals(Arrays.toString(comp))) {
                return false;
            }
        }

        if (!ALLOW_REVCOMP && lhs.equals(SequenceUtil.reverseComplement(rhs))) {
            return false;
        }

        return levenshtein(lhs, rhs, FAR_ENOUGH) >= FAR_ENOUGH &&
                (!COMPARE_REVCOMP || levenshtein(lhs, SequenceUtil.reverseComplement(rhs), FAR_ENOUGH) >= FAR_ENOUGH);
    }

    //returns the number of SEED barcodes that were placed into "barcodes"
    private void openBarcodes() {
        barcodes.clear();
        mustHaveBarcodes.clear();

        BARCODES_MUST_HAVE.forEach(b -> readBarcodesFile(b, mustHaveBarcodes));
        LOG.info("There are " + mustHaveBarcodes.size() + " required barcodes.");
        BARCODES_CHOOSE_FROM.forEach(b -> readBarcodesFile(b, barcodes));
        LOG.info("There are " + barcodes.size() + " (original) barcodes to chose from.");

        if (SEED_BARCODES != null) {
            readBarcodesFile(SEED_BARCODES, seedBarcodes);
            LOG.info("There are " + seedBarcodes.size() + " seed barcodes.");
        }

        //shuffle input barcodes to prevent bias
        Collections.shuffle(barcodes, new Random(51));
    }

    private void readBarcodesFile(final File f, final List<String> addTo) {
        try (BufferedReader br = new BufferedReader(new FileReader(f))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.trim().equals("")) continue;
                if (!addTo.contains(line)) {
                    addTo.add(line);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    /**
     * implements BronKerbosch3 https://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm
     * to find the maximal clique in the graph defined by the edges in {@code graph}. The
     * original algorithm is defined recursively, but here the recursion has been unfolded explicitly
     * to avoid the recursion overhead and to reuse objects.
     * <p>
     * BronKerbosch3(G):
     * 1       P = V(G)
     * 2       R = X = empty
     * 3       for each vertex v in a degeneracy ordering of G:
     * 4           BronKerbosch2({v}, P ⋂ N(v), X ⋂ N(v))
     * 5           P := P \ {v}
     * 6           X := X ⋃ {v}
     * <p>
     * <p>
     * BronKerbosch2(R,P,X):
     * 7       if P and X are both empty:
     * 8           report R as a maximal clique
     * 9       choose a pivot vertex u in P ⋃ X
     * 10      for each vertex v in P \ N(u):
     * 11          BronKerbosch2(R ⋃ {v}, P ⋂ N(v), X ⋂ N(v))
     * 12          P := P \ {v}
     * 13          X := X ⋃ {v}
     * <p>
     * in the code below comments refer to line numbers in this pseudo-code
     *
     * @param graph the input bit-sets defining the edges between barcodes that are compatible
     */
    static BitSet find_cliques(final Map<Integer, BitSet> graph, final BitSet required) {

        final BitSet excluded = new BitSet();
        // any node that isn't a neighbor of all the required nodes should be removed from consideration
        graph.keySet().stream()
                .filter(i -> intersection(graph.get(i), required).cardinality() != required.cardinality())
                .forEach(excluded::set);


        final BitSet presentInGraph = new BitSet();
        graph.keySet().forEach(presentInGraph::set);

        final Integer[] degeneracyOrder = getDegeneracyOrder(subsetGraph(graph, difference(presentInGraph, excluded)));

        final BitSet best_clique = new BitSet();
        int bestCliqueSize = 0;

        //1
        final BitSet pTop = new BitSet();
        graph.keySet().forEach(pTop::set);

        //2 (modified to allow input of required nodes)
        final BitSet rTop = new BitSet();
        rTop.clear();
        rTop.or(required);

        //2 (modified to allow input of required nodes, which, in turn causes some node to be excluded)
        final BitSet xTop = new BitSet();
        xTop.clear();
        xTop.or(excluded);

        //3
        for (final Integer v : degeneracyOrder) {
            recursionLevel = 0;

            LOG.info("examining node " + v);

            //4 {v}
            BitSet r = Rs.get(recursionLevel);
            r.clear();
            r.or(rTop);
            if (r.get(v)) {
                continue;
            }
            r.set(v);

            //4 P ⋂ N(v)
            BitSet p = Ps.get(recursionLevel);
            p.clear();
            p.or(pTop);
            p.and(graph.get(v));

            // 4 X ⋂ N(v)
            BitSet x = Xs.get(recursionLevel);
            x.clear();
            x.or(xTop);
            if (x.get(v)) {
                continue;
            }
            x.and(graph.get(v));

            while (recursionLevel >= 0) {

                p = Ps.get(recursionLevel);
                r = Rs.get(recursionLevel);
                x = Xs.get(recursionLevel);

                //  LOG.info(String.format("in while (%d) %s %s %s (%s)", recursionLevel, r, p, x, best_clique));

                // if there are no more pivots to chose from, we have reached a (locally) maximal clique.
                // or if there is no way we could find a larger clique, stop trying

                // 7
                if (p.isEmpty() || p.cardinality() + r.cardinality() <= bestCliqueSize) {

                    Diffs.remove(recursionLevel);
                    registerClique(r, best_clique);
                    bestCliqueSize = best_clique.cardinality();
                    recursionLevel--;

                    continue;
                }

                if (!Diffs.containsKey(recursionLevel))  {
                    BitSet finalP = p;

                    // 9 choosing a pivot
                    final BitSet pOrX = union(p, x);
                    int u = pOrX.stream()
                            .mapToObj(o-> new Pair<>(o, intersection(finalP, graph.get(o)).cardinality()))
                            .max(Comparator.comparingInt(Pair::getRight))
                            .get()
                            .getLeft();

                    //10
                    Diffs.put(recursionLevel, difference(p, graph.get(u)));
                }

                final int vv = Diffs.get(recursionLevel).nextSetBit(0);
                if (vv == -1) {
                    Diffs.remove(recursionLevel);
                    p.clear();
                    x.clear();
                    r.clear();

                    recursionLevel--;
                    continue;
                }

                Diffs.get(recursionLevel).clear(vv);

                final BitSet recNeighs = graph.get(vv);

                //preparations for recursive call 11
                recursionLevel++;

                final BitSet lowerR = Rs.get(recursionLevel);
                lowerR.clear();
                // 11 R ⋃ {v}
                lowerR.or(r);
                lowerR.set(vv);

                final BitSet lowerX = Xs.get(recursionLevel);
                lowerX.clear();
                // 11 X ⋂ N(v)
                lowerX.or(x);
                lowerX.and(recNeighs);

                final BitSet lowerP = Ps.get(recursionLevel);
                lowerP.clear();
                // 11 P ⋂ N(v)
                lowerP.or(p);
                lowerP.and(recNeighs);

                // p and x are still connected to the previous recursion level
                // 12
                p.clear(vv);
                // 13
                x.set(vv);

            }

            // 5
            pTop.clear(v);
            // 6
            xTop.set(v);
        }
        return best_clique;
    }

    private static void registerClique(final BitSet r, final BitSet bestClique) {
        if (r.cardinality() > bestClique.cardinality()) {
            bestClique.clear();
            bestClique.or(r);
            LOG.info("best cardinality (so far) is " + (bestClique.cardinality() + mustHaveBarcodes.size()));

            output.delete();
            try (PrintWriter writer = new PrintWriter(output)) {

                mustHaveBarcodes.forEach(b -> writer.println(-1 + "\t" + b));
                bestClique.stream().forEach(i -> writer.println(i + "\t" + barcodes.get(i)));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }
    }

    static Map<Integer, BitSet> subsetGraph(final Map<Integer, BitSet> graph, final BitSet mask) {

        final Map<Integer, BitSet> retVal = new HashMap<>();

        mask.stream().forEach(i -> {
            retVal.put(i, intersection(graph.get(i), mask));
        });

        return retVal;
    }

    static BitSet union(final BitSet lhs, final BitSet rhs) {
        BitSet ret = BitSet.valueOf(lhs.toLongArray());
        ret.or(rhs);
        return ret;
    }

    static BitSet intersection(final BitSet lhs, final BitSet rhs) {
        BitSet ret = BitSet.valueOf(lhs.toLongArray());
        ret.and(rhs);
        return ret;
    }

    static BitSet difference(final BitSet lhs, final BitSet rhs) {
        BitSet ret = BitSet.valueOf(lhs.toLongArray());
        ret.andNot(rhs);
        return ret;
    }

    private static Integer[] getDegeneracyOrder(final Map<Integer, BitSet> graph) {
        final List<Integer> ordering = new ArrayList<>();
        final Set<Integer> ordering_set = new HashSet<>();

        // a map from the vertices to their cardinality
        final Map<Integer, Integer> degrees = new CollectionUtil.DefaultingMap<>(0);

        // a map form a given degeneracy to the list of vertices with that degeneracy
        @SuppressWarnings("MismatchedQueryAndUpdateOfCollection")
        Map<Integer, List<Integer>> degen = new CollectionUtil.DefaultingMap<>(i -> new ArrayList<>(), true);
        int max_deg = -1;
        for (int v : graph.keySet()) {
            final int deg = graph.get(v).cardinality();
            degen.get(deg).add(v);
            degrees.put(v, deg);
            if (deg > max_deg) {
                max_deg = deg;
            }
        }
        outter:
        while (true) {
            int i = 0;
            while (true) {
                if (i <= max_deg) {
                    if (degen.get(i).size() != 0) {
                        break;
                    }
                    i += 1;

                } else {
                    break outter;
                }
            }

            final Integer v = degen.get(i).remove(degen.get(i).size() - 1);

            ordering.add(v);
            ordering_set.add(v);
            graph.get(v).stream().forEach(w -> {

                if (!ordering_set.contains(w)) {
                    final int deg = degrees.get(w);
                    degen.get(deg).remove(Integer.valueOf(w));
                    if (deg > 0) {
                        degrees.put(w, degrees.get(w) - 1);
                        degen.get(deg - 1).add(w);
                    }
                }
            });
        }
        Collections.reverse(ordering);
        return ordering.toArray(new Integer[0]);
    }

    private static int argMax(final int[] input) {
        int max = input[0];
        int maxIdx = 0;
        for (int i = 1; i < input.length; i++) {
            if (input[i] > max) {
                max = input[i];
                maxIdx = i;
            }
        }
        return maxIdx;
    }

    private static  class DefaultingArrayList<TYPE> extends ArrayList<TYPE> {

        final private CollectionUtil.DefaultingMap.Factory<TYPE, Integer> generator;

        public DefaultingArrayList(CollectionUtil.DefaultingMap.Factory<TYPE, Integer> defaultGenerator) {
            super();
            this.generator = defaultGenerator;
        }

        @Override
        public TYPE get(int index) {
            while (index >= this.size()) {
                add(generator.make(this.size()));
            }
            return super.get(index);
        }

        @Override
        public TYPE set(int index, TYPE o) {
            while (index - 1 > this.size()) {
                add(generator.make(this.size()));
            }
            if (index > this.size()) {
                add(index, o);
                return super.get(index);
            } else {
                return super.set(index, o);
            }
        }
    }
}