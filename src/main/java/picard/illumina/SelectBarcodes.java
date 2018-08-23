package picard.illumina;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import org.apache.commons.lang3.StringUtils;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by farjoun on 8/21/18.
 */
@CommandLineProgramProperties(
        summary = "blah",
        oneLineSummary = "blah",
        programGroup = OtherProgramGroup.class
)
@DocumentedFeature
public class SelectBarcodes extends CommandLineProgram {

    @Argument
    public List<File> BARCODES_MUST_HAVE;

    @Argument
    public List<File> BARCODES_CHOOSE_FROM;

    @Argument
    public File OUTPUT;

    private static final int FAR_ENOUGH = 3;

    static final List<String> barcodes = new ArrayList<>();

    private final List<BitSet> ajacencyMatrix = new ArrayList<>();

    private int mustHaves = 0;

    static File output = null;
    static Log LOG = Log.getInstance(SelectBarcodes.class);


    @Override
    public int doWork() {
        output = OUTPUT;

        //open files make lists of barcodes
        openBarcodes();
        LOG.info("Opened Barcode files.");
        // calculate distance matrix and ajacency matrix
        calculateAjacencyMatrix();
        LOG.info("Calculated distances");
        // call BronKerbosch2 saving best selections to disk.
        find_cliques(ajacencyMatrix, mustHaves);

        return 0;
    }

    private void calculateAjacencyMatrix() {
        try (final PrintWriter writer = new PrintWriter("distances.txt")) {

            for (int ii = 0; ii < barcodes.size(); ii++) {
                BitSet ajacency = new BitSet(barcodes.size());
                for (int jj = 0; jj < barcodes.size(); jj++) {
                    ajacency.set(jj, areFarEnough(ii, jj));
                    writer.append(ajacency.get(jj) ? "1" : "0").append(String.valueOf('\t'));
                }
                writer.append('\n');
                ajacencyMatrix.add(ii, ajacency);
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    static int hamming(final String lhs, final String rhs) {
        return Math.min(picard.util.SequenceUtil.calculateEditDistance(lhs, rhs),
                picard.util.SequenceUtil.calculateEditDistance(lhs, SequenceUtil.reverseComplement(rhs)));

    }

    static int levenshtein(final String lhs, final String rhs) {
        return Math.min(StringUtils.getLevenshteinDistance(lhs, rhs),
                StringUtils.getLevenshteinDistance(lhs, SequenceUtil.reverseComplement(rhs)));

    }

    private boolean areFarEnough(final int lhs, final int rhs) {

        if (lhs <= mustHaves && rhs <= mustHaves) {
            return lhs != rhs;
        }
        if (lhs > mustHaves && rhs > mustHaves) {
            return levenshtein(barcodes.get(lhs), barcodes.get(rhs)) >= FAR_ENOUGH;
        } else {
            return hamming(barcodes.get(lhs), barcodes.get(rhs)) >= FAR_ENOUGH;
        }
    }

    private void openBarcodes() {
        barcodes.clear();

        BARCODES_MUST_HAVE.forEach(this::readBarcodesFile);

        mustHaves = barcodes.size();

        BARCODES_CHOOSE_FROM.forEach(this::readBarcodesFile);

        try (final PrintWriter writer = new PrintWriter("all.barcodes.txt")) {
            barcodes.forEach(writer::println);

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    private void readBarcodesFile(File f) {
        try (BufferedReader br = new BufferedReader(new FileReader(f))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.trim().equals("")) continue;
                if (!barcodes.contains(line) && !barcodes.contains(SequenceUtil.reverseComplement(line))) {
                    barcodes.add(line);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    static void find_cliques(List<BitSet> graph, int mustHaves) {

        final BitSet r = new BitSet(barcodes.size());
        r.set(0, mustHaves, true);

        final BitSet p = new BitSet(barcodes.size());
        p.set(0, barcodes.size(), true);
        p.andNot(r);

        final BitSet x = new BitSet(barcodes.size());

        final BitSet best_clique = new BitSet();

        for (final Integer v : getDegeneracyOrder(graph)) {

            if (r.get(v)) continue;
            LOG.info("examining node " + v);
            final BitSet neighs = graph.get(v);
            find_cliques_pivot(graph, union(r, v), intersection(p, neighs), intersection(x, neighs), best_clique);
            p.clear(v);
            x.set(v);
        }
    }

    static BitSet union(final BitSet set, final int node) {
        BitSet ret = BitSet.valueOf(set.toLongArray());
        ret.set(node);
        return ret;
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

    private static Integer[] getDegeneracyOrder(final List<BitSet> graph) {
        List<Integer> ordering = new ArrayList<>();
        Set<Integer> ordering_set = new HashSet<>();

        // a map from the vertices to their cardinality
        Map<Integer, Integer> degrees = new CollectionUtil.DefaultingMap<>(0);

        // a map form a given degeneracy to the list of vertices with that degeneracy
        @SuppressWarnings("MismatchedQueryAndUpdateOfCollection")
        Map<Integer, List<Integer>> degen = new CollectionUtil.DefaultingMap<>(i -> new ArrayList<>(), true);
        int max_deg = -1;
        for (int v = 0; v < graph.size(); v++) {
            int deg = graph.get(v).cardinality();
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
        return ordering.toArray(new Integer[0]);

    }

    static void find_cliques_pivot(final List<BitSet> graph, final BitSet r, final BitSet p, final BitSet x, final BitSet bestClique) {

        if (p.isEmpty() && x.isEmpty()) {

            registerClique(r, bestClique);
            return;
        }
        final int u = union(p, x).stream().findFirst().orElse(-1); //should never happen
        for (int v : difference(p, graph.get(u)).stream().toArray()) {
            final BitSet neighs = graph.get(v);
            find_cliques_pivot(graph, union(r, v), intersection(p, neighs), intersection(x, neighs), bestClique);
            p.clear(v);
            x.set(v);
        }
    }

    static private void registerClique(final BitSet r, final BitSet bestClique) {
        if (r.cardinality() > bestClique.cardinality()) {
            bestClique.clear();
            bestClique.or(r);
            System.out.println("best.cardinality()=" + bestClique.cardinality());

            try (PrintWriter writer = new PrintWriter(output)) {
                bestClique.stream().forEach(i -> writer.println(i + '\t' + barcodes.get(i)));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }
    }
}