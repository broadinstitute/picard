package picard.illumina;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import org.apache.commons.lang3.StringUtils;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.*;
import java.util.*;

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

    static final private int FAR_ENOUGH = 3;

    private final List<String> barcodes = new ArrayList<>();

    private List<BitSet> ajacencyMatrix = new ArrayList<>();

    private int mustHaves = 0;

    @Override
    public int doWork() {

        //open files make lists of barcodes
        openBarcodes();
        // calculate distance matrix and ajacency matrix
        calculateAjacencyMatrix();
        // call BronKerbosch2 saving best selections to disk.
        find_cliques();
        return 0;
    }

    private void calculateAjacencyMatrix() {
        for (int ii = 0; ii < barcodes.size(); ii++) {
            BitSet ajacency = new BitSet(barcodes.size());
            for (int jj = 0; jj < barcodes.size(); jj++) {
                ajacency.set(jj, areFarEnough(ii, jj));
            }
            ajacencyMatrix.add(ii, ajacency);
        }
    }

    private boolean areFarEnough(final int lhs, final int rhs) {

        final int distance1;
        final int distance2;
        if (lhs > mustHaves && rhs > mustHaves) {
            distance1 = StringUtils.getLevenshteinDistance(barcodes.get(lhs), barcodes.get(lhs), FAR_ENOUGH);
            distance2 = StringUtils.getLevenshteinDistance(barcodes.get(lhs), SequenceUtil.reverseComplement(barcodes.get(lhs)), FAR_ENOUGH);
        } else {
            distance1 = StringUtils.getLevenshteinDistance(barcodes.get(rhs), barcodes.get(lhs), FAR_ENOUGH);
            distance2 = StringUtils.getLevenshteinDistance(barcodes.get(rhs), SequenceUtil.reverseComplement(barcodes.get(lhs)), FAR_ENOUGH);
        }
        return Math.min(distance1, distance2) < FAR_ENOUGH;
    }

    private void openBarcodes() {
        BARCODES_MUST_HAVE.forEach(f -> {
            try (BufferedReader br = new BufferedReader(new FileReader(f))) {
                String line;
                while ((line = br.readLine()) != null) {
                    if (!barcodes.contains(line)) {
                        barcodes.add(line);
                    }
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        });

        mustHaves = barcodes.size();

        BARCODES_CHOOSE_FROM.forEach(f -> {
            try (BufferedReader br = new BufferedReader(new FileReader(f))) {
                String line;
                while ((line = br.readLine()) != null) {
                    if (!barcodes.contains(line))
                        barcodes.add(line);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        });
    }

    private void find_cliques() {

        final BitSet r = new BitSet(barcodes.size());
        r.set(0, mustHaves, true);

        final BitSet p = new BitSet(barcodes.size());
        p.set(0, barcodes.size(), true);
        p.andNot(r);

        final BitSet x = new BitSet(barcodes.size());

        final BitSet best_clique = new BitSet();

        for (final Integer v : getDegeneracyOrder()) {
            final BitSet neighs = ajacencyMatrix.get(v);
            find_cliques_pivot(union(r, v), intersection(p, neighs), intersection(x, neighs), best_clique);
            p.clear(v);
            x.set(v);
        }
    }

    private BitSet union(final BitSet set, final int node) {
        BitSet ret = new BitSet();
        ret.or(set);
        ret.set(node);
        return ret;
    }

    private BitSet union(final BitSet lhs, final BitSet rhs) {
        BitSet ret = new BitSet();
        ret.or(lhs);
        ret.or(rhs);
        return ret;
    }

    private BitSet intersection(final BitSet lhs, final BitSet rhs) {
        BitSet ret = new BitSet();
        ret.or(lhs);
        ret.and(rhs);
        return ret;
    }

    private BitSet difference(final BitSet lhs, final BitSet rhs) {
        BitSet ret = new BitSet();
        ret.or(lhs);
        ret.andNot(rhs);
        return ret;
    }

    private Integer[] getDegeneracyOrder() {
        List<Integer> ordering = new ArrayList<>();
        Set<Integer> ordering_set = new HashSet<>();

        // a map from the vertices to their cardinality
        Map<Integer, Integer> degrees = new CollectionUtil.DefaultingMap<>(0);

        // a map form a given degeneracy to the list of verticies with that degeneracy
        Map<Integer, Set<Integer>> degen = new CollectionUtil.DefaultingMap<>(i -> new HashSet<>(), true);
        int max_deg = -1;
        for (int v = 0; v < barcodes.size(); v++) {
            int deg = ajacencyMatrix.get(v).cardinality();
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

            final Integer v = degen.get(i).iterator().next();
            degen.get(i).remove(v);

            ordering.add(v);
            ordering_set.add(v);
            ajacencyMatrix.get(v).stream().forEach(w -> {

                if (!ordering_set.contains(w)) {
                    final int deg = degrees.get(w);
                    degen.get(deg).remove(w);
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

    private void find_cliques_pivot(final BitSet r, final BitSet p, final BitSet x, final BitSet bestClique) {
        if (p.isEmpty() && x.isEmpty()) {
            registerClique(r, bestClique);
            return;
        }
        final int u = union(p, x).stream().findFirst().orElse(-1); //should never happen
        for (int v : difference(p, ajacencyMatrix.get(u)).stream().toArray()) {
            final BitSet neighs = ajacencyMatrix.get(v);
            find_cliques_pivot(union(r, v), intersection(p, neighs), intersection(x, neighs), bestClique);
            p.clear(v);
            x.set(v);
        }
    }

    private void registerClique(final BitSet r, final BitSet bestClique) {
        if (r.cardinality() > bestClique.cardinality()) {
            bestClique.clear();
            bestClique.or(r);
            System.out.println("best.cardinality()=" + bestClique.cardinality());

            try (PrintWriter writer = new PrintWriter(OUTPUT)) {
                bestClique.stream().forEach(i -> writer.println(barcodes.get(i)));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        } else {
            System.out.println("backtrack:" + r.cardinality());
        }
    }
}