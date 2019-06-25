package picard.util;

import org.apache.commons.lang3.tuple.Pair;

import java.io.Serializable;
import java.util.*;
import java.util.function.BiFunction;

public class BKTree<T> implements Serializable {
    private T root;
    final private BiFunction<T, T, Integer> distanceMetric;
    final private HashMap<Integer, BKTree<T>> children;

    /**
     * Constructor
     *
     * @param root_value value to initialize as root of tree
     * @param metric distance metric to use for tree.  MUST SATISFY TRIANGLE INEQUALITY, ie d(A,B) + d(B,C) >= d(A,C).  If the function used
     *               does not satisfy the triangle inequality, then results of using BKTree are unreliable, and can easily be incorrect
     */
    public BKTree(final T root_value, final BiFunction<T, T, Integer> metric) {
        root = root_value;
        distanceMetric = metric;
        children = new HashMap<>();
    }

    /**
     *
     * @param metric distance metric to use for tree.  MUST SATISFY TRIANGLE INEQUALITY, ie d(A,B) + d(B,C) >= d(A,C).  If the function used
     *               does not satisfy the triangle inequality, then results of using BKTree are unreliable, and can easily be incorrect
     */
    public BKTree(final BiFunction<T, T, Integer> metric) {
        distanceMetric = metric;
        children = new HashMap<>();
    }

    /**
     * Add entry to tree
     * @param newEntry
     */
    public void insert(final T newEntry) {
        if (root != null) {
            final int distance = distanceMetric.apply(root, newEntry);
            if (children.containsKey(distance)) {
                children.get(distance).insert(newEntry);
            } else {
                children.put(distance, new BKTree<>(newEntry, distanceMetric));
            }
        }
        else {
            root = newEntry;
        }
    }

    /**
     * find all entries within maxDist of query.
     * @param query
     * @param maxDist
     * @return Map from distance to list of entries within that distance of query
     */
    public HashMap<Integer, LinkedList<T>> query(final T query, final int maxDist) {
        final HashMap<Integer, LinkedList<T>> ret = new HashMap<>();

        final int distance = distanceMetric.apply(root, query);
        if (distance <= maxDist) {
            ret.put(distance, new LinkedList<>(Collections.singletonList(root)));
        }

        final int lowLimit = distance - maxDist;
        final int highLimit = maxDist + distance;

        for (int i = lowLimit; i<= highLimit; i++) {
            if (children.containsKey(i)) {
                final HashMap<Integer, LinkedList<T>> retChildren = children.get(i).query(query, maxDist);
                for (HashMap.Entry<Integer, LinkedList<T>> entry : retChildren.entrySet()) {
                    if (ret.containsKey(entry.getKey())) {
                        ret.get(entry.getKey()).addAll(entry.getValue());
                    } else {
                        ret.put(entry.getKey(), entry.getValue());
                    }
                }
            }
        }
        return ret;
    }

    /**
     * get entries with smallest distance, among entries within maxDist
     * @param query
     * @param maxDist
     * @return list of entries with minimum distance <=maxDist.  returns empty list if no entries withing maxDist of query
     */
    public LinkedList<T> queryBest(final T query, final int maxDist) {
        HashMap<Integer, LinkedList<T>> queryRet = query(query, maxDist);
        if (queryRet.isEmpty()) {
            return new LinkedList<>();
        }

        int minDist = Collections.min(queryRet.keySet());
        return queryRet.get(minDist);
    }

    /**
     * get the first entry found within maxDist
     * @param query
     * @param maxDist
     * @return pair of found and entry and distance, or null if no entry found within maxDist
     */
    public Pair<T, Integer> queryFirst(final T query, final int maxDist) {
        final int distance = distanceMetric.apply(root, query);
        if (distance <= maxDist) {
            return Pair.of(root, distance);
        }

        final int lowLimit = distance - maxDist;
        final int highLimit = maxDist + distance;
        for (int i = lowLimit; i<= highLimit; i++)
            if (children.containsKey(i)) {
                final Pair<T, Integer> ret = children.get(i).queryFirst(query, maxDist);
                if (ret != null) {
                    return ret;
                }
            }
        return null;
    }
}
