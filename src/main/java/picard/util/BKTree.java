package picard.util;

import java.io.Serializable;
import java.util.*;
import java.util.function.BiFunction;

public class BKTree<T> implements Serializable {
    private T root;
    private BiFunction<T, T, Integer> distanceMetric;
    private HashMap<Integer, BKTree<T>> children;

    public BKTree(final T root_value, final BiFunction<T, T, Integer> metric) {
        root = root_value;
        distanceMetric = metric;
        children = new HashMap<>();
    }

    public BKTree(final BiFunction<T, T, Integer> metric) {
        distanceMetric = metric;
        children = new HashMap<>();
    }

    public void insert(final T newEntry) {
        if (root != null) {
            final int distance = distanceMetric.apply(newEntry, root);
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

    public HashMap<Integer, LinkedList<T>> query(final T query, final int maxDist) {
        final HashMap<Integer, LinkedList<T>> ret = new HashMap<>();

        final int distance = distanceMetric.apply(query, root);
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

    public LinkedList<T> queryBest(final T query, final int maxDist) {
        HashMap<Integer, LinkedList<T>> queryRet = query(query, maxDist);
        if (queryRet.isEmpty()) {
            return new LinkedList<>();
        }

        int minDist = Collections.min(queryRet.keySet());
        return queryRet.get(minDist);
    }
}
