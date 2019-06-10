package picard.util;
import java.io.Serializable;
import java.util.Arrays;
import java.util.HashMap;

import htsjdk.samtools.util.SequenceUtil;

public class BKNode implements Serializable {
    final private String value;
    final private HashMap<Integer, BKNode> children;

    BKNode(final String value) {
        this.value = value;
        children = new HashMap<>();
    }

    void insert(final String newValue) {
        final int distance = levenshteinDistance(value.getBytes(), newValue.getBytes(), value.length());
        if (children.containsKey(distance)) {
            children.get(distance).insert(newValue);
        }
        else {
            children.put(distance, new BKNode(newValue));
        }
    }

    HashMap<String, Integer> query(final String queryString, final int maxDist) {
        final HashMap<String, Integer> ret = new HashMap<>();

        final int distance = levenshteinDistance(value.getBytes(), queryString.getBytes(), value.length());
        if (distance <= maxDist) {
            ret.put(value, distance);
        }

        final int lowLimit = distance - maxDist;
        final int highLimit = maxDist + distance;

        for (int i = lowLimit; i<= highLimit; i++) {
            if (children.containsKey(i)) {
                ret.putAll(children.get(i).query(queryString, maxDist));
            }
        }
//        children.entrySet().stream().filter(e -> e.getKey() >= lowLimit && e.getKey() <= highLimit)
//                .forEach( e -> ret.putAll(e.getValue().query(queryString, maxDist)));

        return ret;
    }

    public static int levenshteinDistance(final byte[] barcode, final byte[] read, final int threshold) {
        int n = barcode.length; // length of left
        int m = read.length; // length of right

        if (n != m) {
            throw new IllegalArgumentException("This version of levenshteinDistance is speficially made for comparing strings " +
                    "of equal length. found " + n + " and " + m + ".");
        }
        // if one string is empty, the edit distance is necessarily the length
        // of the other
        if (n == 0) {
            return 0;
        }

        // it's easier to ignore indels in the beginning than in the end...so we copy and reverse the arrays
        final byte[] barcodeRev = Arrays.copyOf(barcode, barcode.length);
        final byte[] readRev = Arrays.copyOf(read, read.length);

        SequenceUtil.reverse(barcodeRev, 0, barcodeRev.length);
        SequenceUtil.reverse(readRev, 0, readRev.length);

        int[] previousCost = new int[n + 1]; // 'previous' cost array, horizontally
        int[] cost = new int[n + 1]; // cost array, horizontally
        int[] tempD; // placeholder to assist in swapping previousCost and cost

        // fill in starting table values
        final int boundary = Math.min(n, threshold) + 1;
        for (int i = 0; i < boundary; i++) {
            previousCost[i] = 0;// indels in the beginning do not cost
        }
        // these fills ensure that the value above the rightmost entry of our
        // stripe will be ignored in following loop iterations
        Arrays.fill(previousCost, boundary, previousCost.length, Integer.MAX_VALUE);
        Arrays.fill(cost, Integer.MAX_VALUE);

        // iterates through t
        for (int j = 1; j <= m; j++) {
            final byte readJ = readRev[j - 1]; // jth character of readRev
            cost[0] = 0;

            // compute stripe indices, constrain to array size
            final int min = Math.max(1, j - threshold);
            final int max = j > Integer.MAX_VALUE - threshold ? n : Math.min(
                    n, j + threshold);

            // ignore entry barcodeRev of leftmost ?????
            if (min > 1) {
                cost[min - 1] = Integer.MAX_VALUE-10;
            }

            // iterates through [min, max] in s
            for (int i = min; i <= max; i++) {
                if (barcodeRev[i - 1] == readJ || SequenceUtil.isNoCall(readJ) ) {
                    // diagonally left and up
                    cost[i] = previousCost[i - 1];
                } else {
                    final int snpCost = previousCost[i - 1];
                    final int delCost = cost[i - 1];
                    final int insCost = previousCost[i];

                    cost[i] = 1 + Math.min(Math.min(snpCost, delCost), insCost);
                }
            }

            // copy current distance counts to 'previous row' distance counts
            tempD = previousCost;
            previousCost = cost;
            cost = tempD;
        }

        // if previousCost[n] is greater than the threshold, there's no guarantee on it
        // being the correct
        // distance
        if (previousCost[n] <= threshold) {
            return previousCost[n];
        }
        return threshold + 1;
    }
}
