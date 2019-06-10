package picard.util;
import picard.util.BKNode;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

public class BKTree implements Serializable {
    private BKNode root;

    public BKTree(final String root_string) {
        root = new BKNode(root_string);
    }

    public BKTree() {}

    public void insert(final String newString) {
        if (root != null) {
            root.insert(newString);
        }
        else {
            root = new BKNode(newString);
        }
    }

    public HashMap<String, Integer> query(final String queryString, final int maxDist) {
        return root.query(queryString, maxDist);
    }

    public ArrayList<String> queryBest(final String queryString, final int maxDist) {
        ArrayList<String> ret = new ArrayList<>();

        HashMap<String, Integer> queryRet = query(queryString, maxDist);
        int minDist = maxDist;
        for (HashMap.Entry entry : queryRet.entrySet()) {
            if ((Integer)entry.getValue() <= minDist) {
                minDist = (Integer)entry.getValue();
                ret.add((String)entry.getKey());
            }
        }
        return ret;
    }
}
