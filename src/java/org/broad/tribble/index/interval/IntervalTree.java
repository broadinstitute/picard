/*
 * Copyright (c) 2007-2010 by The Broad Institute, Inc. and the Massachusetts Institute of Technology.
 * All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which
 * is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR WARRANTIES OF
 * ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT
 * OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR
 * RESPECTIVE TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES OF
 * ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES, ECONOMIC
 * DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER THE BROAD OR MIT SHALL
 * BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
 */

package org.broad.tribble.index.interval;


import java.io.DataOutputStream;
import java.io.IOException;
import java.util.*;

/** An implementation of an interval tree, following the explanation.
 * from CLR. For efficiently finding all intervals which overlap a given
 * interval or point.
 *
 * References:
 * http://en.wikipedia.org/wiki/Interval_tree
 *
 * Cormen, Thomas H.; Leiserson, Charles E., Rivest, Ronald L. (1990). Introduction to Algorithms (1st ed.). MIT Press and McGraw-Hill. ISBN 0-262-03141-8
 *
 *
 */
public class IntervalTree {

    Node root;
    Node NIL = Node.NIL;

    /**
     * See {@link #getSize()}
     */
    int size;

    public IntervalTree() {
        this.root = NIL;
        this.size = 0;
    }


    public void insert(Interval interval) {
        Node node = new Node(interval);
        insert(node);
        size++;
    }

    /**
     * The estimated size of the tree. We keep a running count
     * on each insert, this getter returns that count.
     * @see #size()
     * @return
     */
    public int getSize() {
        return size;
    }

    /**
     *
     * @param interval
     * @return all matches as a list of Intervals
     */
    public List<Interval> findOverlapping(Interval interval) {

        if (root().isNull()) {
            return Collections.emptyList();
        }

        List<Interval> results = new ArrayList();
        searchAll(interval, root(), results);
        return results;
    }

    public String toString() {
        return root().toString();
    }

    private List<Interval> searchAll(Interval interval, Node node, List<Interval> results) {
        if (node.interval.overlaps(interval)) {
            results.add(node.interval);
        }
        if (!node.left.isNull() && node.left.max >= interval.start) {
            searchAll(interval, node.left, results);
        }
        if (!node.right.isNull() && node.right.min <= interval.end) {
            searchAll(interval, node.right, results);
        }
        return results;
    }

    /**
     * Return all intervals in tree.
     * TODO: an iterator would be more effecient.
     * @return
     */
    public List<Interval> getIntervals() {
        if (root().isNull()) {
            return Collections.emptyList();
        }
        List<Interval> results = new ArrayList(size);
        getAll(root(), results);
        return results;
    }

    /**
     * Get all nodes which are descendants of {@code node}, inclusive.
     * {@code results} is modified in place
     * @param node
     * @param results
     * @return the total list of descendants, including original {@code results}
     */
    private List<Interval> getAll(Node node, List<Interval> results) {

        results.add(node.interval);
        if (!node.left.isNull()) {
            getAll(node.left, results);
        }
        if (!node.right.isNull()) {
            getAll(node.right, results);
        }
        return results;
    }


    /**
     * Used for testing only.
     *
     * @param node
     * @return
     */
    private int getRealMax(Node node) {
        if (node.isNull())
            return Integer.MIN_VALUE;
        int leftMax = getRealMax(node.left);
        int rightMax = getRealMax(node.right);
        int nodeHigh = (node.interval).end;

        int max1 = (leftMax > rightMax ? leftMax : rightMax);
        return (max1 > nodeHigh ? max1 : nodeHigh);
    }

    /**
     * Used for testing only
     *
     * @param node
     * @return
     */
    private int getRealMin(Node node) {
        if (node.isNull())
            return Integer.MAX_VALUE;

        int leftMin = getRealMin(node.left);
        int rightMin = getRealMin(node.right);
        int nodeLow = (node.interval).start;

        int min1 = (leftMin < rightMin ? leftMin : rightMin);
        return (min1 < nodeLow ? min1 : nodeLow);
    }


    private void insert(Node x) {
        assert (x != null);
        assert (!x.isNull());

        treeInsert(x);
        x.color = Node.RED;
        while (x != this.root && x.parent.color == Node.RED) {
            if (x.parent == x.parent.parent.left) {
                Node y = x.parent.parent.right;
                if (y.color == Node.RED) {
                    x.parent.color = Node.BLACK;
                    y.color = Node.BLACK;
                    x.parent.parent.color = Node.RED;
                    x = x.parent.parent;
                } else {
                    if (x == x.parent.right) {
                        x = x.parent;
                        this.leftRotate(x);
                    }
                    x.parent.color = Node.BLACK;
                    x.parent.parent.color = Node.RED;
                    this.rightRotate(x.parent.parent);
                }
            } else {
                Node y = x.parent.parent.left;
                if (y.color == Node.RED) {
                    x.parent.color = Node.BLACK;
                    y.color = Node.BLACK;
                    x.parent.parent.color = Node.RED;
                    x = x.parent.parent;
                } else {
                    if (x == x.parent.left) {
                        x = x.parent;
                        this.rightRotate(x);
                    }
                    x.parent.color = Node.BLACK;
                    x.parent.parent.color = Node.RED;
                    this.leftRotate(x.parent.parent);
                }
            }
        }
        this.root.color = Node.BLACK;
    }


    private Node root() {
        return this.root;
    }


    private void leftRotate(Node x) {
        Node y = x.right;
        x.right = y.left;
        if (y.left != NIL) {
            y.left.parent = x;
        }
        y.parent = x.parent;
        if (x.parent == NIL) {
            this.root = y;
        } else {
            if (x.parent.left == x) {
                x.parent.left = y;
            } else {
                x.parent.right = y;
            }
        }
        y.left = x;
        x.parent = y;

        applyUpdate(x);
        // no need to apply update on y, since it'll y is an ancestor
        // of x, and will be touched by applyUpdate().
    }


    private void rightRotate(Node x) {
        Node y = x.left;
        x.left = y.right;
        if (y.right != NIL) {
            y.right.parent = x;
        }
        y.parent = x.parent;
        if (x.parent == NIL) {
            this.root = y;
        } else {
            if (x.parent.right == x) {
                x.parent.right = y;
            } else {
                x.parent.left = y;
            }
        }
        y.right = x;
        x.parent = y;


        applyUpdate(x);
        // no need to apply update on y, since it'll y is an ancestor
        // of x, and will be touched by applyUpdate().
    }


    /**
     * Note:  Does not maintain RB constraints,  this is done post insert
     *
     * @param x
     */
    private void treeInsert(Node x) {
        Node node = this.root;
        Node y = NIL;
        while (node != NIL) {
            y = node;
            if (x.interval.start <= node.interval.start) {
                node = node.left;
            } else {
                node = node.right;
            }
        }
        x.parent = y;

        if (y == NIL) {
            this.root = x;
            x.left = x.right = NIL;
        } else {
            if (x.interval.start <= y.interval.start) {
                y.left = x;
            } else {
                y.right = x;
            }
        }

        this.applyUpdate(x);
    }


    // Applies the statistic update on the node and its ancestors.

    private void applyUpdate(Node node) {
        while (!node.isNull()) {
            this.update(node);
            node = node.parent;
        }
    }

    private void update(Node node) {
        node.max = Math.max(Math.max(node.left.max, node.right.max), node.interval.end);
        node.min = Math.min(Math.min(node.left.min, node.right.min), node.interval.start);
    }

    /**
     * @see #getSize()
     * @return Returns the number of nodes in the tree.
     *         Recalculated each call
     */
    public int size() {
        return _size(this.root);
    }


    private int _size(Node node) {
        if (node.isNull())
            return 0;
        return 1 + _size(node.left) + _size(node.right);
    }


    private boolean allRedNodesFollowConstraints(Node node) {
        if (node.isNull())
            return true;

        if (node.color == Node.BLACK) {
            return (allRedNodesFollowConstraints(node.left) &&
                    allRedNodesFollowConstraints(node.right));
        }

        // At this point, we know we're on a RED node.
        return (node.left.color == Node.BLACK &&
                node.right.color == Node.BLACK &&
                allRedNodesFollowConstraints(node.left) &&
                allRedNodesFollowConstraints(node.right));
    }


    // Check that both ends are equally balanced in terms of black height.

    private boolean isBalancedBlackHeight(Node node) {
        if (node.isNull())
            return true;
        return (blackHeight(node.left) == blackHeight(node.right) &&
                isBalancedBlackHeight(node.left) &&
                isBalancedBlackHeight(node.right));
    }


    // The black height of a node should be left/right equal.

    private int blackHeight(Node node) {
        if (node.isNull())
            return 0;
        int leftBlackHeight = blackHeight(node.left);
        if (node.color == Node.BLACK) {
            return leftBlackHeight + 1;
        } else {
            return leftBlackHeight;
        }
    }


    /**
     * Test code: make sure that the tree has all the properties
     * defined by Red Black trees and interval trees
     * <p/>
     * o.  Root is black.
     * <p/>
     * o.  NIL is black.
     * <p/>
     * o.  Red nodes have black children.
     * <p/>
     * o.  Every path from root to leaves contains the same number of
     * black nodes.
     * <p/>
     * o.  getMax(node) is the maximum of any interval rooted at that node..
     * <p/>
     * This code is expensive, and only meant to be used for
     * assertions and testing.
     */
    public boolean isValid() {
        if (this.root.color != Node.BLACK) {
            //logger.warn("root color is wrong");
            return false;
        }
        if (NIL.color != Node.BLACK) {
            //logger.warn("NIL color is wrong");
            return false;
        }
        if (allRedNodesFollowConstraints(this.root) == false) {
            //logger.warn("red node doesn't follow constraints");
            return false;
        }
        if (isBalancedBlackHeight(this.root) == false) {
            //logger.warn("black height unbalanced");
            return false;
        }

        return hasCorrectMaxFields(this.root) &&
                hasCorrectMinFields(this.root);
    }


    private boolean hasCorrectMaxFields(Node node) {
        if (node.isNull())
            return true;
        return (getRealMax(node) == (node.max) &&
                hasCorrectMaxFields(node.left) &&
                hasCorrectMaxFields(node.right));
    }


    private boolean hasCorrectMinFields(Node node) {
        if (node.isNull())
            return true;
        return (getRealMin(node) == (node.min) &&
                hasCorrectMinFields(node.left) &&
                hasCorrectMinFields(node.right));
    }


    static class Node {

        public static boolean BLACK = false;
        public static boolean RED = true;

        Interval interval;
        int min;
        int max;
        Node left;
        Node right;

        // Color and parent are used for inserts.  If tree is immutable these are not required (no requirement
        // to store these persistently).
        boolean color;
        Node parent;


        private Node() {
            this.max = Integer.MIN_VALUE;
            this.min = Integer.MAX_VALUE;
        }

        public void store(DataOutputStream dos) throws IOException {
            dos.writeInt(interval.start);
            dos.writeInt(interval.end);
            dos.writeInt(min);
            dos.writeInt(max);

        }

        public Node(Interval interval) {
            this();
            this.parent = NIL;
            this.left = NIL;
            this.right = NIL;
            this.interval = interval;
            this.color = RED;
        }


        static Node NIL;

        static {
            NIL = new Node();
            NIL.color = BLACK;
            NIL.parent = NIL;
            NIL.left = NIL;
            NIL.right = NIL;
        }


        public boolean isNull() {
            return this == NIL;
        }


        public String toString() {

            // Make some shorthand for the nodes
            Map<Interval, Integer> keys = new LinkedHashMap();

            if (this == NIL) {
                return "nil";
            }

            StringBuffer buf = new StringBuffer();
            _toString(buf, keys);

            buf.append("\n");
            for (Map.Entry<Interval, Integer> entry : keys.entrySet()) {
                buf.append(entry.getValue() + " = " + entry.getKey());
                buf.append("\n");
            }

            return buf.toString();
        }

        public void _toString(StringBuffer buf, Map<Interval, Integer> keys) {
            if (this == NIL) {
                buf.append("nil");
                buf.append("\n");
                return;
            }

            Integer selfKey = keys.get(this.interval);
            if (selfKey == null) {
                selfKey = keys.size();
                keys.put(this.interval, selfKey);
            }
            Integer leftKey = keys.get(this.left.interval);
            if (leftKey == null) {
                leftKey = keys.size();
                keys.put(this.left.interval, leftKey);
            }
            Integer rightKey = keys.get(this.right.interval);
            if (rightKey == null) {
                rightKey = keys.size();
                keys.put(this.right.interval, rightKey);
            }


            buf.append(selfKey + " -> " + leftKey + " , " + rightKey);
            buf.append("\n");
            this.left._toString(buf, keys);
            this.right._toString(buf, keys);
        }
    }
}

