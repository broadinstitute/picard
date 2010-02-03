/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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
package net.sf.picard.util;

import java.util.ConcurrentModificationException;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * A Red-Black tree with intervals for keys.
 * Not thread-safe, and cannot be made so.
 *
 * 7/24/2008: This was copied from the tedUtils package.
 * IMPORTANT!!! It has been modified to use the Reseq way of
 * handling coordinates (end-inclusive).
 *
 * @author tsharpe
 */
public class IntervalTree<V> implements Iterable<IntervalTree.Node<V>>
{
    /**
     * Return the number of intervals in the tree.
     * @return The number of intervals.
     */
    public int size()
    {
        return mRoot == null ? 0 : mRoot.getSize();
    }

    /**
     * Remove all entries.
     */
    public void clear()
    {
        mRoot = null;
    }

    /**
     * Put a new interval into the tree (or update the value associated with an existing interval).
     * If the interval is novel, the special sentinel value is returned.
     * @param start The interval's start.
     * @param end The interval's end.
     * @param value The associated value.
     * @return The old value associated with that interval, or the sentinel.
     */
    @SuppressWarnings("null")
    public V put( final int start, final int end, final V value )
    {
        if ( start > end )
            throw new IllegalArgumentException("Start cannot exceed end.");

        V result = mSentinel;

        if ( mRoot == null )
        {
            mRoot = new Node<V>(start,end,value);
        }
        else
        {
            Node<V> parent = null;
            Node<V> node = mRoot;
            int cmpVal = 0;

            while ( node != null )
            {
                parent = node; // last non-null node
                cmpVal = node.compare(start,end);
                if ( cmpVal == 0 )
                {
                    break;
                }

                node = cmpVal < 0 ? node.getLeft() : node.getRight();
            }

            if ( cmpVal == 0 )
            {
                result = parent.setValue(value);
            }
            else
            {
                if ( cmpVal < 0 )
                {
                    mRoot = parent.insertLeft(start,end,value,mRoot);
                }
                else
                {
                    mRoot = parent.insertRight(start,end,value,mRoot);
                }
            }
        }

        return result;
    }

    /**
     * Remove an interval from the tree.  If the interval does not exist in the tree the
     * special sentinel value is returned.
     * @param start The interval's start.
     * @param end The interval's end.
     * @return The value associated with that interval, or the sentinel.
     */
    public V remove( final int start, final int end )
    {
        V result = mSentinel;
        Node<V> node = mRoot;

        while ( node != null )
        {
            final int cmpVal = node.compare(start,end);
            if ( cmpVal == 0 )
            {
                result = node.getValue();
                mRoot = node.remove(mRoot);
                break;
            }

            node = cmpVal < 0 ? node.getLeft() : node.getRight();
        }

        return result;
    }

    /**
     * Find an interval.
     * @param start The interval's start.
     * @param end The interval's end.
     * @return The Node that represents that interval, or null.
     */
    public Node<V> find( final int start, final int end )
    {
        Node<V> node = mRoot;

        while ( node != null )
        {
            final int cmpVal = node.compare(start,end);
            if ( cmpVal == 0 )
            {
                break;
            }

            node = cmpVal < 0 ? node.getLeft() : node.getRight();
        }

        return node;
    }

    /**
     * Find the nth interval in the tree.
     * @param idx The rank of the interval sought (from 0 to size()-1).
     * @return The Node that represents the nth interval.
     */
    public Node<V> findByIndex( final int idx )
    {
        return Node.findByRank(mRoot,idx+1);
    }

    /**
     * Find the rank of the specified interval.  If the specified interval is not in the
     * tree, then -1 is returned.
     * @param start The interval's start.
     * @param end The interval's end.
     * @return The rank of that interval, or -1.
     */
    public int getIndex( final int start, final int end )
    {
        return Node.getRank(mRoot,start,end) - 1;
    }

    /**
     * Find the least interval in the tree.
     * @return The earliest interval, or null if the tree is empty.
     */
    public Node<V> min()
    {
        Node<V> result = null;
        Node<V> node = mRoot;

        while ( node != null )
        {
            result = node;
            node = node.getLeft();
        }

        return result;
    }

    /**
     * Find the earliest interval in the tree greater than or equal to the specified interval.
     * @param start The interval's start.
     * @param end The interval's end.
     * @return The earliest >= interval, or null if there is none.
     */
    @SuppressWarnings("null")
    public Node<V> min( final int start, final int end )
    {
        Node<V> result = null;
        Node<V> node = mRoot;
        int cmpVal = 0;

        while ( node != null )
        {
            result = node;
            cmpVal = node.compare(start,end);
            if ( cmpVal == 0 )
            {
                break;
            }

            node = cmpVal < 0 ? node.getLeft() : node.getRight();
        }

        if ( cmpVal > 0 )
        {
            result = result.getNext();
        }

        return result;
    }

    /**
     * Find the earliest interval in the tree that overlaps the specified interval.
     * @param start The interval's start.
     * @param end The interval's end.
     * @return The earliest overlapping interval, or null if there is none.
     */
    public Node<V> minOverlapper( final int start, final int end )
    {
        Node<V> result = null;
        Node<V> node = mRoot;

        if ( node != null && node.getMaxEnd() >= start )
        {
            while ( true )
            {
                if ( node.getStart() <= end && start <= node.getEnd() )
                { // this node overlaps.  there might be a lesser overlapper down the left sub-tree.
                  // no need to consider the right sub-tree:  even if there's an overlapper, if won't be minimal
                    result = node;
                    node = node.getLeft();
                    if ( node == null || node.getMaxEnd() < start )
                        break; // no left sub-tree or all nodes end too early
                }
                else
                { // no overlap.  if there might be a left sub-tree overlapper, consider the left sub-tree.
                    final Node<V> left = node.getLeft();
                    if ( left != null && left.getMaxEnd() >= start )
                    {
                        node = left;
                    }
                    else
                    { // left sub-tree cannot contain an overlapper.  consider the right sub-tree.
                        if ( node.getStart() > end )
                            break; // everything in the right sub-tree is past the end of the query interval

                        node = node.getRight();
                        if ( node == null || node.getMaxEnd() < start )
                            break; // no right sub-tree or all nodes end too early
                    }
                }
            }
        }

        return result;
    }

    /**
     * Find the greatest interval in the tree.
     * @return The latest interval, or null if the tree is empty.
     */
    public Node<V> max()
    {
        Node<V> result = null;
        Node<V> node = mRoot;

        while ( node != null )
        {
            result = node;
            node = node.getRight();
        }

        return result;
    }

    /**
     * Find the latest interval in the tree less than or equal to the specified interval.
     * @param start The interval's start.
     * @param end The interval's end.
     * @return The latest >= interval, or null if there is none.
     */
    @SuppressWarnings("null")
    public Node<V> max( final int start, final int end )
    {
        Node<V> result = null;
        Node<V> node = mRoot;
        int cmpVal = 0;

        while ( node != null )
        {
            result = node;
            cmpVal = node.compare(start,end);
            if ( cmpVal == 0 )
            {
                break;
            }

            node = cmpVal < 0 ? node.getLeft() : node.getRight();
        }

        if ( cmpVal < 0 )
        {
            result = result.getPrev();
        }

        return result;
    }

    /**
     * Return an iterator over the entire tree.
     * @return An iterator.
     */
    public Iterator<Node<V>> iterator()
    {
        return new FwdIterator(min());
    }

    /**
     * Return an iterator over all intervals greater than or equal to the specified interval.
     * @param start The interval's start.
     * @param end The interval's end.
     * @return An iterator.
     */
    public Iterator<Node<V>> iterator( final int start, final int end )
    {
        return new FwdIterator(min(start,end));
    }

    /**
     * Return an iterator over all intervals overlapping the specified range.
     * @param start The range start.
     * @param end The range end.
     * @return An iterator.
     */
    public Iterator<Node<V>> overlappers( final int start, final int end )
    {
        return new OverlapIterator(start,end);
    }

    /**
     * Return an iterator over the entire tree that returns intervals in reverse order.
     * @return An iterator.
     */
    public Iterator<Node<V>> reverseIterator()
    {
        return new RevIterator(max());
    }

    /**
     * Return an iterator over all intervals less than or equal to the specified interval, in reverse order.
     * @param start The interval's start.
     * @param end The interval's end.
     * @return An iterator.
     */
    public Iterator<Node<V>> reverseIterator( final int start, final int end )
    {
        return new RevIterator(max(start,end));
    }

    /**
     * Get the special sentinel value that will be used to signal novelty when putting a new interval
     * into the tree, or to signal "not found" when removing an interval.  This is null by default.
     * @return The sentinel value.
     */
    public V getSentinel()
    {
        return mSentinel;
    }

    /**
     * Set the special sentinel value that will be used to signal novelty when putting a new interval
     * into the tree, or to signal "not found" when removing an interval.
     * @param sentinel The new sentinel value.
     * @return The old sentinel value.
     */
    public V setSentinel( final V sentinel )
    {
        final V result = mSentinel;
        mSentinel = sentinel;
        return result;
    }

    /**
     * This method is only for debugging.
     * It verifies whether the tree is internally consistent with respect to the mMaxEnd cached on each node.
     * @throws IllegalStateException If an inconsistency is detected.
     */
    public void checkMaxEnds() {
        if (mRoot != null) mRoot.checkMaxEnd();
    }

    /**
     * This method draws a nested picture of the tree on System.out.
     * Useful for debugging.
     */
    public void printTree() {
        if (mRoot != null) mRoot.printNode();
    }

    void removeNode( final Node<V> node )
    {
        mRoot = node.remove(mRoot);
    }

    private Node<V> mRoot;
    private V mSentinel;

    public static class Node<V1>
    {
        // bit-wise definitions from which the other constants are composed
        public static final int HAS_LESSER_PART = 1;
        public static final int HAS_OVERLAPPING_PART = 2;
        public static final int HAS_GREATER_PART = 4;
        public static final int IS_ADJACENT_AND_EMPTY = 0;
        public static final int IS_STRICTLY_LESS = HAS_LESSER_PART; // 1
        public static final int IS_SUBSET = HAS_OVERLAPPING_PART; // 2
        public static final int IS_LEFT_OVERHANGING_OVERLAPPER = HAS_LESSER_PART | HAS_OVERLAPPING_PART; // 3
        public static final int IS_STRICTLY_GREATER = HAS_GREATER_PART; // 4
        // there is no value that equals 5, since that would imply overhanging on left and right without overlapping
        public static final int IS_RIGHT_OVERHANGING_OVERLAPPER = HAS_GREATER_PART | HAS_OVERLAPPING_PART; // 6
        public static final int IS_SUPERSET = HAS_LESSER_PART | HAS_OVERLAPPING_PART | HAS_GREATER_PART; // 7

        Node( final int start, final int end, final V1 value )
        {
            mStart = start;
            mEnd = end;
            mValue = value;
            mSize = 1;
            mMaxEnd = mEnd;
            mIsBlack = true;
        }

        Node( final Node<V1> parent, final int start, final int end, final V1 value )
        {
            mParent = parent;
            mStart = start;
            mEnd = end;
            mValue = value;
            mMaxEnd = mEnd;
            mSize = 1;
        }

        public int getStart()
        {
            return mStart;
        }

        public int getEnd()
        {
            return mEnd;
        }

        public int getLength()
        {
            return mEnd - mStart;
        }

        public int getRelationship( final Node<V1> interval )
        {
            int result = 0;
            if ( mStart < interval.getStart() )
                result = HAS_LESSER_PART;
            if ( mEnd > interval.getEnd() )
                result |= HAS_GREATER_PART;
            if ( mStart < interval.getEnd() && interval.getStart() < mEnd )
                result |= HAS_OVERLAPPING_PART;
            return result;
        }

        public boolean isAdjacent( final Node<V1> interval )
        {
            return mStart == interval.getEnd() || mEnd == interval.getStart();
        }

        public V1 getValue()
        {
            return mValue;
        }

        public V1 setValue( final V1 value )
        {
            final V1 result = mValue;
            mValue = value;
            return result;
        }

        int getSize()
        {
            return mSize;
        }

        int getMaxEnd()
        {
            return mMaxEnd;
        }

        Node<V1> getLeft()
        {
            return mLeft;
        }

        Node<V1> insertLeft( final int start, final int end, final V1 value, final Node<V1> root )
        {
            mLeft = new Node<V1>(this,start,end,value);
            return insertFixup(mLeft,root);
        }

        Node<V1> getRight()
        {
            return mRight;
        }

        Node<V1> insertRight( final int start, final int end, final V1 value, final Node<V1> root )
        {
            mRight = new Node<V1>(this,start,end,value);
            return insertFixup(mRight,root);
        }

        Node<V1> getNext()
        {
            Node<V1> result;

            if ( mRight != null )
            {
                result = mRight;
                while ( result.mLeft != null )
                {
                    result = result.mLeft;
                }
            }
            else
            {
                Node<V1> node = this;
                result = mParent;
                while ( result != null && node == result.mRight )
                {
                    node = result;
                    result = result.mParent;
                }
            }

            return result;
        }

        Node<V1> getPrev()
        {
            Node<V1> result;

            if ( mLeft != null )
            {
                result = mLeft;
                while ( result.mRight != null )
                {
                    result = result.mRight;
                }
            }
            else
            {
                Node<V1> node = this;
                result = mParent;
                while ( result != null && node == result.mLeft )
                {
                    node = result;
                    result = result.mParent;
                }
            }

            return result;
        }

        boolean wasRemoved()
        {
            return mSize == 0;
        }

        Node<V1> remove( Node<V1> root )
        {
            if ( mSize == 0 )
            {
                throw new IllegalStateException("Entry was already removed.");
            }

            if ( mLeft == null )
            {
                if ( mRight == null )
                { // no children
                    if ( mParent == null )
                    {
                        root = null;
                    }
                    else if ( mParent.mLeft == this )
                    {
                        mParent.mLeft = null;
                        fixup(mParent);

                        if ( mIsBlack )
                            root = removeFixup(mParent,null,root);
                    }
                    else
                    {
                        mParent.mRight = null;
                        fixup(mParent);

                        if ( mIsBlack )
                            root = removeFixup(mParent,null,root);
                    }
                }
                else
                { // single child on right
                    root = spliceOut(mRight,root);
                }
            }
            else if ( mRight == null )
            { // single child on left
                root = spliceOut(mLeft,root);
            }
            else
            { // two children
                final Node<V1> next = getNext();
                root = next.remove(root);

                // put next into tree in same position as this, effectively removing this
                if ( (next.mParent = mParent) == null )
                    root = next;
                else if ( mParent.mLeft == this )
                    mParent.mLeft = next;
                else
                    mParent.mRight = next;

                if ( (next.mLeft = mLeft) != null )
                {
                    mLeft.mParent = next;
                }

                if ( (next.mRight = mRight) != null )
                {
                    mRight.mParent = next;
                }

                next.mIsBlack = mIsBlack;
                next.mSize = mSize;

                // PIC-123 fix
                fixup(next);
            }

            mSize = 0;
            return root;
        }

        // backwards comparison!  compares start+end to this.
        int compare( final int start, final int end )
        {
            int result = 0;

            if ( start > mStart )
                result = 1;
            else if ( start < mStart )
                result = -1;
            else if ( end > mEnd )
                result = 1;
            else if ( end < mEnd )
                result = -1;

            return result;
        }

        @SuppressWarnings("null")
        static <V1> Node<V1> getNextOverlapper( Node<V1> node, final int start, final int end )
        {
            do
            {
                Node<V1> nextNode = node.mRight;
                if ( nextNode != null && nextNode.mMaxEnd >= start )
                {
                    node = nextNode;
                    while ( (nextNode = node.mLeft) != null && nextNode.mMaxEnd >= start )
                        node = nextNode;
                }
                else
                {
                    nextNode = node;
                    while ( (node = nextNode.mParent) != null && node.mRight == nextNode )
                        nextNode = node;
                }

                if ( node != null && node.mStart > end )
                    node = null;
            }
            while ( node != null && !(node.mStart <= end && start <= node.mEnd) );

            return node;
        }

        static <V1> Node<V1> findByRank( Node<V1> node, int rank )
        {
            while ( node != null )
            {
                final int nodeRank = node.getRank();
                if ( rank == nodeRank )
                    break;

                if ( rank < nodeRank )
                {
                    node = node.mLeft;
                }
                else
                {
                    node = node.mRight;
                    rank -= nodeRank;
                }
            }

            return node;
        }

        static <V1> int getRank( Node<V1> node, final int start, final int end )
        {
            int rank = 0;

            while ( node != null )
            {
                final int cmpVal = node.compare(start,end);
                if ( cmpVal < 0 )
                {
                    node = node.mLeft;
                }
                else
                {
                    rank += node.getRank();
                    if ( cmpVal == 0 )
                        return rank; // EARLY RETURN!!!

                    node = node.mRight;
                }
            }

            return 0;
        }

        private int getRank()
        {
            int result = 1;
            if ( mLeft != null )
                result = mLeft.mSize + 1;
            return result;
        }

        private Node<V1> spliceOut( final Node<V1> child, Node<V1> root )
        {
            if ( (child.mParent = mParent) == null )
            {
                root = child;
                child.mIsBlack = true;
            }
            else
            {
                if ( mParent.mLeft == this )
                    mParent.mLeft = child;
                else
                    mParent.mRight = child;
                fixup(mParent);

                if ( mIsBlack )
                    root = removeFixup(mParent,child,root);
            }

            return root;
        }

        private Node<V1> rotateLeft( Node<V1> root )
        {
            final Node<V1> child = mRight;

            final int childSize = child.mSize;
            child.mSize = mSize;
            mSize -= childSize;

            if ( (mRight = child.mLeft) != null )
            {
                mRight.mParent = this;
                mSize += mRight.mSize;
            }

            if ( (child.mParent = mParent) == null )
                root = child;
            else if ( this == mParent.mLeft )
                mParent.mLeft = child;
            else
                mParent.mRight = child;

            child.mLeft = this;
            mParent = child;

            setMaxEnd();
            child.setMaxEnd();

            return root;
        }

        private Node<V1> rotateRight( Node<V1> root )
        {
            final Node<V1> child = mLeft;

            final int childSize = child.mSize;
            child.mSize = mSize;
            mSize -= childSize;

            if ( (mLeft = child.mRight) != null )
            {
                mLeft.mParent = this;
                mSize += mLeft.mSize;
            }

            if ( (child.mParent = mParent) == null )
                root = child;
            else if ( this == mParent.mLeft )
                mParent.mLeft = child;
            else
                mParent.mRight = child;

            child.mRight = this;
            mParent = child;

            setMaxEnd();
            child.setMaxEnd();

            return root;
        }

        private void setMaxEnd()
        {
            mMaxEnd = mEnd;
            if ( mLeft != null )
                mMaxEnd = Math.max(mMaxEnd,mLeft.mMaxEnd);
            if ( mRight != null )
                mMaxEnd = Math.max(mMaxEnd,mRight.mMaxEnd);
        }

        private static <V1> void fixup( Node<V1> node )
        {
            do
            {
                node.mSize = 1;
                node.mMaxEnd = node.mEnd;
                if ( node.mLeft != null )
                {
                    node.mSize += node.mLeft.mSize;
                    node.mMaxEnd = Math.max(node.mMaxEnd,node.mLeft.mMaxEnd);
                }
                if ( node.mRight != null )
                {
                    node.mSize += node.mRight.mSize;
                    node.mMaxEnd = Math.max(node.mMaxEnd,node.mRight.mMaxEnd);
                }
            }
            while ( (node = node.mParent) != null );
        }

        private static <V1> Node<V1> insertFixup( Node<V1> daughter, Node<V1> root )
        {
            Node<V1> mom = daughter.mParent;
            fixup(mom);

            while( mom != null && !mom.mIsBlack )
            {
                final Node<V1> gramma = mom.mParent;
                Node<V1> auntie = gramma.mLeft;
                if ( auntie == mom )
                {
                    auntie = gramma.mRight;
                    if ( auntie != null && !auntie.mIsBlack )
                    {
                        mom.mIsBlack = true;
                        auntie.mIsBlack = true;
                        gramma.mIsBlack = false;
                        daughter = gramma;
                    }
                    else
                    {
                        if ( daughter == mom.mRight )
                        {
                            root = mom.rotateLeft(root);
                            mom = daughter;
                        }
                        mom.mIsBlack = true;
                        gramma.mIsBlack = false;
                        root = gramma.rotateRight(root);
                        break;
                    }
                }
                else
                {
                    if ( auntie != null && !auntie.mIsBlack )
                    {
                        mom.mIsBlack = true;
                        auntie.mIsBlack = true;
                        gramma.mIsBlack = false;
                        daughter = gramma;
                    }
                    else
                    {
                        if ( daughter == mom.mLeft )
                        {
                            root = mom.rotateRight(root);
                            mom = daughter;
                        }
                        mom.mIsBlack = true;
                        gramma.mIsBlack = false;
                        root = gramma.rotateLeft(root);
                        break;
                    }
                }
                mom = daughter.mParent;
            }
            root.mIsBlack = true;
            return root;
        }

        private static <V1> Node<V1> removeFixup( Node<V1> parent, Node<V1> node, Node<V1> root )
        {
            do
            {
                if ( node == parent.mLeft )
                {
                    Node<V1> sister = parent.mRight;
                    if ( !sister.mIsBlack )
                    {
                        sister.mIsBlack = true;
                        parent.mIsBlack = false;
                        root = parent.rotateLeft(root);
                        sister = parent.mRight;
                    }
                    if ( (sister.mLeft == null || sister.mLeft.mIsBlack) && (sister.mRight == null || sister.mRight.mIsBlack) )
                    {
                        sister.mIsBlack = false;
                        node = parent;
                    }
                    else
                    {
                        if ( sister.mRight == null || sister.mRight.mIsBlack )
                        {
                            sister.mLeft.mIsBlack = true;
                            sister.mIsBlack = false;
                            root = sister.rotateRight(root);
                            sister = parent.mRight;
                        }
                        sister.mIsBlack = parent.mIsBlack;
                        parent.mIsBlack = true;
                        sister.mRight.mIsBlack = true;
                        root = parent.rotateLeft(root);
                        node = root;
                    }
                }
                else
                {
                    Node<V1> sister = parent.mLeft;
                    if ( !sister.mIsBlack )
                    {
                        sister.mIsBlack = true;
                        parent.mIsBlack = false;
                        root = parent.rotateRight(root);
                        sister = parent.mLeft;
                    }
                    if ( (sister.mLeft == null || sister.mLeft.mIsBlack) && (sister.mRight == null || sister.mRight.mIsBlack) )
                    {
                        sister.mIsBlack = false;
                        node = parent;
                    }
                    else
                    {
                        if ( sister.mLeft == null || sister.mLeft.mIsBlack )
                        {
                            sister.mRight.mIsBlack = true;
                            sister.mIsBlack = false;
                            root = sister.rotateLeft(root);
                            sister = parent.mLeft;
                        }
                        sister.mIsBlack = parent.mIsBlack;
                        parent.mIsBlack = true;
                        sister.mLeft.mIsBlack = true;
                        root = parent.rotateRight(root);
                        node = root;
                    }
                }
                parent = node.mParent;
            }
            while ( parent != null && node.mIsBlack );

            node.mIsBlack = true;
            return root;
        }

        public void checkMaxEnd() {
            if (mMaxEnd != calcMaxEnd()) {
                throw new IllegalStateException("Max end mismatch " + mMaxEnd + " vs " + calcMaxEnd() + ": " + this);
            }
            if (mLeft != null) mLeft.checkMaxEnd();
            if (mRight != null) mRight.checkMaxEnd();
        }

        private int calcMaxEnd() {
            int end = mEnd;
            if (mLeft != null) end = Math.max(end, mLeft.mMaxEnd);
            if (mRight != null) end = Math.max(end, mRight.mMaxEnd);
            return end;
        }

        public void printNode() {
            this.printNodeInternal("", "root: ");
        }

        private void printNodeInternal(final String padding, final String tag) {
            System.out.println(padding + tag + " " + this);
            if (mLeft != null) mLeft.printNodeInternal(padding + "  ", "left: ");
            if (mRight != null) mRight.printNodeInternal(padding + "  ", "right:");
        }

        public String toString() {
            return "Node(" + mStart + "," + mEnd + "," + mValue + "," + mSize + "," + mMaxEnd + "," + mIsBlack + ")";
        }

        private Node<V1> mParent;
        private Node<V1> mLeft;
        private Node<V1> mRight;
        private final int mStart;
        private final int mEnd;
        private V1 mValue;
        private int mSize;
        private int mMaxEnd;
        private boolean mIsBlack;
    }

    public class FwdIterator
        implements Iterator<Node<V>>
    {
        public FwdIterator( final Node<V> node )
        {
            mNext = node;
        }

        public boolean hasNext()
        {
            return mNext != null;
        }

        public Node<V> next()
        {
            if ( mNext == null )
            {
                throw new NoSuchElementException("No next element.");
            }

            if ( mNext.wasRemoved() )
            {
                mNext = min(mNext.getStart(),mNext.getEnd());
                if ( mNext == null )
                    throw new ConcurrentModificationException("Current element was removed, and there are no more elements.");
            }
            mLast = mNext;
            mNext = mNext.getNext();
            return mLast;
        }

        public void remove()
        {
            if ( mLast == null )
            {
                throw new IllegalStateException("No entry to remove.");
            }

            removeNode(mLast);
            mLast = null;
        }

        private Node<V> mNext;
        private Node<V> mLast;
    }

    public class RevIterator
        implements Iterator<Node<V>>
    {
        public RevIterator( final Node<V> node )
        {
            mNext = node;
        }

        public boolean hasNext()
        {
            return mNext != null;
        }

        public Node<V> next()
        {
            if ( mNext == null )
                throw new NoSuchElementException("No next element.");
            if ( mNext.wasRemoved() )
            {
                mNext = max(mNext.getStart(),mNext.getEnd());
                if ( mNext == null )
                    throw new ConcurrentModificationException("Current element was removed, and there are no more elements.");
            }
            mLast = mNext;
            mNext = mNext.getPrev();
            return mLast;
        }

        public void remove()
        {
            if ( mLast == null )
            {
                throw new IllegalStateException("No entry to remove.");
            }

            removeNode(mLast);
            mLast = null;
        }

        private Node<V> mNext;
        private Node<V> mLast;
    }

    public class OverlapIterator
        implements Iterator<Node<V>>
    {
        public OverlapIterator( final int start, final int end )
        {
            mNext = minOverlapper(start,end);
            mStart = start;
            mEnd = end;
        }

        public boolean hasNext()
        {
            return mNext != null;
        }

        public Node<V> next()
        {
            if ( mNext == null )
            {
                throw new NoSuchElementException("No next element.");
            }

            if ( mNext.wasRemoved() )
            {
                throw new ConcurrentModificationException("Current element was removed.");
            }

            mLast = mNext;
            mNext = Node.getNextOverlapper(mNext,mStart,mEnd);
            return mLast;
        }

        public void remove()
        {
            if ( mLast == null )
            {
                throw new IllegalStateException("No entry to remove.");
            }

            removeNode(mLast);
            mLast = null;
        }

        private Node<V> mNext;
        private Node<V> mLast;
        private final int mStart;
        private final int mEnd;
    }

    public static class ValuesIterator<V1>
        implements Iterator<V1>
    {
        public ValuesIterator( final Iterator<Node<V1>> itr )
        {
            mItr = itr;
        }

        public boolean hasNext()
        {
            return mItr.hasNext();
        }

        public V1 next()
        {
            return mItr.next().getValue();
        }

        public void remove()
        {
            mItr.remove();
        }

        private final Iterator<Node<V1>> mItr;
    }
}
