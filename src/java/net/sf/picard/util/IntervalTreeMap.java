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

import java.util.*;

/**
 * Utility class that implements an interval map.
 * This class functions as a java map but also supports efficient interval overlap queries.
 *
 * @author Bob Handsaker
 */
public class IntervalTreeMap<T>
    extends AbstractMap<Interval, T>
{
    private final Map<String, IntervalTree<T>> mSequenceMap = new HashMap<String, IntervalTree<T>>();
    private final EntrySet mEntrySet = new EntrySet();

    public IntervalTree<T> debugGetTree(final String sequence) {
        return mSequenceMap.get(sequence);
    }

    public IntervalTreeMap() {
    }

    public IntervalTreeMap(final Map<? extends Interval, ? extends T> map) {
        for (final Map.Entry<? extends Interval, ? extends T> entry : map.entrySet()) {
            put(entry.getKey(), entry.getValue());
        }
    }

    public void clear() {
        mSequenceMap.clear();
    }

    public boolean containsKey(final Object object) {
        if (!(object instanceof Interval)) {
            return false;
        }
        return containsKey((Interval) object);
    }

    public boolean containsKey(final Interval key) {
        final IntervalTree<T> tree = mSequenceMap.get(key.getSequence());
        if (tree == null) {
            return false;
        }
        return (tree.find(key.getStart(), key.getEnd()) != null);
    }

    public Set<Entry<Interval, T>> entrySet() {
        return mEntrySet;
    }

    public boolean equals(final Object o) {
        if (!(o instanceof IntervalTreeMap)) {
            return false;
        }
        return mSequenceMap.equals(((IntervalTreeMap)o).mSequenceMap);
    }

    public int hashCode() {
        return mSequenceMap.hashCode();
    }

    public T get(final Object object) {
        if (!(object instanceof Interval)) {
            return null;
        }
        return get((Interval) object);
    }

    public T get(final Interval key) {
        final IntervalTree<T> tree = mSequenceMap.get(key.getSequence());
        if (tree == null) {
            return null;
        }
        final IntervalTree.Node<T> node = tree.find(key.getStart(), key.getEnd());
        if (node == null) {
            return null;
        }
        return node.getValue();
    }

    public boolean isEmpty() {
        for (final IntervalTree<T> tree : mSequenceMap.values()) {
            if (tree.size() > 0) {
                return false;
            }
        }
        return true;
    }

    public T put(final Interval key, final T value) {
        IntervalTree<T> tree = mSequenceMap.get(key.getSequence());
        if (tree == null) {
            tree = new IntervalTree<T>();
            mSequenceMap.put(key.getSequence(), tree);
        }
        return tree.put(key.getStart(), key.getEnd(), value);
    }

    public T remove(final Object object) {
        if (!(object instanceof Interval)) {
            return null;
        }
        return remove((Interval)object);
    }

    public T remove(final Interval key) {
        final IntervalTree<T> tree = mSequenceMap.get(key.getSequence());
        if (tree == null) {
            return null;
        }
        return tree.remove(key.getStart(), key.getEnd());
    }

    public int size() {
        // Note: We should think about caching the size to avoid having to recompute it.
        int size = 0;
        for (final IntervalTree<T> tree : mSequenceMap.values()) {
            size += tree.size();
        }
        return size;
    }

    public Collection<T> getOverlapping(final Interval key) {
        final List<T> result = new ArrayList<T>();
        final IntervalTree<T> tree = mSequenceMap.get(key.getSequence());
        if (tree != null) {
            final Iterator<IntervalTree.Node<T>> iterator = tree.overlappers(key.getStart(), key.getEnd());
            while (iterator.hasNext()) {
                result.add(iterator.next().getValue());
            }
        }
        return result;
    }

    public Collection<T> getContained(final Interval key) {
        final List<T> result = new ArrayList<T>();
        final IntervalTree<T> tree = mSequenceMap.get(key.getSequence());
        if (tree != null) {
            final Iterator<IntervalTree.Node<T>> iterator = tree.overlappers(key.getStart(), key.getEnd());
            while (iterator.hasNext()) {
                final IntervalTree.Node<T> node = iterator.next();
                if (node.getStart() >= key.getStart() && node.getEnd() <= key.getEnd()) {
                    result.add(node.getValue());
                }
            }
        }
        return result;
    }

    private class EntrySet
        extends AbstractSet<Map.Entry<Interval,T>> {

        public void clear() {
           IntervalTreeMap.this.clear();
        }

        public boolean contains(final Map.Entry<Interval,T> entry) {
            if (entry == null) {
                return false;
            }
            return entry.getValue().equals(IntervalTreeMap.this.get(entry.getKey()));
        }

        public boolean isEmpty() {
            return IntervalTreeMap.this.isEmpty();
        }

        public Iterator<Map.Entry<Interval,T>> iterator() {
            return new EntryIterator();
        }

        @SuppressWarnings("unchecked")
        public boolean remove(final Object object) {
            // Note: Could not figure out how to eliminate the unchecked cast.
            if (!(object instanceof Map.Entry)) {
                return false;
            }
            return remove((Map.Entry<Interval,T>)object);
        }

        public boolean remove(final Map.Entry<Interval,T> entry) {
            if (this.contains(entry)) {
                IntervalTreeMap.this.remove(entry.getKey());
                return true;
            } else {
                return false;
            }
        }

        public int size() {
            return IntervalTreeMap.this.size();
        }
    }

    private class EntryIterator
        implements Iterator<Map.Entry<Interval, T>> {

        private String mSequence = null;
        private Iterator<String> mSequenceIterator = null;
        private Iterator<IntervalTree.Node<T>> mTreeIterator = null;

        EntryIterator() {
            mSequenceIterator = mSequenceMap.keySet().iterator();
            advanceSequence();
        }

        public boolean hasNext() {
            return (mTreeIterator != null && mTreeIterator.hasNext());
        }

        public Map.Entry<Interval,T> next() {
            if (!hasNext()) {
                throw new NoSuchElementException("Iterator exhausted");
            }
            final IntervalTree.Node<T> node = mTreeIterator.next();
            if (!mTreeIterator.hasNext()) {
                advanceSequence();
            }
            final Interval key = new Interval(mSequence, node.getStart(), node.getEnd());
            final T value = node.getValue();
            return new MapEntry(key, value);
        }

        public void remove() {
            if (mTreeIterator == null) {
                throw new IllegalStateException("Iterator.next() has not been called");
            }
            mTreeIterator.remove();
        }

        private void advanceSequence() {
            while (mSequenceIterator.hasNext()) {
                mSequence = mSequenceIterator.next();
                mTreeIterator = mSequenceMap.get(mSequence).iterator();
                if (mTreeIterator.hasNext()) {
                    break;
                }
            }
        }
    }

    private class MapEntry
        implements Map.Entry<Interval,T> {

        private final Interval mKey;
        private T mValue;

        MapEntry(final Interval key, final T value) {
            mKey = key;
            mValue = value;
        }

        public Interval getKey() {
            return mKey;
        }

        public T getValue() {
            return mValue;
        }

        public T setValue(final T value) {
            mValue = value;
            return IntervalTreeMap.this.put(mKey, mValue);
        }
    }
}
