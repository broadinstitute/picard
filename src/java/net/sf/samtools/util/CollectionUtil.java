/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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
package net.sf.samtools.util;

import java.util.*;

/**
 * Small utility methods for dealing with collection classes.
 * @author mccowan
 */
public class CollectionUtil {

    public static <T> List<T> makeList (final T... list) {
        final List<T> result = new ArrayList<T>();
        Collections.addAll(result, list);

        return result;
    }
    
    public static <T> Set<T> makeSet (final T... list) {
        final Set<T> result = new HashSet<T>();
        Collections.addAll(result, list);
        return result;
    }
    
    public static <T> Collection<T> makeCollection (final Iterator<T> i) {
        final List<T> list = new LinkedList<T>();
        while (i.hasNext()) {
            list.add(i.next());
        }
        return list;
    }
    
    /** Construct a string by toString()ing each item in the collection with inBetween between each item. */
    public static String join(final Collection<?> items, final String inBetween) {
        final StringBuilder builder = new StringBuilder();
        for (final Object item : items) {
            if (builder.length() > 0) builder.append(inBetween);
            builder.append(item);
        }

        return builder.toString();
    }

    public static <T> T getSoleElement(final Collection<T> items) {
        if (items.size() != 1)
            throw new IllegalArgumentException(String.format("Expected a single element in %s, but found %s.", items, items.size()));
        return items.iterator().next();
    }
    
    /** Simple multi-map for convenience of storing collections in map values. */
    public static class MultiMap<K, V> extends HashMap<K, Collection<V>> {
        public void append(final K k, final V v) {
            this.initializeKeyIfUninitialized(k);
            this.get(k).add(v);
        }

        public void appendAll(final K k, final Collection<? extends V> v) {
            this.initializeKeyIfUninitialized(k);
            this.get(k).addAll(v);
        }

        private void initializeKeyIfUninitialized(final K k) {
            if (!this.containsKey(k))
                this.put(k, new LinkedList<V>());
        }
    }

    /** 
     * Partitions a collection into groups based on a characteristics of that group.  Partitions are embodied in a map, whose keys are the
     * value of that characteristic, and the values are the partition of elements whose characteristic evaluate to that key.
     */
    public static <K, V> Map<K,Collection<V>> partition(final Collection<V> collection,  final Partitioner<V, K> p) {
        final MultiMap<K, V> partitionToValues = new MultiMap<K, V>();
        for (final V entry : collection) {
            partitionToValues.append(p.getPartition(entry), entry);
        }
        return partitionToValues;
    }
    public static abstract class Partitioner<V, K> {
        public abstract K getPartition(final V v);
    }
    
    /**
     * A defaulting map, which returns a default value when a value that does not exist in the map is looked up.
     * 
     * This map supports two modes: injecting-on-default, and not injecting-on-default.  When injecting on default, when a lookup is
     * performed and a default value is returned, the default value is injected at that key, so that it now lives in the underlying map.
     * Without this mode, the value is simply returned and the underlying map is unaffected.
     * 
     * Note: When using injecting-on-default mode, and performing a lookup with a non-key type (the get method accepts any object), a 
     * class cast exception will be thrown because a non-key type cannot be added to the map.
     * @param <K>
     * @param <V>
     */
    public static class DefaultingMap<K, V> extends HashMap<K, V> {
        final Factory<V, K> defaultGenerator;
        final boolean injectValueOnDefault;
        
        /** Creates a defaulting map which defaults to the provided value and with injecting-on-default disabled. */
        public DefaultingMap(final V defaultValue) {
            this(new Factory<V, K>() {
                @Override
                public V make(final K k) {
                    return defaultValue;
                }
            }, false);
        }
        
        /**
         * Creates a defaulting map that generates defaults from the provided factory. This is useful when the default is non-static, or
         * the default is mutable, and the client wishes to get a value and mutate it and persist those changes in the map.
         */
        public DefaultingMap(final Factory<V, K> defaultGenerator, final boolean injectValueOnDefaulting) {
            this.defaultGenerator = defaultGenerator;
            this.injectValueOnDefault = injectValueOnDefaulting;
        }

        @Override
        @SuppressWarnings("unchecked") // Expect that the cast is successful; otherwise, client is breaking contract.
        public V get(final Object key) {
            if (!this.containsKey(key)) {
                final V val = this.defaultGenerator.make((K) key);
                if (this.injectValueOnDefault) {
                    this.put((K) key, val); 
                }
                return val;
            } else {
                return super.get(key);
            }
        }
        
        public interface Factory<V, K> {
            /**
             * @param k
             */
            V make(K k);
        }
    }

}
