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

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * LRU collection class for managing objects that place some resource burden such that not too many of them
 * can existing in the VM at one time, but they can be reconstructed ias necessary.
 * Original motivation was for an LRU cache of FileOutputStreams.
 *
 * @author alecw@broadinstitute.org
 */
public class ResourceLimitedMap<Key, Value> {
    private static final float   hashTableLoadFactor = 0.75f;

    private final LinkedHashMap<Key,Value> map;
    private final int cacheSize;
    private final ResourceLimitedMapFunctor<Key, Value> functor;

    /**
     * Create LRU cache
     * @param cacheSize Max number of objects to be stored in the cache.
     * @param functor Encapsulates methods for creating a new object if it isn't in the cache, and
     * for finalizing an object that is getting LRU'ed out of the cache.
     */
    public ResourceLimitedMap(final int cacheSize, final ResourceLimitedMapFunctor<Key, Value> functor) {
        this.cacheSize = cacheSize;
        this.functor = functor;
        // Make hash table big enough so that it never gets resized.
        final int hashTableCapacity = (int)Math.ceil(cacheSize / hashTableLoadFactor) + 1;

        // Created LinkedHashMap in LRU mode
        map = new LinkedHashMap<Key,Value>(hashTableCapacity, hashTableLoadFactor, true) {
            @Override protected boolean removeEldestEntry (final Map.Entry<Key,Value> eldest) {
               if (size() > ResourceLimitedMap.this.cacheSize) {
                   ResourceLimitedMap.this.functor.finalizeValue(eldest.getKey(), eldest.getValue());
                   return true;
               } else {
                   return false;
               }
            }
        };
    }

    /**
     * Return an existing value, or create a new one if necessary.  If creating a new one and the
     * cache is full, the eldest object is pushed out of the cache.
     * @param key Key of desired value.
     * @return Either existing value, or new value created from key and inserted into map.
     */
    public Value get(final Key key) {
        if (!map.containsKey(key)) {
            map.put(key, functor.makeValue(key));
        }
        return map.get(key);
    }

    public Value remove(final Key key) {
        return map.remove(key);
    }

    /**
     * Determine if the map contains the given key.  Note that even if the map does not contain
     * the key, get(key) will return a value, because one will be created.
     * @param key
     * @return true if the map currently contains the given key.  It is unknown whether the map
     * may have contained the key in the past.
     */
    public boolean containsKey(final Key key) {
        return map.containsKey(key);
    }

    /**
     * Remove all the values from the map, and call functory.finalizeValue() on each of them.
     */
    public void finalizeAll() {
        for (final Map.Entry<Key, Value> entry : map.entrySet()) {
            functor.finalizeValue(entry.getKey(), entry.getValue());
        }
        map.clear();
    }
}
