package net.sf.samtools.util;

/**
 * Simple utility for building an on-demand (lazy) object-initializer.
 * 
 * Works by accepting an initializer describing how to build the on-demand object, which is only called once and only after the first
 * invocation of {@link #get()} (or it may not be called at all).
 * 
 * @author mccowan
 */
public class Lazy<T> {
    private final LazyInitializer<T> initializer;
    private boolean isInitialized = false;
    private T instance;

    /** Simple cons */
    public Lazy(final LazyInitializer<T> initializer) {
        this.initializer = initializer;
    }

    /** Returns the instance associated with this {@link Lazy}, initializing it if necessary. */
    public synchronized T get() {
        if (!isInitialized) {
            this.instance = initializer.make();
            isInitialized = true;
        }
        return instance;
    }

    /** Describes how to build the instance of the lazy object. */
    public interface LazyInitializer<T> {
        /** Returns the desired object instance. */
        T make();
    }

    public boolean isInitialized() {
        return isInitialized;
    }
}
