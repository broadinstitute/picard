package net.sf.samtools.util;

/**
 * A simple tuple class.
 *
 * @author mccowan
 */
public class Tuple<A, B> {
    public final A a;
    public final B b;

    public Tuple(final A a, final B b) {
        this.a = a;
        this.b = b;
    }
}
