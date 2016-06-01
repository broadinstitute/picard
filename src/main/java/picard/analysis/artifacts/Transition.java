package picard.analysis.artifacts;

import htsjdk.samtools.util.SequenceUtil;

import java.util.Arrays;

enum Transition {
    AtoA('A','A'), AtoC('A','C'), AtoG('A','G'), AtoT('A','T'),
    CtoA('C','A'), CtoC('C','C'), CtoG('C','G'), CtoT('C','T'),
    GtoA('G','A'), GtoC('G','C'), GtoG('G','G'), GtoT('G','T'),
    TtoA('T','A'), TtoC('T','C'), TtoG('T','G'), TtoT('T','T');

    private static final Transition[] ALT_VALUES = new Transition[]{
        AtoC, AtoG, AtoT, CtoA, CtoG, CtoT, GtoA, GtoC, GtoT, TtoA, TtoC, TtoG
    };

    private final char ref;
    private final char call;

    Transition(final char ref, final char call) {
        this.ref = ref;
        this.call = call;
    }

    public static Transition transitionOf(final char ref, final char call) {
        return transitionIndexMap[baseIndexMap[ref]][baseIndexMap[call]];
    }

    /**
     * Like values(), but ignores the ref:ref "transitions".
     */
    public static Transition[] altValues() { return ALT_VALUES; }

    /**
     * Return the complementary transition. Both ref and call must be complemented.
     */
    public Transition complement() {
        return transitionOf((char) SequenceUtil.complement((byte) this.ref), (char) SequenceUtil.complement((byte) this.call));
    }

    public char ref() { return this.ref; }

    public char call() { return this.call; }

    @Override
    public String toString() {
        return this.ref + ">" + this.call;
    }

    protected enum Base {
        A ('A'),
        C ('C'),
        G ('G'),
        T ('T');

        public byte base;

        private Base(final char base) {
            this.base = (byte)base;
        }
    }

    // only 4 of the values will be used but we want to optimize for speed by accessing via the int value of a char
    static protected final int[] baseIndexMap = new int[256];
    static {
        Arrays.fill(baseIndexMap, -1);
        for (final Base b : Base.values()) {
            baseIndexMap[b.base] = b.ordinal();
        }
    }

    static private final Transition[][] transitionIndexMap = new Transition[Base.values().length][Base.values().length];
    static {
        for (final Base b1 : Base.values()) {
            for (final Base b2 : Base.values()) {
                transitionIndexMap[b1.ordinal()][b2.ordinal()] = Transition.valueOf(b1.toString() + "to" + b2.toString());
            }
        }
    }
}