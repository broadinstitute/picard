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
package net.sf.samtools;

/**
 * The operators that can appear in a cigar string, and information about their disk representations.
 */
public enum CigarOperator {
    /** Match or mismatch */
    M(true, true,   'M'),
    /** Insertion vs. the reference. */
    I(true, false,  'I'),
    /** Deletion vs. the reference. */
    D(false, true,  'D'),
    /** Skipped region from the reference. */
    N(false, true,  'N'),
    /** Soft clip. */
    S(true, false,  'S'),
    /** Hard clip. */
    H(false, false, 'H'),
    /** Padding. */
    P(false, false, 'P'),
    /** Matches the reference. */
    EQ(true, true,  '='),
    /** Mismatches the reference. */
    X(true, true,   'X')
    ;

    // Representation of CigarOperator in BAM file
    private static final byte OP_M = 0;
    private static final byte OP_I = 1;
    private static final byte OP_D = 2;
    private static final byte OP_N = 3;
    private static final byte OP_S = 4;
    private static final byte OP_H = 5;
    private static final byte OP_P = 6;
    private static final byte OP_EQ = 7;
    private static final byte OP_X = 8;

    private final boolean consumesReadBases;
    private final boolean consumesReferenceBases;
    private final byte character;
    private final String string;

    // Readable synonyms of the above enums
    public static final CigarOperator MATCH_OR_MISMATCH = M;
    public static final CigarOperator INSERTION = I;
    public static final CigarOperator DELETION = D;
    public static final CigarOperator SKIPPED_REGION = N;
    public static final CigarOperator SOFT_CLIP = S;
    public static final CigarOperator HARD_CLIP = H;
    public static final CigarOperator PADDING = P;

    /** Default constructor. */
    CigarOperator(boolean consumesReadBases, boolean consumesReferenceBases, char character) {
        this.consumesReadBases = consumesReadBases;
        this.consumesReferenceBases = consumesReferenceBases;
        this.character = (byte) character;
        this.string = new String(new char[] {character}).intern();
    }

    /** If true, represents that this cigar operator "consumes" bases from the read bases. */
    public boolean consumesReadBases() { return consumesReadBases; }

    /** If true, represents that this cigar operator "consumes" bases from the reference sequence. */
    public boolean consumesReferenceBases() { return consumesReferenceBases; }

    /**
     * @param b CIGAR operator in character form as appears in a text CIGAR string
     * @return CigarOperator enum value corresponding to the given character.
     */
    public static CigarOperator characterToEnum(final int b) {
        switch (b) {
        case 'M':
            return M;
        case 'I':
            return I;
        case 'D':
            return D;
        case 'N':
            return N;
        case 'S':
            return S;
        case 'H':
            return H;
        case 'P':
            return P;
        case '=':
            return EQ;
        case 'X':
            return X;
        default:
            throw new IllegalArgumentException("Unrecognized CigarOperator: " + b);
        }
    }

    /**
     * @param i CIGAR operator in binary form as appears in a BAMRecord.
     * @return CigarOperator enum value corresponding to the given int value.
     */
    public static CigarOperator binaryToEnum(final int i) {
        switch(i) {
            case OP_M:
                return M;
            case OP_I:
                return I;
            case OP_D:
                return D;
            case OP_N:
                return N;
            case OP_S:
                return S;
            case OP_H:
                return H;
            case OP_P:
                return P;
            case OP_EQ:
                return EQ;
            case OP_X:
                return X;
            default:
                throw new IllegalArgumentException("Unrecognized CigarOperator: " + i);
        }
    }

    /**
     *
     * @param e CigarOperator enum value.
     * @return CIGAR operator corresponding to the enum value in binary form as appears in a BAMRecord.
     */
    public static int enumToBinary(final CigarOperator e) {
        switch(e) {
            case M:
                return OP_M;
            case I:
                return OP_I;
            case D:
                return OP_D;
            case N:
                return OP_N;
            case S:
                return OP_S;
            case H:
                return OP_H;
            case P:
                return OP_P;
            case EQ:
                return OP_EQ;
            case X:
                return OP_X;
            default:
                throw new IllegalArgumentException("Unrecognized CigarOperator: " + e);
        }
    }

    /** Returns the character that should be used within a SAM file. */
    public static byte enumToCharacter(final CigarOperator e) {
        return e.character;
    }

    /** Returns the cigar operator as it would be seen in a SAM file. */
    @Override public String toString() {
        return this.string;
    }
}
