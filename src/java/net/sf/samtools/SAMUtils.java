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

import net.sf.samtools.util.StringUtil;

import java.util.ArrayList;
import java.util.List;


/**
 * Utilty methods.
 */
public final class SAMUtils
{
    // Representation of bases, one for when in low-order nybble, one for when in high-order nybble.
    private static final byte COMPRESSED_EQUAL_LOW = 0;
    private static final byte COMPRESSED_A_LOW = 1;
    private static final byte COMPRESSED_C_LOW = 2;
    private static final byte COMPRESSED_G_LOW = 4;
    private static final byte COMPRESSED_T_LOW = 8;
    private static final byte COMPRESSED_N_LOW = 15;
    private static final byte COMPRESSED_EQUAL_HIGH = COMPRESSED_EQUAL_LOW << 4;
    private static final byte COMPRESSED_A_HIGH = COMPRESSED_A_LOW << 4;
    private static final byte COMPRESSED_C_HIGH = COMPRESSED_C_LOW << 4;
    private static final byte COMPRESSED_G_HIGH = COMPRESSED_G_LOW << 4;
    private static final byte COMPRESSED_T_HIGH = (byte)(COMPRESSED_T_LOW << 4);
    private static final byte COMPRESSED_N_HIGH = (byte)(COMPRESSED_N_LOW << 4);
    public static final int MAX_PHRED_SCORE = 93;

    /**
     * Convert from a byte array containing =AaCcGgTtNn represented as ASCII, to a byte array half as long,
     * with =, A, C, G, T converted to 0, 1, 2, 4, 8, 15.
     * @param readBases Bases as ASCII bytes.
     * @return New byte array with bases represented as nybbles, in BAM binary format.
     */
    static byte[] bytesToCompressedBases(final byte[] readBases) {
        final byte[] compressedBases = new byte[(readBases.length + 1)/2];
        int i;
        for (i = 1; i < readBases.length; i+=2) {
            compressedBases[i/2] = (byte)(charToCompressedBaseHigh(readBases[i-1]) |
                                    charToCompressedBaseLow(readBases[i]));
        }
        // Last nybble
        if (i == readBases.length) {
            compressedBases[i/2] = charToCompressedBaseHigh((char)readBases[i-1]);
        }
        return compressedBases;
    }

    /**
     * Convert from a byte array with basese stored in nybbles, with =, A, C, G, T represented as 0, 1, 2, 4, 8, 15,
     * to a a byte array containing =AaCcGgTtNn represented as ASCII.
     * @param length Number of bases (not bytes) to convert.
     * @param compressedBases Bases represented as nybbles, in BAM binary format.
     * @param compressedOffset Byte offset in compressedBases to start.
     * @return New byte array with bases as ASCII bytes.
     */
    public static byte[] compressedBasesToBytes(final int length, final byte[] compressedBases, final int compressedOffset) {
        final byte[] ret = new byte[length];
        int i;
        for (i = 1; i < length; i+=2) {
            final int compressedIndex = i / 2 + compressedOffset;
            ret[i-1] = compressedBaseToByteHigh(compressedBases[compressedIndex]);
            ret[i] = compressedBaseToByteLow(compressedBases[compressedIndex]);
        }
        // Last nybble
        if (i == length) {
            ret[i-1] = compressedBaseToByteHigh(compressedBases[i/2 + compressedOffset]);
        }
        return ret;
    }

    /**
     * Convert from ASCII byte to BAM nybble representation of a base in low-order nybble.
     * @param base One of =AaCcGgTtNn.
     * @return Low-order nybble-encoded equivalent.
     */
    private static byte charToCompressedBaseLow(final int base) {
        switch (base) {
            case '=':
                return COMPRESSED_EQUAL_LOW;
            case 'a':
            case 'A':
                return COMPRESSED_A_LOW;
            case 'c':
            case 'C':
                return COMPRESSED_C_LOW;
            case 'g':
            case 'G':
                return COMPRESSED_G_LOW;
            case 't':
            case 'T':
                return COMPRESSED_T_LOW;
            case 'n':
            case 'N':
            case '.':
                return COMPRESSED_N_LOW;
            default:
                throw new IllegalArgumentException("Bad  byte passed to charToCompressedBase: " + base);
        }
    }

    /**
     * Convert from ASCII byte to BAM nybble representation of a base in high-order nybble.
     * @param base One of =AaCcGgTtNn.
     * @return High-order nybble-encoded equivalent.
     */
    private static byte charToCompressedBaseHigh(final int base) {
        switch (base) {
            case '=':
                return COMPRESSED_EQUAL_HIGH;
            case 'a':
            case 'A':
                return COMPRESSED_A_HIGH;
            case 'c':
            case 'C':
                return COMPRESSED_C_HIGH;
            case 'g':
            case 'G':
                return COMPRESSED_G_HIGH;
            case 't':
            case 'T':
                return COMPRESSED_T_HIGH;
            case 'n':
            case 'N':
            case '.':
                return COMPRESSED_N_HIGH;
            default:
                throw new IllegalArgumentException("Bad  byte passed to charToCompressedBase: " + base);
        }
    }

    /**
     * Convert from BAM nybble representation of a base in low-order nybble to ASCII byte.
     * @param base One of COMPRESSED_*_LOW, a low-order nybble encoded base.
     * @return ASCII base, one of ACGTN=.
     */
    private static byte compressedBaseToByteLow(final int base) {
        switch (base & 0xf) {
            case COMPRESSED_EQUAL_LOW:
                return '=';
            case COMPRESSED_A_LOW:
                return 'A';
            case COMPRESSED_C_LOW:
                return 'C';
            case COMPRESSED_G_LOW:
                return 'G';
            case COMPRESSED_T_LOW:
                return 'T';
            case COMPRESSED_N_LOW:
                return 'N';
            default:
                throw new IllegalArgumentException("Bad  byte passed to charToCompressedBase: " + base);
        }
    }

    /**
     * Convert from BAM nybble representation of a base in high-order nybble to ASCII byte.
     * @param base One of COMPRESSED_*_HIGH, a high-order nybble encoded base.
     * @return ASCII base, one of ACGTN=.
     */
    private static byte compressedBaseToByteHigh(final int base) {
        switch ((byte)(base & 0xf0)) {
            case COMPRESSED_EQUAL_HIGH:
                return '=';
            case COMPRESSED_A_HIGH:
                return 'A';
            case COMPRESSED_C_HIGH:
                return 'C';
            case COMPRESSED_G_HIGH:
                return 'G';
            case COMPRESSED_T_HIGH:
                return 'T';
            case COMPRESSED_N_HIGH:
                return 'N';
            default:
                throw new IllegalArgumentException("Bad  byte passed to charToCompressedBase: " + base);
        }
    }

    /**
     * Convert bases in place into canonical form, upper case, and with no-call represented as N.
     * @param bases
     */
    static void normalizeBases(final byte[] bases) {
        for (int i = 0; i < bases.length; ++i) {
            bases[i] = StringUtil.toUpperCase(bases[i]);
            if (bases[i] == '.') {
                bases[i] = 'N';
            }
        }
    }

    /**
     * Convert an array of bytes, in which each byte is a binary phred quality score, to
     * printable ASCII representation of the quality scores, ala FASTQ format.
     *
     * Equivalent to phredToFastq(data, 0, data.length)
     *
     * @param data Array of bytes in which each byte is a binar phred score.
     * @return String with ASCII representation of those quality scores.
     */
    public static String phredToFastq(final byte[] data) {
        if (data == null) {
            return null;
        }
        return phredToFastq(data, 0, data.length);
    }

    /**
     * Convert an array of bytes, in which each byte is a binary phred quality score, to
     * printable ASCII representation of the quality scores, ala FASTQ format.
     * @param buffer Array of bytes in which each byte is a binar phred score.
     * @param offset Where in buffer to start conversion.
     * @param length How many bytes of buffer to convert.
     * @return String with ASCII representation of those quality scores.
     */
    static String phredToFastq(final byte[] buffer, final int offset, final int length) {
        final char[] chars = new char[length];
        for (int i = 0; i < length; i++) {
            chars[i] = phredToFastq(buffer[offset+i] & 0xFF);
        }
        return new String(chars);
    }

    /**
     * Convert a single binary phred score to printable ASCII representation, ala FASTQ format.
     * @param phredScore binary phred score.
     * @return Printable ASCII representation of phred score.
     */
    static char phredToFastq(final int phredScore) {
        if (phredScore < 0 || phredScore > MAX_PHRED_SCORE) {
            throw new IllegalArgumentException("Cannot encode phred score: " + phredScore);
        }
        return (char) (33 + phredScore);
    }

    /**
     * Convert a string with phred scores in printable ASCII FASTQ format to an array
     * of binary phred scores.
     * @param fastq Phred scores in FASTQ printable ASCII format.
     * @return byte array of binary phred scores in which each byte corresponds to a character in the input string.
     */
    static byte[] fastqToPhred(final String fastq) {
        if (fastq == null) {
            return null;
        }
        final int length = fastq.length();
        final byte[] scores = new byte[length];
        for (int i = 0; i < length; i++) {
            scores[i] = (byte) fastqToPhred(fastq.charAt(i));
        }
        return scores;
    }

    /**
     * Converts printable qualities in Sanger fastq format to binary phred scores.
     */
    public static void fastqToPhred(final byte[] fastq) {
        for (int i = 0; i < fastq.length; ++i) {
            fastq[i] = (byte)fastqToPhred((char)(fastq[i] & 0xff));
        }
    }

    /**
     * Convert a single printable ASCII FASTQ format phred score to binary phred score.
     * @param ch Printable ASCII FASTQ format phred score.
     * @return Binary phred score.
     */
    static int fastqToPhred(final char ch) {
        if (ch < 33 || ch > 126) {
            throw new IllegalArgumentException("Invalid fastq character: " + ch);
        }
        return (ch - 33);
    }

    /**
     * calculate the bin given an alignment in [beg,end)
     * Copied from SAM spec.
     * @param beg 0-based start of read (inclusive)
     * @param end 0-based end of read (exclusive)
     */
    static int reg2bin(final int beg, int end)
    {

        --end;

        if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
        if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
        if (beg>>20 == end>>20) return  ((1<<9)-1)/7 + (beg>>20);
        if (beg>>23 == end>>23) return  ((1<<6)-1)/7 + (beg>>23);
        if (beg>>26 == end>>26) return  ((1<<3)-1)/7 + (beg>>26);
        return 0;
    }

    /**
     * Handle a list of validation errors according to the validation stringency.
     * @param validationErrors List of errors to report, or null if there are no errors.
     * @param samRecordIndex Record number of the SAMRecord corresponding to the validation errors, or -1 if
     * the record number is not known.
     * @param validationStringency If STRICT, throw a RuntimeException.  If LENIENT, print the validation
     * errors to stderr.  If SILENT, do nothing.
     */
    static void processValidationErrors(final List<SAMValidationError> validationErrors,
                                        final long samRecordIndex,
                                        final SAMFileReader.ValidationStringency validationStringency) {
        if (validationErrors != null && validationErrors.size() > 0) {
            for (final SAMValidationError validationError : validationErrors) {
                validationError.setRecordNumber(samRecordIndex);
            }
            if (validationStringency == SAMFileReader.ValidationStringency.STRICT) {
                throw new RuntimeException("SAM validation error: " + validationErrors.get(0));
            }
            else if (validationStringency == SAMFileReader.ValidationStringency.LENIENT) {
                for (final SAMValidationError error : validationErrors) {
                    System.err.println("Ignoring SAM validation error: " + error);
                }
            }
        }
    }

    static void processValidationError(final SAMValidationError validationError,
                                       SAMFileReader.ValidationStringency validationStringency) {
        if (validationStringency == SAMFileReader.ValidationStringency.STRICT) {
            throw new RuntimeException("SAM validation error: " + validationError);
        }
        else if (validationStringency == SAMFileReader.ValidationStringency.LENIENT) {
            System.err.println("Ignoring SAM validation error: " + validationError);
        }
        
    }


    /**
     * Chains <code>program</code> in front of the first "head" item in the list of
     * SAMProgramRecords in <code>header</code>.  This method should not be used
     * when there are multiple chains of program groups in a header, only when
     * it can safely be assumed that there is only one chain.  It correctly handles
     * the case where <code>program</code> has already been added to the header, so
     * it can be used whether creating a SAMProgramRecord with a constructor or when
     * calling SAMFileHeader.createProgramRecord().
     */
    public static void chainSAMProgramRecord(SAMFileHeader header, SAMProgramRecord program) {

        List<SAMProgramRecord> pgs = header.getProgramRecords();
        if (pgs.size() > 0) {
            List<String> referencedIds = new ArrayList<String>();
            for (SAMProgramRecord pg : pgs) {
                if (pg.getPreviousProgramGroupId() != null) {
                    referencedIds.add(pg.getPreviousProgramGroupId());
                }
            }
            for (SAMProgramRecord pg : pgs) {
                // if record being chained has already been added, ignore it
                if (pg.getProgramGroupId().equals(program.getProgramGroupId())) {
                    continue;
                }
                if (!referencedIds.contains(pg.getProgramGroupId())) {
                    program.setPreviousProgramGroupId(pg.getProgramGroupId());
                    break;
                }
            }
        }
    }
}
