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
package net.sf.samtools.util;

import net.sf.samtools.Defaults;

/**
 * Constants shared by BlockCompressed{Input,Output}Stream classes
 */
public class BlockCompressedStreamConstants {
    // Number of bytes in the gzip block before the deflated data.
    // This is not the standard header size, because we include one optional subfield,
    // but it is the standard for us.
    public static final int BLOCK_HEADER_LENGTH = 18;

    // Location in the gzip block of the total block size (actually total block size - 1)
    public static final int BLOCK_LENGTH_OFFSET = 16;

    // Number of bytes that follow the deflated data
    public static final int BLOCK_FOOTER_LENGTH = 8;

    // We require that a compressed block (including header and footer, be <= this)
    public static final int MAX_COMPRESSED_BLOCK_SIZE = 64 * 1024;

    // Gzip overhead is the header, the footer, and the block size (encoded as a short).
    public static final int GZIP_OVERHEAD = BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH + 2;

    // If Deflater has compression level == NO_COMPRESSION, 10 bytes of overhead (determined experimentally).
    public static final int NO_COMPRESSION_OVERHEAD = 10;

    // Push out a gzip block when this many uncompressed bytes have been accumulated.
    // This size is selected so that if data is not compressible,  if Deflater is given
    // compression level == NO_COMPRESSION, compressed size is guaranteed to be <= MAX_COMPRESSED_BLOCK_SIZE.
    public static final int DEFAULT_UNCOMPRESSED_BLOCK_SIZE = 64 * 1024 - (GZIP_OVERHEAD + NO_COMPRESSION_OVERHEAD);

    // Magic numbers
    public static final byte GZIP_ID1 = 31;
    public static final int GZIP_ID2 = 139;

    // FEXTRA flag means there are optional fields
    public static final int GZIP_FLG = 4;

    // extra flags
    public static final int GZIP_XFL = 0;

    // length of extra subfield
    public static final short GZIP_XLEN = 6;

    // The deflate compression, which is customarily used by gzip
    public static final byte GZIP_CM_DEFLATE = 8;

    public static final int DEFAULT_COMPRESSION_LEVEL = Defaults.COMPRESSION_LEVEL;

    // We don't care about OS because we're not doing line terminator translation
    public static final int GZIP_OS_UNKNOWN = 255;

    // The subfield ID
    public static final byte BGZF_ID1 = 66;
    public static final byte BGZF_ID2 = 67;

    // subfield length in bytes
    public static final byte BGZF_LEN = 2;

    public static final byte[] EMPTY_GZIP_BLOCK = {
            BlockCompressedStreamConstants.GZIP_ID1,
            (byte)BlockCompressedStreamConstants.GZIP_ID2,
            BlockCompressedStreamConstants.GZIP_CM_DEFLATE,
            BlockCompressedStreamConstants.GZIP_FLG,
            0, 0, 0, 0, // Modification time
            BlockCompressedStreamConstants.GZIP_XFL,
            (byte)BlockCompressedStreamConstants.GZIP_OS_UNKNOWN,
            BlockCompressedStreamConstants.GZIP_XLEN, 0, // Little-endian short
            BlockCompressedStreamConstants.BGZF_ID1,
            BlockCompressedStreamConstants.BGZF_ID2,
            BlockCompressedStreamConstants.BGZF_LEN, 0, // Little-endian short
            // Total block size - 1
            BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH +
                    BlockCompressedStreamConstants.BLOCK_FOOTER_LENGTH - 1 + 2, 0, // Little-endian short
            // Dummy payload?
            3, 0,
            0, 0, 0, 0, // crc
            0, 0, 0, 0, // uncompressedSize
    };
    public static final byte[] GZIP_BLOCK_PREAMBLE = {
            BlockCompressedStreamConstants.GZIP_ID1,
            (byte)BlockCompressedStreamConstants.GZIP_ID2,
            BlockCompressedStreamConstants.GZIP_CM_DEFLATE,
            BlockCompressedStreamConstants.GZIP_FLG,
            0, 0, 0, 0, // Modification time
            BlockCompressedStreamConstants.GZIP_XFL,
            (byte)BlockCompressedStreamConstants.GZIP_OS_UNKNOWN,
            BlockCompressedStreamConstants.GZIP_XLEN, 0, // Little-endian short
            BlockCompressedStreamConstants.BGZF_ID1,
            BlockCompressedStreamConstants.BGZF_ID2,
            BlockCompressedStreamConstants.BGZF_LEN, 0, // Little-endian short
    };
}
