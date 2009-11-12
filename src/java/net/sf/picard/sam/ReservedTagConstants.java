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
package net.sf.picard.sam;

import net.sf.samtools.SAMTag;

/**
 * Constants for tags used in our SAM/BAM files
 */
public class ReservedTagConstants {
    public static final String READ_GROUP_ID = SAMTag.RG.name(); // Specified in the SAM spec doc
    public static final String PROGRAM_GROUP_ID =  SAMTag.PG.name(); // Specified in the SAM spec doc

    /** Present and set to 1 if a read is a noise read. */
    public static final String XN = "XN";

    /** Number of nucleotide differences (Specified in the SAM spec doc) */
    public static final String NM = SAMTag.NM.name();

    /** The sum of the mismatched qualities. */
    public static final String XQ = "XQ";

    /**
     * The name of an attribute which stores the 1-based index of the start of
     * sequence within a read (in original orientation) that should be clipped
     * or trimmed before alignment and downstream use.
     */
    public static final String XT = "XT";
}
