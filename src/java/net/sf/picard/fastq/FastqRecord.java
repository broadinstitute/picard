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
package net.sf.picard.fastq;

/**
 * Represents a fastq record, fairly literally, i.e. without any conversion.
 */
public class FastqRecord {

    private String seqHeaderPrefix;
    private String seqLine;
    private String qualHeaderPrefix;
    private String qualLine;

    public FastqRecord(final String seqHeaderPrefix, final String seqLine, final String qualHeaderPrefix, final String qualLine) {
        if (seqHeaderPrefix.length() > 0) this.seqHeaderPrefix = seqHeaderPrefix;
        if (qualHeaderPrefix.length() > 0) this.qualHeaderPrefix = qualHeaderPrefix;
        this.seqLine = seqLine ;
        this.qualLine = qualLine ;
    }

    public String getReadHeader() { return seqHeaderPrefix; }
    public String getReadString() { return seqLine; }
    public String getBaseQualityHeader() { return qualHeaderPrefix; }
    public String getBaseQualityString() { return qualLine; }
}
