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
 */
package net.sf.samtools.apps;

import java.io.RandomAccessFile;
import java.io.File;
import java.io.IOException;

/**
 * @author alecw@broadinstitute.org
 */
public class TimeRandomAccessFile {
    public static void main(String[] args) throws Exception {
        RandomAccessFile raf = new RandomAccessFile(new File(args[0]), "r");
        byte[] buf = new byte[64 * 1024];
        long totalBytesRead = 0;
        int bytesRead;
        while ((bytesRead = readBytes(raf, buf, 0, buf.length)) > 0) {
            totalBytesRead += bytesRead;
        }
        System.out.println("Total bytes: " + totalBytesRead);
    }
    private static int readBytes(final RandomAccessFile file, final byte[] buffer, final int offset, final int length)
        throws IOException {
        int bytesRead = 0;
        while (bytesRead < length) {
            final int count = file.read(buffer, offset + bytesRead, length - bytesRead);
            if (count <= 0) {
                break;
            }
            bytesRead += count;
        }
        return bytesRead;
    }
}
