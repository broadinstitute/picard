/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This is copyright (2007-2009) by the Broad Institute/Massachusetts Institute
 * of Technology.  It is licensed to You under the Gnu Public License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 *  the License.  You may obtain a copy of the License at
 *
 *    http://www.opensource.org/licenses/gpl-2.0.php
 *
 * This software is supplied without any warranty or guaranteed support
 * whatsoever. Neither the Broad Institute nor MIT can be responsible for its
 * use, misuse, or functionality.
 */
package net.sf.samtools.seekablestream;

import java.io.IOException;
import java.io.InputStream;

public abstract class SeekableStream extends InputStream {

    public abstract long length();

    public abstract long position() throws IOException;

    public abstract void seek(long position) throws IOException;

    public abstract int read(byte[] buffer, int offset, int length) throws IOException;

    public abstract void close() throws IOException;

    public abstract boolean eof() throws IOException;

    /**
     * @return String representation of source (e.g. URL, file path, etc.), or null if not available.
     */
    public abstract String getSource();

    /**
     * Not sure why this is here, or what uses it
     */
//    public void readFully(byte b[]) throws IOException {
//        int len = b.length;
//        if (len < 0)
//            throw new IndexOutOfBoundsException();
//        int n = 0;
//        while (n < len) {
//            int count = read(b, n, len - n);
//            if (count < 0)
//                throw new EOFException();
//            n += count;
//        }
//    }

}
