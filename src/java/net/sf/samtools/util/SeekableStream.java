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
package net.sf.samtools.util;

import java.io.IOException;
import java.io.InputStream;

public abstract class SeekableStream extends InputStream {

    public abstract long length();

    public abstract void seek(long position) throws IOException;

    public abstract int read(byte[] buffer, int offset, int length) throws IOException;

    public abstract void close() throws IOException;

    public abstract boolean eof() throws IOException;

    /**
     * @return String representation of source (e.g. URL, file path, etc.), or null if not available.
     * Should end with .bam if not null.
     */
    public abstract String getSource();
}