/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.variant.bcf2;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Arrays;

/**
 * Simple holder for BCF version information
 *
 * User: depristo
 * Date: 8/2/12
 * Time: 2:16 PM
 */
public class BCFVersion {
    /**
     * BCF2 begins with the MAGIC info BCF_M_m where M is the major version (currently 2)
     * and m is the minor version, currently 1
     */
    public static final byte[] MAGIC_HEADER_START = "BCF".getBytes();

    final int majorVersion;
    final int minorVersion;

    public BCFVersion(int majorVersion, int minorVersion) {
        this.majorVersion = majorVersion;
        this.minorVersion = minorVersion;
    }

    /**
     * @return the major version number of this BCF file
     */
    public int getMajorVersion() {
        return majorVersion;
    }

    /**
     * @return the minor version number of this BCF file
     */
    public int getMinorVersion() {
        return minorVersion;
    }

    /**
     * Return a new BCFVersion object describing the major and minor version of the BCF file in stream
     *
     * Note that stream must be at the very start of the file.
     *
     * @param stream
     * @return a BCFVersion object, or null if stream doesn't contain a BCF file
     * @throws IOException
     */
    public static BCFVersion readBCFVersion(final InputStream stream) throws IOException {
        final byte[] magicBytes = new byte[MAGIC_HEADER_START.length];
        stream.read(magicBytes);
        if ( Arrays.equals(magicBytes, MAGIC_HEADER_START) ) {
            // we're a BCF file
            final int majorByte = stream.read();
            final int minorByte = stream.read();
            return new BCFVersion( majorByte, minorByte );
        } else
            return null;
    }

    /**
     * Write out the BCF magic information indicating this is a BCF file with corresponding major and minor versions
     * @param out
     * @throws IOException
     */
    public void write(final OutputStream out) throws IOException {
        out.write(MAGIC_HEADER_START);
        out.write(getMajorVersion() & 0xFF);
        out.write(getMinorVersion() & 0xFF);
    }

    @Override
    public String toString() {
        return String.format("BCF%d.%d", getMajorVersion(), getMinorVersion());
    }
}
