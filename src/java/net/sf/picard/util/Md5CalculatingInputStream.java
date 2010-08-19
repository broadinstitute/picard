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
 * THE SOFTWARE.
 */
package net.sf.picard.util;

import net.sf.picard.PicardException;

import java.io.*;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

/**
 * Class to generate an MD5 string for a file as it is being read
 *
 * @author ktibbett@broadinstitue.org
 */
public class Md5CalculatingInputStream extends InputStream {

    private final InputStream is;
    private final MessageDigest md5;
    private final File digestFile;

    /**
     * Constructor that takes in the InputStream that we are wrapping
     * and creates the MD5 MessageDigest
     */
    public Md5CalculatingInputStream(InputStream is, File digestFile) {
        super();
        this.is = is;
        this.digestFile = digestFile;
        try {
            md5 = MessageDigest.getInstance("MD5");
            md5.reset();
        }
        catch (NoSuchAlgorithmException e) {
            throw new PicardException("MD5 algorithm not found", e);
        }
    }

    public int read() throws IOException {
        int result = is.read();
        if (result != -1) md5.update((byte)result);
        return result;
    }

    public int read(byte[] b) throws IOException {
        int result = is.read(b);
        if (result != -1) md5.update(b, 0, result);
        return result;
    }


    public int read(byte[] b, int off, int len) throws IOException {
        int result = is.read(b, off, len);
        if (result != -1) md5.update(b, off, result);
        return result;
    }

    private String getDigest() {
        String hash = new BigInteger(1, md5.digest()).toString(16);
        if (hash.length() != 32) {
            final String zeros = "00000000000000000000000000000000";
            hash = zeros.substring(0, 32 - hash.length()) + hash;
        }
        return hash;
    }

    public void close() throws IOException {
        is.close();
        BufferedWriter writer = new BufferedWriter(new FileWriter(digestFile));
        writer.write(getDigest());
        writer.close();
    }

    // Methods not supported or overridden because they would not result in a valid hash
    public boolean markSupported() { return false; }
    public void	mark(int readlimit) {
        throw new UnsupportedOperationException("mark() is not supported by the MD5CalculatingInputStream");
    }
    public void	reset() throws IOException {
        throw new UnsupportedOperationException("reset() is not supported by the MD5CalculatingInputStream");
    }
    public long skip(long n) throws IOException {
        throw new UnsupportedOperationException("skip() is not supported by the MD5CalculatingInputStream");
    }

    // Methods delegated to the wrapped InputStream
    public int available() throws IOException { return is.available(); }

}
