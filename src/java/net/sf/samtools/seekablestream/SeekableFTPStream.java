/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package net.sf.samtools.seekablestream;

import java.io.IOException;
import java.net.URL;

/**
 * Unfortunately the seekable stream classes exist for both Tribble and Picard, and we need both.  This class
 * is for use with Tribble and delegates all the work to a helper.
 *
 * @author jrobinso
 * @date Oct 27, 2010
 */
public class SeekableFTPStream extends SeekableStream {

    SeekableFTPStreamHelper helper;
    public SeekableFTPStream(URL url) throws IOException {
        this(url, null);
    }

    public SeekableFTPStream(URL url, UserPasswordInput userPasswordInput) throws IOException {
        helper = new SeekableFTPStreamHelper(url, userPasswordInput);
    }

    public void seek(long position) {
        helper.seek(position);
    }

    public long position() {
        return helper.position();
    }

    @Override
    public boolean eof() throws IOException {
        return helper.eof();
    }

    @Override
    public String getSource() {
        return null; //TODO
    }

    @Override
    public long length() {
        return helper.length();
    }


    @Override
    public long skip(long n) throws IOException {
        return helper.skip(n);
    }


    @Override
    public int read(byte[] buffer, int offset, int len) throws IOException {
        return helper.read(buffer, offset, len);
    }


    public void close() throws IOException {
        helper.close();
    }

    public int read() throws IOException {
        return helper.read();
    }

//    private static final String EXPECTED = "Apache Software Foundation";
    private static final String EXPECTED1 = "\u00cf\u00ac\u00c9\u0075\u0043\u00d4\u00d5\u0079";
    private static final String EXPECTED2 = "\u00e4\u006c\u0077\u000c\u0016\u00f1\u0030\u008f";
    public static void main(String[] args) throws IOException {
//    	String testURL = (args.length < 1) ? "ftp://apache.cs.utah.edu/apache.org/HEADER.html" : args[0];
    	String testURL = (args.length < 1) ? "ftp://hgdownload.cse.ucsc.edu/goldenPath/panTro3/vsHg19/panTro3.hg19.all.chain.gz" : args[0];
        long startPosition = (args.length < 2) ? 0x0b66c78l : Long.parseLong(args[1]);
        int len = (args.length < 3) ? 8 : Integer.parseInt(args[2]);
        int skipLen = (args.length < 4) ? 0x18 : Integer.parseInt(args[3]);
        SeekableStream s = SeekableStreamFactory.getStreamFor(testURL);
        byte[] buffer = new byte[len];
        s.seek(startPosition);
        s.read(buffer, 0, len);
        if (s.position() != startPosition + len && s.position() != s.length()) {
        	System.out.println("1) updated position is incorrect");
        }
        String data = new String(buffer);
        System.out.println("1) read:" + data);
        if (args.length == 0) {
            System.out.println("1) expected:" + EXPECTED1);
            System.out.println("1) values do" + (EXPECTED1.equals(data) ? "" : " not") + " match");
        }
        s.skip(skipLen);
        s.read(buffer, 0, len);
        if (s.position() != startPosition + 2 * len + skipLen && s.position() != s.length()) {
        	System.out.println("2) updated position is incorrect");
        }
        String data2 = new String(buffer);
        System.out.println("2) read:" + data2);
        if (args.length == 0) {
            System.out.println("2) expected:" + EXPECTED2);
            System.out.println("2) values do" + (EXPECTED2.equals(data2) ? "" : " not") + " match");
        }
    }
}
