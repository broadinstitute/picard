/*
 * Copyright (c) 2007-2009 by The Broad Institute, Inc. and the Massachusetts Institute of Technology.
 * All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which
 * is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR WARRANTIES OF
 * ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT
 * OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR
 * RESPECTIVE TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES OF
 * ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES, ECONOMIC
 * DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER THE BROAD OR MIT SHALL
 * BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
 */
package org.broad.tribble.readers;

import org.broad.tribble.TribbleException;

import java.io.*;

/**
 * A simple class that provides readLine() functionality around a PositionalBufferedStream
 *
 * @author jrobinso
 */
public class AsciiLineReader implements LineReader {
    private static final int BUFFER_OVERFLOW_INCREASE_FACTOR = 2;
    private static final byte LINEFEED = (byte) ('\n' & 0xff);
    private static final byte CARRIAGE_RETURN = (byte) ('\r' & 0xff);

    PositionalBufferedStream is;
    char[] lineBuffer;

    public AsciiLineReader() {
        this(null);
    }

    public AsciiLineReader(InputStream is){
        this(new PositionalBufferedStream(is));
    }

    public AsciiLineReader(PositionalBufferedStream is) {
        this.is = is;
        // Allocate this only once, even though it is essentially a local variable of
        // readLine.  This makes a huge difference in performance
        lineBuffer = new char[10000];
    }

    public long getPosition(){
        if(is == null){
            throw new TribbleException("getPosition() called but no default stream was provided to the class on creation");
        }
        return is.getPosition();
    }

    /**
     * Read a line of text.  A line is considered to be terminated by any one
     * of a line feed ('\n'), a carriage return ('\r'), or a carriage return
     * followed immediately by a linefeed.
     *
     * @param stream the stream to read the next line from
     * @return A String containing the contents of the line or null if the
     *         end of the stream has been reached
     */
    public final String readLine(final PositionalBufferedStream stream) throws IOException {
        int linePosition = 0;

        while (true) {
            final int b = stream.read();

            if (b == -1) {
                // eof reached.  Return the last line, or null if this is a new line
                if (linePosition > 0) {
                    return new String(lineBuffer, 0, linePosition);
                } else {
                    return null;
                }
            }

            final char c = (char) (b & 0xFF);
            if (c == LINEFEED || c == CARRIAGE_RETURN) {
                if (c == CARRIAGE_RETURN && stream.peek() == LINEFEED) {
                    stream.read(); // <= skip the trailing \n in case of \r\n termination
                }

                return new String(lineBuffer, 0, linePosition);
            } else {
                // Expand line buffer size if neccessary.  Reserve at least 2 characters
                // for potential line-terminators in return string

                if (linePosition > (lineBuffer.length - 3)) {
                    char[] temp = new char[BUFFER_OVERFLOW_INCREASE_FACTOR * lineBuffer.length];
                    System.arraycopy(lineBuffer, 0, temp, 0, lineBuffer.length);
                    lineBuffer = temp;
                }

                lineBuffer[linePosition++] = c;
            }
        }
    }

    /**
     * Same as readLine(stream) but uses the stream provided in the constructure
     *
     * @return
     * @throws IOException
     */
    public final String readLine() throws IOException {
        if ( is == null ){
            throw new TribbleException("readLine() called without an explicit stream argument but no default stream was provided to the class on creation");
        }
        return readLine(is);
    }

    @Override
    public void close() {
        if ( is != null ) is.close();
        lineBuffer = null;
    }

    public static void main(String[] args) throws Exception {
        File testFile = new File(args[0]);
        final int iterations = Integer.valueOf(args[1]);
        final boolean includeBufferedReader = Boolean.valueOf(args[2]);
        long t0, lineCount, dt;
        double rate;

        System.out.printf("Testing %s%n", args[0]);
        for (int i = 0; i < iterations; i++) {
            if ( includeBufferedReader ) {
                BufferedReader reader2 = new BufferedReader(new FileReader(testFile));
                t0 = System.currentTimeMillis();
                lineCount = 0;
                while (reader2.readLine() != null) {
                    lineCount++;
                }
                dt = System.currentTimeMillis() - t0;
                rate = ((double) lineCount) / dt;
                printStatus("BufferedReader", lineCount, rate, dt);
                reader2.close();
            }

            PositionalBufferedStream pbs = new PositionalBufferedStream(new FileInputStream(testFile));
            LineReader reader = new AsciiLineReader(pbs);
            t0 = System.currentTimeMillis();
            lineCount = 0;
            while (reader.readLine() != null) {
                lineCount++;
            }
            dt = System.currentTimeMillis() - t0;
            rate = ((double) lineCount) / dt;
            printStatus("PositionalBufferedStream", lineCount, rate, dt);
            pbs.close();
        }
    }

    private static final void printStatus(final String name, long lineCount, double rate, long dt) {
        System.out.printf("%30s: %d lines read.  Rate = %.2e lines per second.  DT = %d%n", name, lineCount, rate, dt);
        System.out.flush();
    }
}

