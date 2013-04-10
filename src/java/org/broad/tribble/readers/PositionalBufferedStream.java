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
 * A wrapper around an {@code InputStream} which performs it's own buffering, and keeps track
 * of the position
 * @author Mark DePristo
 */
public final class PositionalBufferedStream extends InputStream implements Positional {
    final InputStream is;
    final byte[] buffer;
    int nextChar;
    int nChars;
    long position;

    public PositionalBufferedStream(final InputStream is) {
        this(is, 512000);
    }

    public PositionalBufferedStream(final InputStream is, final int bufferSize) {
        this.is = is;
        buffer = new byte[bufferSize];
        nextChar = nChars = 0;
    }

    public final long getPosition() {
        return position;
    }

    @Override
    public final int read() throws IOException {
        final int c = peek();
        if ( c >= 0 ) {
            // update position and buffer offset if peek says we aren't yet done
            position++;
            nextChar++;
        }
        return c;
    }

    @Override
    public final int read(final byte[] bytes, final int start, final int len) throws IOException {
        if ( len == 0 ) // If len is zero, then no bytes are read and 0 is returned
            return 0;
        if (nChars < 0) // If no byte is available because the stream is at end of file, the value -1 is returned
            return -1;
        else {
            int nRead = 0;
            int remaining = len;

            while ( remaining > 0 ) {
                // Try to Refill buffer if at the end of current buffer
                if ( nChars == nextChar )
                    if ( fill() < 0 )
                        break;

                // we copy as many bytes from the buffer as possible, up to the number of need
                final int nCharsToCopy = Math.min(nChars - nextChar, remaining);
                System.arraycopy(buffer, nextChar, bytes, start + nRead, nCharsToCopy);

                // update nextChar (pointer into buffer) and keep track of nRead and remaining
                nextChar  += nCharsToCopy;
                nRead     += nCharsToCopy;
                remaining -= nCharsToCopy;
            }

            // make sure we update our position tracker to reflect having advanced by nRead bytes
            position += nRead;
            return nRead;
        }
    }

    @Override
    public final int read(final byte[] bytes) throws IOException {
        return read(bytes, 0, bytes.length);
    }

    @Override
    public final boolean isDone() throws IOException {
        return nChars == -1 || peek() == -1;
    }

    @Override
    public final int peek() throws IOException {
        // Check for EOF
        if (nChars < 0) {
            return -1;
        } else if (nextChar == nChars){
            //Try to Refill buffer if at the end of current buffer
            if ( fill() < 0 ){
                return -1;
            }
        }

        return byteToInt(buffer[nextChar]);
    }

    private final int fill() throws IOException {
        nChars = is.read(buffer);
        nextChar = 0;
        return nChars;
    }

    public final long skip(final long nBytes) throws IOException {
        long remainingToSkip = nBytes;

        // because we have this buffer, that may be shorter than nBytes
        // we loop while there are bytes to skip, filling the buffer
        // When the buffer contains enough data that we have less than
        // its less left to skip we increase nextChar by the remaining
        // amount
        while ( remainingToSkip > 0 && ! isDone() ) {
            final long bytesLeftInBuffer = nChars - nextChar;
            if ( remainingToSkip > bytesLeftInBuffer ) {
                // we need to refill the buffer and continue our skipping
                remainingToSkip -= bytesLeftInBuffer;
                fill();
            } else {
                // there are enough bytes in the buffer to not read again
                // we just push forward the pointer nextChar
                nextChar += remainingToSkip;
                remainingToSkip = 0;
            }
        }

        final long actuallySkipped = nBytes - remainingToSkip;
        position += actuallySkipped;
        return actuallySkipped;
    }

    public final void close() {
        try {
            is.close();
        } catch (IOException ex) {
            new TribbleException("Failed to close PositionalBufferedStream", ex);
        }
    }

    private final static int byteToInt(byte b) {
        return b & 0xFF;
    }

    public static void main(String[] args) throws Exception {
        final File testFile = new File(args[0]);
        final int iterations = Integer.valueOf(args[1]);
        final boolean includeInputStream = Boolean.valueOf(args[2]);
        final boolean doReadFileInChunks = Boolean.valueOf(args[3]);

        System.out.printf("Testing %s%n", args[0]);
        for (int i = 0; i < iterations; i++) {
            if ( includeInputStream ) {
                final InputStream is = new FileInputStream(testFile);
                if ( doReadFileInChunks )
                    readFileInChunks("InputStream", is);
                else
                    readFileByLine("InputStream", is);
                is.close();
            }

            final PositionalBufferedStream pbs = new PositionalBufferedStream(new FileInputStream(testFile));
            if ( doReadFileInChunks )
                readFileInChunks("PositionalBufferedStream", pbs);
            else
                readFileByLine("PositionalBufferedStream", pbs);
            pbs.close();
        }
    }

    private static void readFileByLine(final String name, final InputStream is) throws IOException {
        final BufferedReader reader2 = new BufferedReader(new InputStreamReader(is));
        final long t0 = System.currentTimeMillis();
        long lineCount = 0;
        while (reader2.readLine() != null) {
            lineCount++;
        }
        final long dt = System.currentTimeMillis() - t0;
        final double rate = ((double) lineCount) / dt;
        printStatus(name, lineCount, rate, dt);
        reader2.close();
    }

    private static void readFileInChunks(final String name, final InputStream is) throws IOException {
        final long t0 = System.currentTimeMillis();
        long chunk = 0;
        final byte[] bytes = new byte[4096];
        while (is.read(bytes) != -1) {
            chunk++;
        }
        final long dt = System.currentTimeMillis() - t0;
        final double rate = ((double) chunk) / dt;
        printStatus(name, chunk, rate, dt);
        is.close();
    }


    private static final void printStatus(final String name, long lineCount, double rate, long dt) {
        System.out.printf("%30s: %d lines read.  Rate = %.2e lines per second.  DT = %d%n", name, lineCount, rate, dt);
        System.out.flush();
    }
}

