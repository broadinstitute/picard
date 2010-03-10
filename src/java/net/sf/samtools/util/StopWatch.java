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
package net.sf.samtools.util;

/**
 * Utility to help in performance testing.
 */
public class StopWatch {

    private long startTime = 0;
    private long elapsedTime = 0;
    private boolean running = false;


    public void start() {
        this.startTime = System.currentTimeMillis();
        this.running = true;
    }


    public void stop() {
        long stopTime = System.currentTimeMillis();
        elapsedTime += stopTime - startTime;
        this.running = false;
    }

    public void reset() {
        stop();
        elapsedTime = 0;
        startTime = 0;
    }

    /**
     * Returns the cumulative time between all the start() and stop() calls made to this object.  If the
     * StopWatch is currently running, also includes the time between the most recent start() and now.
     * @return elapsedTime in milliseconds
     */
    public long getElapsedTime() {
        final long currentElapsed;
        if (running) {
             currentElapsed = (System.currentTimeMillis() - startTime);
        } else {
            currentElapsed = 0;
        }
        return currentElapsed + elapsedTime;
    }


    /**
     * @return same as getElapsedTime(), but truncated to seconds.
     */
    public long getElapsedTimeSecs() {
        return getElapsedTime() / 1000;
    }
}
