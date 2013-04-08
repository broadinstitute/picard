/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broad.tribble.example;

import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;

public class ProfileIndexReading {

    /**
     * This class times the loading of an index file a number of times,
     * and prints the result of each trial
     * @param args Two parameters:
     *             1) The number of trials to run
     *             2) Index file to load
     */
    public static void main(String[] args) {

        // check yourself before you wreck yourself - we require one arg, the input file
        if (args.length < 2)
            printUsage();

        int iterations = Integer.valueOf(args[0]);
        for ( int j = 1; j < args.length; j++  ) {
            String indexFile = args[j];
            System.out.printf("Reading %s%n", indexFile);
            long startTime = System.currentTimeMillis();
            for ( int i = 0; i < iterations; i++ ) {
                System.out.printf("  iteration %d%n", i);
                Index index = IndexFactory.loadIndex(indexFile);
            }
            long stopTime = System.currentTimeMillis();
            System.out.printf("Runtime %s %.2f%n", indexFile, (stopTime - startTime) / 1000.0);
        }
    }

    /**
     * print usage information
     */
    public static void printUsage() {
        System.err.println("Usage: java -jar ReadIndices.jar iterations index.file...");
        System.exit(1);
    }
}