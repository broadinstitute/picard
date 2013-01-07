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

import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.index.linear.LinearIndex;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.io.FileNotFoundException;


public class IndexToTable {

    /**
     * this class:
     *  1) checks to see that the feature file exists
     *  2) loads an index from disk, if one doesn't exist, it creates it and writes it to disk
     *  3) creates a FeatureSource
     *  4) iterates over the records, emitting a final tally for the number of features seen
     *
     * @param args a single parameter, the file name to load
     */
    public static void main(String[] args) {

        // check yourself before you wreck yourself - we require one arg, the input file
        if (args.length != 2)
            printUsage();

        //LinearIndex.enableAdaptiveIndexing = false;
        LinearIndex idx = (LinearIndex)IndexFactory.loadIndex(new File(args[0]).getAbsolutePath());
        try {
            idx.writeTable(new PrintStream(new FileOutputStream(new File(args[1]))));
        } catch ( FileNotFoundException e ) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * print usage information
     */
    public static void printUsage() {
        System.err.println("Usage: java -jar IndexToTable.jar index.file output.table");
        System.exit(1);
    }
}