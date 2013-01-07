/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
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
package org.broad.tribble.example;

import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;

import java.io.File;

/**
 * a quick example of how to index a feature file, and then count all the records in the file.  This is also useful
 * for testing the feature reader
 */
public class IndicesAreEqual {

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
        if ( args.length != 2 )
            printUsage();
        else {
            Index index1 = loadIndex(args[0]);
            Index index2 = loadIndex(args[1]);
            System.out.printf("%n");
            System.out.printf("index1: %s%n", args[0]);
            System.out.printf("index2: %s%n", args[1]);
            boolean eq = index1.equals(index2);
            System.out.printf("  equals() = %b%n", eq);
        }
    }

    /**
     * print usage information
     */
    public static void printUsage() {
        System.err.println("Usage: java -jar IndicesAreEqual.jar index1 index2");
        System.err.println("    Prints out true / false if index1 and index2 are equal");
        System.exit(1);
    }

    /**
     * @return an index instance
     */
    public static Index loadIndex(String filename) {
        //System.err.println("Loading index from disk for index file -> " + filename);
        File file = new File(filename);
        if (file.canRead()) {
            return IndexFactory.loadIndex(file.getAbsolutePath());
        } else {
            printUsage();
            return null;
        }
    }
}
