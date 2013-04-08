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

import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.Tribble;
import org.broad.tribble.bed.BEDCodec;
import org.broad.tribble.dbsnp.OldDbSNPCodec;
import org.broad.tribble.gelitext.GeliTextCodec;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.index.linear.LinearIndex;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Iterator;

/**
 * An example of how to index a feature file, and then count all the records in the file.
 * This is also useful for testing the feature reader
 */
public class CountRecords {

    /**
     *  The main method of this class:
     *  1) checks to see that the feature file exists
     *  2) loads an index from disk, if one doesn't exist, it creates it and writes it to disk
     *  3) creates a FeatureSource
     *  4) iterates over the records, emitting a final tally for the number of features seen
     *
     * @param args a single parameter, the file name to load
     */
    public static void main(String[] args) {
        
        // check yourself before you wreck yourself - we require one arg, the input file
        if (args.length > 2)
            printUsage();

        // our feature file
        File featureFile = new File(args[0]);
        if (!featureFile.exists()) {
            System.err.println("File " + featureFile.getAbsolutePath() + " doesnt' exist");
            printUsage();
        }

        int optimizeIndex = args.length == 2 ? Integer.valueOf(args[1]) : -1;

        // determine the codec
        FeatureCodec codec = getFeatureCodec(featureFile);

        runWithIndex(featureFile, codec, optimizeIndex);

    }

    /**
     *
     * @see org.broad.tribble.index.linear.LinearIndex#optimize(double)
     * @param featureInput  File containing features
     * @param codec  Codec used to read the features
     * @param optimizeThreshold Threshold used to optimize the linear index
     * @return
     */
    public static long runWithIndex(File featureInput, FeatureCodec codec, int optimizeThreshold) {
        // get an index
        Index index = loadIndex(featureInput, codec);
        if ( optimizeThreshold != -1 )
            ((LinearIndex)index).optimize(optimizeThreshold);

        // get a reader
        AbstractFeatureReader reader = null;
        try {
            reader =   AbstractFeatureReader.getFeatureReader(featureInput.getAbsolutePath(), codec, index);

            // now read iterate over the file
            long recordCount = 0l;

            // this call could be replaced with a query
            Iterator<Feature> iter = reader.iterator();

            // cycle through the iterators
            while (iter.hasNext()) {
                Feature feat = iter.next();
                ++recordCount;
            }

            System.err.println("We saw " + recordCount + " record in file " + featureInput);
            return recordCount;

        } catch (IOException e) {
            throw new RuntimeException("Something went wrong while reading feature file " + featureInput, e);
        }
    }

    /**
     * print usage information
     */
    public static void printUsage() {
        System.err.println("Usage: java -jar CountRecords.jar <inputFile>");
        System.err.println("    Where input can be of type: VCF (ends in .vcf or .VCF");
        System.err.println("                                Bed (ends in .bed or .bed");
        System.err.println("                                OldDbSNP (ends in .snp or .rod");
        /**
         * you could add others here; also look in the GATK code-base for an example of a dynamic way
         * to load Tribble codecs.
         */
        System.exit(1);
    }

    /**
     *
     * @param featureFile the feature file
     * @param codec the codec to decode features with
     * @return an index instance
     */
    public static Index loadIndex(File featureFile, FeatureCodec codec) {
        // lets setup a index file name
        File indexFile = Tribble.indexFile(featureFile);

        // our index instance;
        Index index = null;

        // can we read the index file
        if (indexFile.canRead()) {
            System.err.println("Loading index from disk for index file -> " + indexFile);
            index = IndexFactory.loadIndex(indexFile.getAbsolutePath());
        // else we want to make the index, and write it to disk if possible
        } else {
            System.err.println("Creating the index and memory, then writing to disk for index file -> " + indexFile);
            index = createAndWriteNewIndex(featureFile,indexFile,codec);
        }

        return index;
    }

    /**
     * creates a new index, given the feature file and the codec
     * @param featureFile the feature file (i.e. .vcf, .bed)
     * @param indexFile the index file; the location we should be writing the index to
     * @param codec the codec to read features with
     * @return an index instance
     */
    public static Index createAndWriteNewIndex(File featureFile, File indexFile, FeatureCodec codec) {
        try {
            Index index = IndexFactory.createLinearIndex(featureFile, codec);

            // try to write it to disk
            LittleEndianOutputStream stream = new LittleEndianOutputStream(new FileOutputStream(indexFile));
            index.write(stream);
            stream.close();

            return index;
        } catch (IOException e) {
            throw new RuntimeException("Unable to create index from file " + featureFile,e);
        }
    }

    /**
     * Return a {@code FeatureCodec} instance appropriate for the given
     * {@code featureFile}. Codec is generated based on file extension
     * @param featureFile
     * @return
     * @throws IllegalArgumentException If a codec cannot be found
     */
    public static FeatureCodec getFeatureCodec(File featureFile) {
        // quickly determine the codec type
        //if (featureFile.getName().endsWith(".vcf") || featureFile.getName().endsWith(".VCF") )
        //    return new VCFCodec();
        if (featureFile.getName().endsWith(".bed") || featureFile.getName().endsWith(".BED") )
            return new BEDCodec();
        if (featureFile.getName().endsWith(".snp") || featureFile.getName().endsWith(".rod") )
            return new OldDbSNPCodec();
        if (featureFile.getName().endsWith(".geli.calls") || featureFile.getName().endsWith(".geli") )
            return new GeliTextCodec();
        //if (featureFile.getName().endsWith(".txt") || featureFile.getName().endsWith(".TXT") )
        //    return new SoapSNPCodec();
        throw new IllegalArgumentException("Unable to determine correct file type based on the file name, for file -> " + featureFile);
    }
}
