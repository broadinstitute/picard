/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
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

package picard.arrays.illumina;

import picard.PicardException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/**
 * A class to parse the contents of an Illumina Infinium Normalization Manifest file
 *
 * An Illumina Infinium Normalization Manifest file contains a subset of the information contained in the Illumina
 * Manifest file in addition to the normalization ID which is needed for normalizating intensities in GtcToVcf
 *
 */
public class InfiniumNormalizationManifest {
    private long[] positions;
    private byte[] chromosomes;
    private int[] normIds;
    private Integer[] allNormIds;

    private static final Map<String, Byte> CHROM_TO_BYTE = new HashMap<>();

    static {
        // Note.  Illumina uses Chromosome 0 (position 0) to flag probes they have failed in the manifest.
        for (int i = 0; i <= 22; i++) {
            String chrom = Integer.toString(i);
            CHROM_TO_BYTE.put(chrom, new Byte(chrom));
        }
        CHROM_TO_BYTE.put("X", new Byte("23"));
        CHROM_TO_BYTE.put("XY", new Byte("23"));        // XY is Illumina notation for PAR, the pseudo-autosomal region
        CHROM_TO_BYTE.put("Y", new Byte("24"));
        CHROM_TO_BYTE.put("MT", new Byte("25"));
    }

    public InfiniumNormalizationManifest(File illuminaNormalizationManifest) {
        parse(illuminaNormalizationManifest);
    }

    private void parse(File illuminaNormalizationManifest) {
        try (LineNumberReader lnr = new LineNumberReader(new FileReader(illuminaNormalizationManifest))) {
            Set<Integer> allNormIdSet = new TreeSet<>();
            lnr.skip(Long.MAX_VALUE);
            int numberOfSnps = lnr.getLineNumber() - 1;

            normIds = new int[numberOfSnps];
            positions = new long[numberOfSnps];
            chromosomes = new byte[numberOfSnps];
            try (BufferedReader reader = new BufferedReader(new FileReader(illuminaNormalizationManifest))) {
                boolean headerRead = false;
                int count = 0;
                while (reader.ready()) {
                    if (!headerRead) {
                        reader.readLine();
                        headerRead = true;
                    }

                    String line = reader.readLine();
                    String[] tokens = line.split(",");

                    allNormIdSet.add(new Integer(tokens[8].trim()));
                    normIds[count] = new Integer(tokens[8].trim());
                    String chrom = tokens[2].trim();
                    Byte chromByte = CHROM_TO_BYTE.get(chrom);
                    chromosomes[count] = chromByte;
                    positions[count] = new Long(tokens[3].trim());
                    count++;
                }
                allNormIds = allNormIdSet.toArray(new Integer[allNormIdSet.size()]);
            }
        } catch (IOException e) {
            throw new PicardException("Error parsing Infinium normalization manifest", e);
        }
    }

    Integer[] getAllNormIds() {
        return allNormIds;
    }

    int[] getNormIds() {
        return normIds;
    }

    public byte[] getChromosomes() {
        return chromosomes;
    }

    public long[] getPositions() {
        return positions;
    }
}

