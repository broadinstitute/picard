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

package picard.arrays;

import htsjdk.samtools.util.IOUtil;
import picard.PicardException;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;

/**
 * Represents a single-sample .ped file output from ZCall.
 */
class ZCallPedFile {
    private static final int OFFSET = 6;
    private final Map<String, String> snpToAlleleMap = new HashMap<>();

    private void addAllele(final String snp, final String allele) {
        snpToAlleleMap.put(snp, allele);
    }

    String getAlleles(final String snp) {
        return snpToAlleleMap.get(snp);
    }

    /**
     *
     * @param pedFile  .ped files are whitespace-separated files.
     *                The .ped format is defined at http://zzz.bwh.harvard.edu/plink/data.shtml#ped.
     *                The first six columns, in the following order, are mandatory:
     *                  Family ID
     *                  Individual ID
     *                  Paternal ID
     *                  Maternal ID
     *                  Sex (1 = male, 2 = female, other = unknown)
     *                  Phenotype
     *                The seventh column onward should be biallelic genotype data. ZCall outputs these as A or B,
     *                representing which cluster an allele falls in. Each element of the biallelic pairs should still
     *                be tab-separated.
     *                This file should be a single line, representing a single sample.
     * @param mapFile .map files are whitespace-separated files.
     *                The .map format is defined at http://zzz.bwh.harvard.edu/plink/data.shtml#map.
     *                It has exactly four columns in the following order:
     *                  Chromosome (1-22, X, Y, or 0 if unknown)
     *                  rs# or SNP identifier
     *                  Genetic distance in morgans
     *                  Base-pair position in bp units
     * @return A ZCallPedFile representing the input .ped and .map files.
     * @throws FileNotFoundException
     */
    public static ZCallPedFile fromFile(final File pedFile,
                                        final File mapFile) throws FileNotFoundException {
        final String[] pedFileLines = IOUtil.slurpLines(pedFile).toArray(new String[0]);
        if (pedFileLines.length > 1) {
            throw new PicardException("Only single-sample .ped files are supported.");
        }
        final String[] pedFileFields = IOUtil.slurp(pedFile).split("\\s");
        final String[] mapFileLines = IOUtil.slurpLines(mapFile).toArray(new String[0]);

        final ZCallPedFile zCallPedFile = new ZCallPedFile();

        /* first six fields are ignored
            Family ID
            Individual ID
            Paternal ID
            Maternal ID
            Sex (1=male; 2=female; other=unknown)
            Phenotype
         */
        //two fields for each snp (each allele)
        for (int i = 0; i < mapFileLines.length; i++) {
            final int index = (i * 2) + OFFSET;
            // The fields are supposed to be one character each
            if (pedFileFields[index].length() != 1 || pedFileFields[index + 1].length() != 1) {
                throw new PicardException("Malformed file: each allele should be a single character.");
            }
            final String alleles = pedFileFields[index] + pedFileFields[index + 1];
            zCallPedFile.addAllele(mapFileLines[i].split("\\s")[1], alleles);
        }
        return zCallPedFile;
    }
}
