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
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;

class ZCallPedFile {
    private static final Log log = Log.getInstance(ZCallPedFile.class);
    private static final ProgressLogger logger = new ProgressLogger(log, 10000);

    private static final int OFFSET = 6;
    private final Map<String, String> snpToAlleleMap = new HashMap<>();

    private void addAllele(final String snp, final String allele) {
        logger.record("0", 0);
        snpToAlleleMap.put(snp, allele);
    }

    String getAlleles(final String snp) {
        return snpToAlleleMap.get(snp);
    }

    public static ZCallPedFile fromFile(final File pedFile,
                                        final File mapFile) throws FileNotFoundException {
        final String[] pedFileFields = IOUtil.slurp(pedFile).split(" ");
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
            final String alleles = pedFileFields[index] + pedFileFields[index + 1];
            zCallPedFile.addAllele(mapFileLines[i].split("\\s")[1], alleles);
        }
        return zCallPedFile;
    }
}
